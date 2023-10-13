/*
 * @Description: jacobian_test
 * @Author: Liu Chenxi
 * @Date: 2023-06-03 21:58:10
 * @LastEditTime: 2023-06-04 15:44:45
 */
#include <iostream>
#include "common/eigen_types.h"
#include "common/nav_state.h"
#include "ch3/static_imu_init.h"
#include "ch3/utm_convert.h"
#include "ch4/gins_pre_integ.h"


#include <gflags/gflags.h>
#include <glog/logging.h>
#include <fstream>
#include <iomanip>

using namespace sad;

Vec3d grav(0, 0, -9.8);      // Z 向上，重力方向为负
double dt = 0.05;            // dt_ij
double angular_rad = 20 * sad::math::kDEG2RAD;  // 弧度制角度
Vec3d omega(0, 0, angular_rad);   // 角度矢量

// Compute Error
// imu 9x1
//崔：根据i时刻和j时刻IMU的姿态：包括旋转和位置，计算IMU预积分残差 公式4.41
//崔：形参说明：i时刻位姿 j时刻位姿 i时刻速度 j时刻速度 IMU预计分量 Rij v_ij p_ij
Vec9d Eimu(SE3& posei, SE3& posej, Vec3d& vi, Vec3d& vj, SO3& dR, Vec3d& dv, Vec3d& dp) {
    const Vec3d er = (dR.inverse() * posei.so3().inverse() * posej.so3()).log();
    Mat3d RiT = posei.so3().inverse().matrix();
    const Vec3d ev = RiT * (vj - vi - grav * dt) - dv;
    const Vec3d ep = RiT * (posej.translation() - posei.translation() - vi * dt
                     - grav * dt * dt / 2) - dp;
    Vec9d error;
    //写成9维残差变量
    error << er, ev, ep;
    return error;
}
//GNSS和先验的残差计算就是简单的两个状态相减（应用时为传感器测量的状态和待优化的状态之间相减），旋转部分广义上的相减
// GNSS 6x1
Vec6d Egnss(SE3& posej, SE3& measurement) {
    Vec6d error;
    error.head<3>() = (measurement.so3().inverse() * posej.so3()).log();
    error.tail<3>() = posej.translation() - measurement.translation();
    return error;    
}
// prior 15x1
Eigen::Matrix<double, 15, 1> Eprior(SE3& posei, Vec3d& vi, Vec3d& bgi,Vec3d& bai, NavStated& state) {
    const Vec3d er = SO3(state.R_.matrix().transpose() * posei.so3().matrix()).log();
    const Vec3d ep = posei.translation() - state.p_;
    const Vec3d ev = vi - state.v_;
    const Vec3d ebg = bgi - state.bg_;
    const Vec3d eba = bai - state.ba_;

    Eigen::Matrix<double, 15, 1> error;
    error << er, ep, ev, ebg, eba;
    return error;
}


int main() {
    // parameters
    SE3 posei, posej;   // qp1  qp2
    Vec3d vi, vj;       // vi   vj
    Vec3d bai, bgi;     // bai  bgi

    // Random initial value
    posej.so3() = posej.so3() * SO3::exp(omega);
    posej.translation() << 0.5, 0.5, 0.5;
    vi << 1, 2, 3;
    vj << 2, 2, 2;
    bai << 0.02, 0.03, 0.04;
    bgi << 0.005, 0.005, 0.005;
    SO3 dR;
    Vec3d dv(0, 0, 0);
    Vec3d dp(0, 0, 0);

    Mat3d dR_dbg = Mat3d::Identity();
    Mat3d dv_dbg = Mat3d::Identity();
    Mat3d dp_dbg = Mat3d::Identity();
    Mat3d dv_dba = Mat3d::Identity();
    Mat3d dp_dba = Mat3d::Identity();

    // Small update value
    double delta = 1e-5;




    //FIXME: IMU  ******************************************************************** // 
    //       | R1 | p1 | v1 | bg1 | ba1 | R2 | p2 | v2 |
    // 1. Numerical jacobian
    // J1 = ( f(x+dx) -f(x-dx) ) / (2 * dx)
    //崔：IMU预积分残差的雅可比矩阵为9X24维，初始化为0 9维是IMU预积分残差的九个维度，旋转预积分残差，速度预积分残差，位置预积分残差 24维是i时刻IMU的旋转 位置 速度 角速度零偏 加速度零偏 j时刻的旋转 位置 速度
    Eigen::Matrix<double, 9, 24> J1 = Eigen::Matrix<double, 9, 24>::Zero();
    // posei
    //崔：设置i时刻IMU旋转和位置的更新量 均在SE3里表示，用于计算，但是数值求导时，R和p的更新都需要单独考虑
    //崔：先只更新旋转R
    for (int i = 0; i < 3; ++i) {
        //设置旋转微小更新量，用于进行数值求导
        Vec3d update_se3 = Vec3d::Zero();
        //更新量包含旋转的三个维度，设置为相同值
        update_se3(i) = delta;
        //计算更新后的旋转姿态，一个作广义加法，一个作广义减法
        SE3 SE3_updated1 = posei;
        SE3 SE3_updated2 = posei;
        //广义加法
        SE3_updated1.so3() = SE3_updated1.so3() * SO3::exp(update_se3);
        //广义减法
        SE3_updated2.so3() = SE3_updated2.so3() * SO3::exp(-update_se3);
        Vec9d err1 = Eimu(SE3_updated1, posej, vi, vj, dR, dv, dp);
        Vec9d err2 = Eimu(SE3_updated2, posej, vi, vj, dR, dv, dp);
        //将数值求导的结果存在对应的9X24维的前三列，为预积分对i时刻旋转的求导
        J1.col(i) = (err1 - err2) / (2 * delta);
    }
    //崔：只更新位置p 同样用se3表示，对应se3中的translation
    for (int i = 0; i < 3; ++i) {
        Vec3d update_se3 = Vec3d::Zero();
        update_se3(i) = delta;
        //计算更新后的位置
        SE3 SE3_updated1 = posei;
        SE3 SE3_updated2 = posei;
        //一个加，一个减     
        SE3_updated1.translation() += update_se3;
        SE3_updated2.translation() -= update_se3;
        Vec9d err1 = Eimu(SE3_updated1, posej, vi, vj, dR, dv, dp);
        Vec9d err2 = Eimu(SE3_updated2, posej, vi, vj, dR, dv, dp);
        //将数值求导的结果存在9X24维矩阵的第4-6列，为预积分对i时刻位置的求导
        J1.col(3+i) = (err1 - err2) / (2 * delta);
    }
    // vi
    // 数值求导方法计算预积分残差对i时刻速度的导数
    for (int i = 0; i < 3; ++i) {
        Vec3d update_v = Vec3d::Zero();
        //设置速度微小更新量
        update_v(i) = delta;
        //计算i时刻更新的速度，一个加 一个减
        Vec3d update_vi1 = vi + update_v;
        Vec3d update_vi2 = vi - update_v;
        //计算需要对应的两个imu残差
        Vec9d err1 = Eimu(posei, posej, update_vi1, vj, dR, dv, dp);
        Vec9d err2 = Eimu(posei, posej, update_vi2, vj, dR, dv, dp);
        //将数值求导结果存入9X24维矩阵的第7-9列，为预积分对i时刻速度的求导
        J1.col(6+i) = (err1 - err2) / (2 * delta);
    }
    // bgi
    // 数值求导方法计算预积分 残差对i时刻imu角速度零偏的导数
    for (int i = 0; i < 3; ++i) {
        Vec3d update_bg = Vec3d::Zero();
        //设置角速度零偏微小更新量
        update_bg(i) = delta;
        //当角速度零偏更新时，利用公式4.32修正预积分观测量，旋转 速度 位置预积分都要更新
        //其中的偏导数由公式4.38更新，其中右乘近似雅可比可根据十四讲4.33公式计算
        SO3 dR1 = dR * SO3::exp(dR_dbg * (update_bg));
        SO3 dR2 = dR * SO3::exp(dR_dbg * (-update_bg));
        Vec3d dv1 = dv + dv_dbg * (update_bg);
        Vec3d dv2 = dv + dv_dbg * (-update_bg);
        Vec3d dp1 = dp + dp_dbg * (update_bg);
        Vec3d dp2 = dp + dp_dbg * (-update_bg);
        //计算更新后的残差值
        Vec9d err1 = Eimu(posei, posej, vi, vj, dR1, dv1, dp1);
        Vec9d err2 = Eimu(posei, posej, vi, vj, dR2, dv2, dp2);
        //数值求导结果存入9X24维矩阵的10-12列
        J1.col(9+i) = (err1 - err2) / (2 * delta);
    }      
    // bai
    // 数值求导方法计算预积分 残差对i时刻imu加速度零偏的导数
    for (int i = 0; i < 3; ++i) {
        Vec3d update_ba = Vec3d::Zero();
        update_ba(i) = delta;
        //加速度零偏更新时，速度和位置预积分根据4.32更新
        Vec3d dv1 = dv + dv_dba * (update_ba);
        Vec3d dv2 = dv + dv_dba * (-update_ba);
        Vec3d dp1 = dp + dp_dba * (update_ba);
        Vec3d dp2 = dp + dp_dba * (-update_ba);
        //计算残差     
        Vec9d err1 = Eimu(posei, posej, vi, vj, dR, dv1, dp1);
        Vec9d err2 = Eimu(posei, posej, vi, vj, dR, dv2, dp2);
        //存数值求导结果
        J1.col(12+i) = (err1 - err2) / (2 * delta);
    }
    // posej     
    for (int i = 0; i < 3; ++i) {
        Vec3d update_se3 = Vec3d::Zero();
        update_se3(i) = delta;
        SE3 SE3_updated1 = posej;
        SE3 SE3_updated2 = posej;
        SE3_updated1.so3() = SE3_updated1.so3() * SO3::exp(update_se3);
        SE3_updated2.so3() = SE3_updated2.so3() * SO3::exp(-update_se3);
        Vec9d err1 = Eimu(posei, SE3_updated1, vi, vj, dR, dv, dp);
        Vec9d err2 = Eimu(posei, SE3_updated2, vi, vj, dR, dv, dp);
        J1.col(15+i) = (err1 - err2) / (2 * delta);
    }
    for (int i = 0; i < 3; ++i) {
        Vec3d update_se3 = Vec3d::Zero();
        update_se3(i) = delta;
        SE3 SE3_updated1 = posej;
        SE3 SE3_updated2 = posej;
        SE3_updated1.translation() += update_se3;
        SE3_updated2.translation() -= update_se3;        
        Vec9d err1 = Eimu(posei, SE3_updated1, vi, vj, dR, dv, dp);
        Vec9d err2 = Eimu(posei, SE3_updated2, vi, vj, dR, dv, dp);
        J1.col(18+i) = (err1 - err2) / (2 * delta);
    }
    // vj
    for (int i = 0; i < 3; ++i) {
        Vec3d update_v = Vec3d::Zero();
        update_v(i) = delta;
        Vec3d update_vj1 = vj + update_v;
        Vec3d update_vj2 = vj - update_v;
        Vec9d err1 = Eimu(posei, posej, vi, update_vj1, dR, dv, dp);
        Vec9d err2 = Eimu(posei, posej, vi, update_vj2, dR, dv, dp);
        J1.col(21+i) = (err1 - err2) / (2 * delta);
    }


    // 2. jacobian in SAD
    //公式求导 J2是IMU预积分残差的雅可比矩阵，为9X24维，初始化为0 9维是IMU预积分残差的九个维度，旋转预积分残差，速度预积分残差，位置预积分残差 24维是i时刻IMU的旋转 位置 速度 角速度零偏 加速度零偏 j时刻的旋转 位置 速度
    Eigen::Matrix<double, 9, 24> J2 = Eigen::Matrix<double, 9, 24>::Zero();

    // 一些中间符号
    const SO3 R1 = posei.so3();
    const SO3 R1T = R1.inverse();
    const SO3 R2 = posej.so3();

    // const SO3 dR = preint_->GetDeltaRotation(bg);
    const SO3 eR = SO3(dR).inverse() * R1T * R2;
    const Vec3d er = eR.log();
    const Mat3d invJr = SO3::jr_inv(eR);

    Vec3d pi = posei.translation();
    Vec3d pj = posej.translation();

    Eigen::Matrix<double, 9, 6> J_pose1;
    /// 残差对R1, 9x3
    J_pose1.setZero();
    // dR/dR1, 4.42
    J_pose1.block<3, 3>(0, 0) = -invJr * (R2.inverse() * R1).matrix();
    // dv/dR1, 4.47
    J_pose1.block<3, 3>(3, 0) = SO3::hat(R1T * (vj - vi - grav * dt));
    // dp/dR1, 4.48d
    J_pose1.block<3, 3>(6, 0) = SO3::hat(R1T * (pj - pi - vi * dt - 0.5 * grav * dt * dt));

    /// 残差对p1, 9x3
    // dp/dp1, 4.48a
    J_pose1.block<3, 3>(6, 3) = -R1T.matrix();

    Eigen::Matrix<double, 9, 3> J_v1;
    /// 残差对v1, 9x3
    J_v1.setZero();
    // dv/dv1, 4.46a
    J_v1.block<3, 3>(3, 0) = -R1T.matrix();
    // dp/dv1, 4.48c
    J_v1.block<3, 3>(6, 0) = -R1T.matrix() * dt;

    Eigen::Matrix<double, 9, 3> J_bg1;
    /// 残差对bg1
    J_bg1.setZero();
    // dR/dbg1, 4.45
    J_bg1.block<3, 3>(0, 0) = -invJr * eR.inverse().matrix() * SO3::jr((dR_dbg * Vec3d::Zero()).eval()) * dR_dbg;
    // dv/dbg1
    J_bg1.block<3, 3>(3, 0) = -dv_dbg;
    // dp/dbg1
    J_bg1.block<3, 3>(6, 0) = -dp_dbg;

    Eigen::Matrix<double, 9, 3> J_ba1;
    /// 残差对ba1
    J_ba1.setZero();
    // dv/dba1
    J_ba1.block<3, 3>(3, 0) = -dv_dba;
    // dp/dba1
    J_ba1.block<3, 3>(6, 0) = -dp_dba;

    Eigen::Matrix<double, 9, 6> J_pose2;
    /// 残差对pose2
    J_pose2.setZero();
    // dr/dr2, 4.43
    J_pose2.block<3, 3>(0, 0) = invJr;
    // dp/dp2, 4.48b
    J_pose2.block<3, 3>(6, 3) = R1T.matrix();

    Eigen::Matrix<double, 9, 3> J_v2;
    /// 残差对v2
    J_v2.setZero();
    // dv/dv2, 4,46b
    J_v2.block<3, 3>(3, 0) = R1T.matrix();  // OK


    J2.block<9, 6>(0, 0) = J_pose1;
    J2.block<9, 3>(0, 6) = J_v1;
    J2.block<9, 3>(0, 9) = J_bg1;
    J2.block<9, 3>(0, 12) = J_ba1;
    J2.block<9, 6>(0, 15) = J_pose2;
    J2.block<9, 3>(0, 21) = J_v2;

    //对比数值求导和公式求导的结果
    Eigen::Matrix<double, 9, 24> delta_J_IMU = J1 - J2;
    LOG(INFO) << "delta_J_IMU: ";
    LOG(INFO) << delta_J_IMU;
    LOG(INFO) << "           ";
    LOG(INFO) << "------------------------------------------------------------ ";




    //FIXME: GNSS  ******************************************************************** // 
    // Random initial value
    //随机设置传感器测量数值 用于和待优化变量相减作残差
    SE3 measurement = posej;
    measurement.so3() = posej.so3() * SO3::exp(omega);
    measurement.translation() << 0.45, 0.55, 0.5; 

    // 1. Numerical jacobian
    // J1 = ( f(x+dx) -f(x-dx) ) / (2 * dx)
    // 和IMU预积分一样思路 GNSS残差雅可比矩阵为6X6维 
    Eigen::Matrix<double, 6, 6> J3 = Eigen::Matrix<double, 6, 6>::Zero();
    // posej
    for (int i = 0; i < 3; ++i) {
        Vec3d update_se3 = Vec3d::Zero();
        update_se3(i) = delta;
        SE3 SE3_updated1 = posej;
        SE3 SE3_updated2 = posej;
        SE3_updated1.so3() = SE3_updated1.so3() * SO3::exp(update_se3);
        SE3_updated2.so3() = SE3_updated2.so3() * SO3::exp(-update_se3);
        Vec6d err1 = Egnss(SE3_updated1, measurement);
        Vec6d err2 = Egnss(SE3_updated2, measurement);
        J3.col(i) = (err1 - err2) / (2 * delta);
    }
    for (int i = 0; i < 3; ++i) {
        Vec3d update_se3 = Vec3d::Zero();
        update_se3(i) = delta;
        SE3 SE3_updated1 = posej;
        SE3 SE3_updated2 = posej;        
        SE3_updated1.translation() += update_se3;
        SE3_updated2.translation() -= update_se3;
        Vec6d err1 = Egnss(SE3_updated1, measurement);
        Vec6d err2 = Egnss(SE3_updated2, measurement);
        J3.col(3+i) = (err1 - err2) / (2 * delta);
    }

    // 2. jacobian in SAD
    //公式求导 参考课程代码即可
    Eigen::Matrix<double, 6, 6> J4 = Eigen::Matrix<double, 6, 6>::Zero();

    // jacobian 6x6
    J4.setZero();
    J4.block<3, 3>(0, 0) = (measurement.so3().inverse() * posej.so3()).jr_inv();  // dR/dR
    J4.block<3, 3>(3, 3) = Mat3d::Identity();          

    Eigen::Matrix<double, 6, 6> delta_J_GNSS = J3 - J4;
    LOG(INFO) << "         "; 
    LOG(INFO) << "deta_J_GNSS: " ;
    LOG(INFO) << delta_J_GNSS;    
    LOG(INFO) << "           ";
    LOG(INFO) << "------------------------------------------------------------ ";




    //FIXME: Prior  ******************************************************************** // 
    // Random initial value
    //随机设置先验数值 用于和待优化变量相减作残差
    NavStated state;
    state.timestamp_ = 0;
    state.p_ = posei.translation();
    state.R_ = posei.so3();
    state.v_.setZero();
    state.bg_ = bgi;
    state.ba_ = bai;

    //       | R1 | p1 | v1 | bg1 | ba1 |
    // 1. Numerical jacobian
    // J1 = ( f(x+dx) -f(x-dx) ) / (2 * dx)
    Eigen::Matrix<double, 15, 15> J5 = Eigen::Matrix<double, 15, 15>::Zero();
    // posei
    for (int i = 0; i < 3; ++i) {
        Vec3d update_se3 = Vec3d::Zero();
        update_se3(i) = delta;
        SE3 SE3_updated1 = posei;
        SE3 SE3_updated2 = posei;
        SE3_updated1.so3() = SE3_updated1.so3() * SO3::exp(update_se3);
        SE3_updated2.so3() = SE3_updated2.so3() * SO3::exp(-update_se3);
        Eigen::Matrix<double, 15, 1> err1 = Eprior(SE3_updated1, vi, bgi, bai, state);
        Eigen::Matrix<double, 15, 1> err2 = Eprior(SE3_updated2, vi, bgi, bai, state);
        J5.col(i) = (err1 - err2) / (2 * delta);
    }
    for (int i = 0; i < 3; ++i) {
        Vec3d update_se3 = Vec3d::Zero();
        update_se3(i) = delta;
        SE3 SE3_updated1 = posei;
        SE3 SE3_updated2 = posei;        
        SE3_updated1.translation() += update_se3;
        SE3_updated2.translation() -= update_se3;
        Eigen::Matrix<double, 15, 1> err1 = Eprior(SE3_updated1, vi, bgi, bai, state);
        Eigen::Matrix<double, 15, 1> err2 = Eprior(SE3_updated2, vi, bgi, bai, state);
        J5.col(3+i) = (err1 - err2) / (2 * delta);
    }
    // vi
    for (int i = 0; i < 3; ++i) {
        Vec3d update_v = Vec3d::Zero();
        update_v(i) = delta;
        Vec3d update_vi1 = vi + update_v;
        Vec3d update_vi2 = vi - update_v;
        Eigen::Matrix<double, 15, 1> err1 = Eprior(posei, update_vi1, bgi, bai, state);
        Eigen::Matrix<double, 15, 1> err2 = Eprior(posei, update_vi2, bgi, bai, state);
        J5.col(6+i) = (err1 - err2) / (2 * delta);
    }
    // bgi
    for (int i = 0; i < 3; ++i) {
        Vec3d update_bg = Vec3d::Zero();
        update_bg(i) = delta;
        Vec3d update_bgi1 = bgi + update_bg;
        Vec3d update_bgi2 = bgi - update_bg;
        Eigen::Matrix<double, 15, 1> err1 = Eprior(posei, vi, update_bgi1, bai, state);
        Eigen::Matrix<double, 15, 1> err2 = Eprior(posei, vi, update_bgi2, bai, state);        
        J5.col(9+i) = (err1 - err2) / (2 * delta);
    }      
    // bai
    for (int i = 0; i < 3; ++i) {
        Vec3d update_ba = Vec3d::Zero();
        update_ba(i) = delta;
        Vec3d update_bai1 = bgi + update_ba;
        Vec3d update_bai2 = bgi - update_ba;
        Eigen::Matrix<double, 15, 1> err1 = Eprior(posei, vi, bgi, update_bai1, state);
        Eigen::Matrix<double, 15, 1> err2 = Eprior(posei, vi, bgi, update_bai2, state);    
        J5.col(12+i) = (err1 - err2) / (2 * delta);
    }    


    // 2. jacobian in SAD
    //公式求导 参考课程代码即可
    Eigen::Matrix<double, 15, 15> J6 = Eigen::Matrix<double, 15, 15>::Zero();

    // jacobian 15x15
    const Vec3d er2 = SO3(state.R_.matrix().transpose() * posei.so3().matrix()).log();

    /// 注意有3个index, 顶点的，自己误差的，顶点内部变量的
    Eigen::Matrix<double, 15, 6> Jpp;
    Jpp.setZero();
    Jpp.block<3, 3>(0, 0) = SO3::jr_inv(er2);    // dr/dr
    Jpp.block<3, 3>(3, 3) = Mat3d::Identity();  // dp/dp

    Eigen::Matrix<double, 15, 3> Jvv;
    Jvv.setZero();
    Jvv.block<3, 3>(6, 0) = Mat3d::Identity();  // dv/dv

    Eigen::Matrix<double, 15, 3> Jbgg;
    Jbgg.setZero();
    Jbgg.block<3, 3>(9, 0) = Mat3d::Identity();  // dbg/dbg

    Eigen::Matrix<double, 15, 3> Jbaa;
    Jbaa.setZero();
    Jbaa.block<3, 3>(12, 0) = Mat3d::Identity();  // dba/dba

    J6.block<15, 6>(0, 0) = Jpp;
    J6.block<15, 3>(0, 6) = Jvv;
    J6.block<15, 3>(0, 9) = Jbgg;
    J6.block<15, 3>(0, 12) = Jbaa;


    Eigen::Matrix<double, 15, 15> delta_J_prior = J5 - J6;
    LOG(INFO) << "         "; 
    LOG(INFO) << "delta_J_prior: " ;
    LOG(INFO) << delta_J_prior;    
    LOG(INFO) << "           ";
    LOG(INFO) << "------------------------------------------------------------ ";


    // FIXME:
    // bg   ba   odom
    // These parameters all have good linearity.
    // (Jacobian == I or -I)
    // So there's no need to use numerical differentiation to verify accuracy.
    // Undoubtedly, their Jacobian matrix(I or -I) are correct.
    // Moreover, the gnss/prior/imu examples have been verified many times on identity matrix.

}