//
// Created by xiang on 2022/3/15.
//

#include "ch6/icp_2d.h"
#include "common/math_utils.h"

#include <glog/logging.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/search/impl/kdtree.hpp>

namespace sad {
//作业
typedef g2o::BlockSolver<g2o::BlockSolverTraits<3, 3>> BlockSolverType;
typedef g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType> LinearSolverType;
typedef g2o::OptimizationAlgorithmLevenberg OptimizationAlgorithmType;

bool Icp2d::AlignG2OPoint2Point(SE2& init_pose) {
    // 创建优化器
    g2o::SparseOptimizer optimizer;
    BlockSolverType* solver_ptr = new BlockSolverType(std::unique_ptr<LinearSolverType>(new LinearSolverType()));
    OptimizationAlgorithmType* solver = new OptimizationAlgorithmType(solver_ptr);
    optimizer.setAlgorithm(solver);

    // 创建顶点 待优化的当前位姿
    VertexSE2* v = new VertexSE2();
    v->setId(0);
    v->setEstimate(init_pose);
    optimizer.addVertex(v);

    const float max_dis2 = 0.01;    // 最近邻时的最远距离（平方）
    const int min_effect_pts = 20;  // 最小有效点数
    int effective_num = 0;  // 有效点数

    // 遍历source
    for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
        //遍历2D激光雷达扫描数据的测量数据
        float r = source_scan_->ranges[i];
        //数据不在测量范围内，跳过
        if (r < source_scan_->range_min || r > source_scan_->range_max) {
            continue;
        }
        //获取当前scan激光点对应的角度（当前的机器人坐标系）
        float angle = source_scan_->angle_min + i * source_scan_->angle_increment;
        //获取当前机器人位姿角度
        SE2 current_pose = init_pose;
        float theta = current_pose.so2().log();
        //测试后和公式结果一样 计算前一帧scan坐标系下的检测点 当前scan坐标系的激光点坐标转换到前一帧scan坐标系下
        Vec2d pw = current_pose * Vec2d(r * std::cos(angle), r * std::sin(angle));
        Point2d pt;
        pt.x = pw.x();
        pt.y = pw.y();

        // 最近邻
        std::vector<int> nn_idx;
        std::vector<float> dis;
        //搜索检测点的最近邻 nn_idx最近邻点的索引 dis最近邻点到检测点的平方距离
        //可是建树的时候近邻 nn_idx最近邻点的索引 dis最近邻点到检测点的平方距离
        kdtree_.nearestKSearch(pt, 1, nn_idx, dis);

        if (nn_idx.size() > 0 && dis[0] < max_dis2) {
            // 有效点数+1
            effective_num++;
            // 创建一元边，测量值是上一scan中对应最近邻点的坐标 qi_w
            EdgeSE2_Point2Point* edge = new EdgeSE2_Point2Point(angle, theta, r);
            // edge->setId(i);
            edge->setVertex(0, v);
            //测量值是target点云中pt对应的最近邻点 数据类型是Point2d
            //测量值是上一scan最近邻点坐标
            edge->setMeasurement(target_cloud_->points[nn_idx[0]]);
            edge->setInformation(Eigen::Matrix2d::Identity());  // 设置边的信息矩阵

            optimizer.addEdge(edge);
        }
    }

    // 有效点不够多
    if (effective_num < min_effect_pts) {
        return false;
    }

    // 解决优化问题
    optimizer.initializeOptimization();
    optimizer.optimize(10);
    // 优化后的位姿存入传入的init_pose（引用）
    init_pose = v->estimate();

    return true;
}
    
bool Icp::AlignG2OPoint2Plane(SE2& init_pose){

}


//基于高斯牛顿迭代的 2D ICP 形参：初始位姿 SE2
bool Icp2d::AlignGaussNewton(SE2& init_pose) {
    //迭代次数
    int iterations = 10;
    double cost = 0, lastCost = 0;
    //初始位姿存入当前位姿
    SE2 current_pose = init_pose;
    const float max_dis2 = 0.01;    // 最近邻时的最远距离（平方）
    const int min_effect_pts = 20;  // 最小有效点数
    //遍历迭代次数 ，一次迭代是一次高斯牛顿
    for (int iter = 0; iter < iterations; ++iter) {
        //3X3double类型矩阵
        Mat3d H = Mat3d::Zero();
        //3X1double类型向量
        Vec3d b = Vec3d::Zero();
        cost = 0;

        int effective_num = 0;  // 有效点数

        // 遍历source
        for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
            //遍历2D激光雷达扫描数据的测量数据
            float r = source_scan_->ranges[i];
            //数据不在测量范围内，跳过
            if (r < source_scan_->range_min || r > source_scan_->range_max) {
                continue;
            }
            //获取当前观测对应的角度（机器人坐标系）
            float angle = source_scan_->angle_min + i * source_scan_->angle_increment;
            //获取当前机器人位姿角度
            float theta = current_pose.so2().log();
            //测试后和公式结果一样 计算世界坐标系下的检测点 机器人坐标系的激光点坐标转换到世界坐标系下
            Vec2d pw = current_pose * Vec2d(r * std::cos(angle), r * std::sin(angle));
            Point2d pt;
            pt.x = pw.x();
            pt.y = pw.y();

            // 最近邻
            std::vector<int> nn_idx;
            std::vector<float> dis;
            //搜索检测点的最近邻 nn_idx最近邻点的索引 dis最近邻点到检测点的平方距离
            //可是建树的时候近邻 nn_idx最近邻点的索引 dis最近邻点到检测点的平方距离
            kdtree_.nearestKSearch(pt, 1, nn_idx, dis);

            if (nn_idx.size() > 0 && dis[0] < max_dis2) {
                //有效点数+1
                effective_num++;
                //3X2double矩阵 公式6.9 雅可比矩阵
                Mat32d J;
                J << 1, 0, 0, 1, -r * std::sin(angle + theta), r * std::cos(angle + theta);
                //高斯牛顿法
                //3X3矩阵
                H += J * J.transpose();
                //残差
                Vec2d e(pt.x - target_cloud_->points[nn_idx[0]].x, pt.y - target_cloud_->points[nn_idx[0]].y);
                b += -J * e;
                //计算e的点积
                cost += e.dot(e);
            }
        }
        //有效点不够多 返回false
        if (effective_num < min_effect_pts) {
            return false;
        }
        //这行代码使用 LDLT 分解方法解算出向量 dx。H 是一个矩阵，ldlt() 是对 H 进行 LDLT 分解的方法，solve(b) 是使用 LDLT 分解解算方程 H * dx = b 的方法。
        // solve for dx
        Vec3d dx = H.ldlt().solve(b);
        //这行代码检查向量 dx 的第一个元素是否为非数值（NaN）。如果是非数值，则跳出循环。
        if (isnan(dx[0])) {
            break;
        }
        //这行代码将 cost 的值除以 effective_num，用于计算每个数据点的平均成本。effective_num 是一个表示有效数据点数量的变量。
        cost /= effective_num;
        //这行代码用于检查当前迭代的成本 cost 是否大于或等于上一次迭代的成本 lastCost。如果是这样，那么代码会跳出循环，即结束迭代过程。
        if (iter > 0 && cost >= lastCost) {
            break;
        }
        //显示当前已迭代次数 和有效点
        LOG(INFO) << "iter " << iter << " cost = " << cost << ", effect num: " << effective_num;
        //更新当前优化位置
        current_pose.translation() += dx.head<2>();
        current_pose.so2() = current_pose.so2() * SO2::exp(dx[2]);
        lastCost = cost;
    }
    //优化后的位姿存入 传入的init_pose（引用）
    init_pose = current_pose;
    LOG(INFO) << "estimated pose: " << current_pose.translation().transpose()
              << ", theta: " << current_pose.so2().log();

    return true;
}

bool Icp2d::AlignGaussNewtonPoint2Plane(SE2& init_pose) {
    int iterations = 10;
    double cost = 0, lastCost = 0;
    SE2 current_pose = init_pose;
    const float max_dis = 0.3;      // 最近邻时的最远距离
    const int min_effect_pts = 20;  // 最小有效点数

    for (int iter = 0; iter < iterations; ++iter) {
        Mat3d H = Mat3d::Zero();
        Vec3d b = Vec3d::Zero();
        cost = 0;

        int effective_num = 0;  // 有效点数

        // 遍历source
        for (size_t i = 0; i < source_scan_->ranges.size(); ++i) {
            float r = source_scan_->ranges[i];
            if (r < source_scan_->range_min || r > source_scan_->range_max) {
                continue;
            }

            float angle = source_scan_->angle_min + i * source_scan_->angle_increment;
            float theta = current_pose.so2().log();
            Vec2d pw = current_pose * Vec2d(r * std::cos(angle), r * std::sin(angle));
            Point2d pt;
            pt.x = pw.x();
            pt.y = pw.y();

            // 查找5个最近邻
            //对一个激光点，找五个最近邻点，拟合出一条直线，把该激光点到直线的垂直距离当作残差
            std::vector<int> nn_idx;
            std::vector<float> dis;
            kdtree_.nearestKSearch(pt, 5, nn_idx, dis);

            std::vector<Vec2d> effective_pts;  // 有效点
            for (int j = 0; j < nn_idx.size(); ++j) {
                if (dis[j] < max_dis) {
                    effective_pts.emplace_back(
                        Vec2d(target_cloud_->points[nn_idx[j]].x, target_cloud_->points[nn_idx[j]].y));
                }
            }

            if (effective_pts.size() < 3) {
                continue;
            }

            // 拟合直线，组装J、H和误差
            Vec3d line_coeffs;
            //effective_pts 有效点容器  line_coeffs拟合出的直线参数
            if (math::FitLine2D(effective_pts, line_coeffs)) {
                //该激光点有效
                effective_num++;
                Vec3d J;
                //公式6.17
                J << line_coeffs[0], line_coeffs[1],
                    -line_coeffs[0] * r * std::sin(angle + theta) + line_coeffs[1] * r * std::cos(angle + theta);
                H += J * J.transpose();
                //公式6.14
                double e = line_coeffs[0] * pw[0] + line_coeffs[1] * pw[1] + line_coeffs[2];
                b += -J * e;

                cost += e * e;
            }
        }
        //如果某一次高斯牛顿迭代，有效的激光点数目不足，迭代失败
        if (effective_num < min_effect_pts) {
            return false;
        }

        // solve for dx
        //使用LDLT分解求解线性方程组，得到位姿更新量 dx。
        Vec3d dx = H.ldlt().solve(b);
        //如果 dx 中包含NaN值（表示求解失败），则跳出迭代循环。
        if (isnan(dx[0])) {
            break;
        }
        //将代价 cost 除以有效点的数量 effective_num，得到平均代价。
        cost /= effective_num;
        //如果当前迭代次数大于0且当前代价 cost 大于等于上一次代价 lastCost，则跳出迭代循环。
        if (iter > 0 && cost >= lastCost) {
            break;
        }
        //打印当前迭代的次数、代价和有效点的数量的日志信息。
        LOG(INFO) << "iter " << iter << " cost = " << cost << ", effect num: " << effective_num;
        //更新当前位姿的平移部分，将 dx 的前两个元素加到当前位姿的平移向量上。
        current_pose.translation() += dx.head<2>();
        //更新当前位姿的旋转部分，将 dx 的第三个元素作为旋转的增量，通过乘以一个旋转变换 SO2::exp(dx[2])，将其应用到当前位姿的旋转上
        current_pose.so2() = current_pose.so2() * SO2::exp(dx[2]);
        //将当前代价 cost 赋值给上一次代价 lastCost。
        lastCost = cost;
    }

    init_pose = current_pose;
    LOG(INFO) << "estimated pose: " << current_pose.translation().transpose()
              << ", theta: " << current_pose.so2().log();

    return true;
}

void Icp2d::BuildTargetKdTree() {
    //目标的scan 首先检查 target_scan_ 是否为空指针，如果为空则输出错误信息并返回
    if (target_scan_ == nullptr) {
        LOG(ERROR) << "target is not set";
        return;
    }

    target_cloud_.reset(new Cloud2d);
    //使用循环遍历目标扫描数据中的每个点
    for (size_t i = 0; i < target_scan_->ranges.size(); ++i) {
        //在循环中，首先检查目标扫描数据的某个点是否在指定的范围内（range_min 和 range_max 之间），如果不在范围内则跳过该点。
        if (target_scan_->ranges[i] < target_scan_->range_min || target_scan_->ranges[i] > target_scan_->range_max) {
            continue;
        }
        //计算目标点的实际角度(机器人坐标系)，根据扫描数据的 angle_min、i 和 angle_increment 计算得到
        double real_angle = target_scan_->angle_min + i * target_scan_->angle_increment;
        //body系下的xy坐标
        Point2d p;
        p.x = target_scan_->ranges[i] * std::cos(real_angle);
        p.y = target_scan_->ranges[i] * std::sin(real_angle);
        //将计算得到的点 p 添加到 target_cloud_->points 中，即将点云数据添加到目标点云中
        target_cloud_->points.push_back(p);
    }
    //设置目标点云的宽度为 target_cloud_->points.size()，即点云中点的数量。
    //将目标点云标记为非稠密（is_dense = false），表示点云中可能存在无效点。
    target_cloud_->width = target_cloud_->points.size();
    target_cloud_->is_dense = false;
    //设置kd树点云
    kdtree_.setInputCloud(target_cloud_);
}

}  // namespace sad