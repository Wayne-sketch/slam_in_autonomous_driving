//
// Created by xiang on 2022/3/18.
//

#include "ch6/g2o_types.h"
#include "ch6/likelihood_filed.h"

#include <glog/logging.h>

#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/robust_kernel.h>
#include <g2o/core/robust_kernel_impl.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>
#include <g2o/solvers/dense/linear_solver_dense.h>

namespace sad {
//设置目标扫描数据，并生成对应的似然场函数。
void LikelihoodField::SetTargetScan(Scan2d::Ptr scan) {
    //将输入的目标扫描数据指针（scan）赋值给成员变量 target_
    target_ = scan;

    // 在target点上生成场函数
    //创建一个大小为 1000x1000 的 CV_32F 类型的图像 field_，用于存储似然场函数。
    field_ = cv::Mat(1000, 1000, CV_32F, 30.0);
    //遍历目标扫描数据的每个激光点
    for (size_t i = 0; i < scan->ranges.size(); ++i) {
        if (scan->ranges[i] < scan->range_min || scan->ranges[i] > scan->range_max) {
            continue;
        }
        //计算激光点在似然场图像中的坐标（x 和 y），根据激光点的距离、角度和分辨率进行计算，同时将坐标原点移到图像中心。
        double real_angle = scan->angle_min + i * scan->angle_increment;
        double x = scan->ranges[i] * std::cos(real_angle) * resolution_ + 500;
        double y = scan->ranges[i] * std::sin(real_angle) * resolution_ + 500;

        // 在(x,y)附近填入场函数
        //针对每个模板点（model_pt），在以激光点坐标为中心的区域内更新似然场函数。
        for (auto& model_pt : model_) {
            int xx = int(x + model_pt.dx_);
            int yy = int(y + model_pt.dy_);
            if (xx >= 0 && xx < field_.cols && yy >= 0 && yy < field_.rows &&
                field_.at<float>(yy, xx) > model_pt.residual_) {
                field_.at<float>(yy, xx) = model_pt.residual_;
            }
        }
    }
}

void LikelihoodField::BuildModel() {
    const int range = 20;  // 生成多少个像素的模板
    for (int x = -range; x <= range; ++x) {
        for (int y = -range; y <= range; ++y) {
            model_.emplace_back(x, y, std::sqrt((x * x) + (y * y)));
        }
    }
}

void LikelihoodField::SetSourceScan(Scan2d::Ptr scan) { source_ = scan; }

bool LikelihoodField::AlignGaussNewton(SE2& init_pose) {
    int iterations = 10;
    double cost = 0, lastCost = 0;
    SE2 current_pose = init_pose;
    const int min_effect_pts = 20;  // 最小有效点数
    const int image_boarder = 20;   // 预留图像边界

    has_outside_pts_ = false;
    for (int iter = 0; iter < iterations; ++iter) {
        Mat3d H = Mat3d::Zero();
        Vec3d b = Vec3d::Zero();
        cost = 0;

        int effective_num = 0;  // 有效点数

        // 遍历source
        for (size_t i = 0; i < source_->ranges.size(); ++i) {
            float r = source_->ranges[i];
            if (r < source_->range_min || r > source_->range_max) {
                continue;
            }

            float angle = source_->angle_min + i * source_->angle_increment;
            if (angle < source_->angle_min + 30 * M_PI / 180.0 || angle > source_->angle_max - 30 * M_PI / 180.0) {
                continue;
            }

            float theta = current_pose.so2().log();
            Vec2d pw = current_pose * Vec2d(r * std::cos(angle), r * std::sin(angle));

            // 在field中的图像坐标
            Vec2i pf = (pw * resolution_ + Vec2d(500, 500)).cast<int>();

            if (pf[0] >= image_boarder && pf[0] < field_.cols - image_boarder && pf[1] >= image_boarder &&
                pf[1] < field_.rows - image_boarder) {
                effective_num++;

                // 图像梯度
                float dx = 0.5 * (field_.at<float>(pf[1], pf[0] + 1) - field_.at<float>(pf[1], pf[0] - 1));
                float dy = 0.5 * (field_.at<float>(pf[1] + 1, pf[0]) - field_.at<float>(pf[1] - 1, pf[0]));
                //公式6.22
                Vec3d J;
                J << resolution_ * dx, resolution_ * dy,
                    -resolution_ * dx * r * std::sin(angle + theta) + resolution_ * dy * r * std::cos(angle + theta);
                H += J * J.transpose();

                float e = field_.at<float>(pf[1], pf[0]);
                b += -J * e;

                cost += e * e;
            } else {
                has_outside_pts_ = true;
            }
        }

        if (effective_num < min_effect_pts) {
            return false;
        }

        // solve for dx
        Vec3d dx = H.ldlt().solve(b);
        if (isnan(dx[0])) {
            break;
        }

        cost /= effective_num;
        if (iter > 0 && cost >= lastCost) {
            break;
        }

        LOG(INFO) << "iter " << iter << " cost = " << cost << ", effect num: " << effective_num;

        current_pose.translation() += dx.head<2>();
        current_pose.so2() = current_pose.so2() * SO2::exp(dx[2]);
        lastCost = cost;
    }

    init_pose = current_pose;
    return true;
}

//用于获取似然场的图像表示
cv::Mat LikelihoodField::GetFieldImage() {
    //创建一个大小与似然场函数 field_ 相同的图像 image，类型为 CV_8UC3。
    cv::Mat image(field_.rows, field_.cols, CV_8UC3);
    for (int x = 0; x < field_.cols; ++x) {
        for (int y = 0; y < field_.rows; ++y) {
            //计算当前像素的值 r，将其映射到0-255的范围，以便进行图像显示。
            //将像素值 r 转换为 CV_8UC3 类型的颜色值，并赋值给图像 image 的对应像素位置。
            float r = field_.at<float>(y, x) * 255.0 / 30.0;
            image.at<cv::Vec3b>(y, x) = cv::Vec3b(uchar(r), uchar(r), uchar(r));
        }
    }

    return image;
}

//使用g2o库进行扫描数据的配准，并返回配准结果。
bool LikelihoodField::AlignG2O(SE2& init_pose) {
    //定义g2o中所需的块求解器（BlockSolverType）和线性求解器（LinearSolverType），并创建相应的实例。
    //每个误差项优化变量维度为3，误差值维度为1
    using BlockSolverType = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
    using LinearSolverType = g2o::LinearSolverCholmod<BlockSolverType::PoseMatrixType>;
    //创建一个稀疏优化器（optimizer）并设置使用的优化算法为Levenberg-Marquardt算法。
    auto* solver = new g2o::OptimizationAlgorithmLevenberg(
        g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(solver);
    //创建一个SE2类型的顶点（v），设置其初始位姿为传入的初始位姿（init_pose），并将其添加到优化器中。
    auto* v = new VertexSE2();
    v->setId(0);
    v->setEstimate(init_pose);
    optimizer.addVertex(v);

    const double range_th = 15.0;  // 不考虑太远的scan，不准
    const double rk_delta = 0.8;

    has_outside_pts_ = false;
    // 遍历source
    for (size_t i = 0; i < source_->ranges.size(); ++i) {
        float r = source_->ranges[i];
        if (r < source_->range_min || r > source_->range_max) {
            continue;
        }

        if (r > range_th) {
            continue;
        }

        float angle = source_->angle_min + i * source_->angle_increment;
        if (angle < source_->angle_min + 30 * M_PI / 180.0 || angle > source_->angle_max - 30 * M_PI / 180.0) {
            continue;
        }

        auto e = new EdgeSE2LikelihoodFiled(field_, r, angle, resolution_);
        e->setVertex(0, v);

        if (e->IsOutSide()) {
            has_outside_pts_ = true;
            delete e;
            continue;
        }

        e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());
        auto rk = new g2o::RobustKernelHuber;
        rk->setDelta(rk_delta);
        e->setRobustKernel(rk);
        optimizer.addEdge(e);
    }

    optimizer.setVerbose(false);
    optimizer.initializeOptimization();
    optimizer.optimize(10);

    init_pose = v->estimate();
    return true;
}

//用于从占据栅格地图生成似然场图像。
void LikelihoodField::SetFieldImageFromOccuMap(const cv::Mat& occu_map) {
    //定义边界大小（boarder）用于避免处理地图边缘
    const int boarder = 25;
    //创建一个大小为 1000x1000 的 CV_32F 类型的图像 field_，并将初始值设置为 30.0。
    field_ = cv::Mat(1000, 1000, CV_32F, 30.0);
    //遍历占据栅格地图的每个像素
    for (int x = boarder; x < occu_map.cols - boarder; ++x) {
        for (int y = boarder; y < occu_map.rows - boarder; ++y) {
            //检查当前像素的灰度值是否小于 127，表示该像素为占据栅格。
            if (occu_map.at<uchar>(y, x) < 127) {
                // 在该点生成一个model
                //如果是占据栅格，根据当前像素的位置生成似然场模板。
                for (auto& model_pt : model_) {
                    int xx = int(x + model_pt.dx_);
                    int yy = int(y + model_pt.dy_);
                    if (xx >= 0 && xx < field_.cols && yy >= 0 && yy < field_.rows &&
                        field_.at<float>(yy, xx) > model_pt.residual_) {
                        field_.at<float>(yy, xx) = model_pt.residual_;
                    }
                }
            }
        }
    }
}

}  // namespace sad