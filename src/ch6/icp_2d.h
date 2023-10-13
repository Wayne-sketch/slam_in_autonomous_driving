//
// Created by xiang on 2022/3/15.
//

#ifndef SLAM_IN_AUTO_DRIVING_ICP_2D_H
#define SLAM_IN_AUTO_DRIVING_ICP_2D_H

#include "common/eigen_types.h"
#include "common/lidar_utils.h"

#include <pcl/search/kdtree.h>

namespace sad {

/**
 * 第六章谈到的各种类型的ICP代码实现
 * 用法：先SetTarget, 此时构建target点云的KD树；再SetSource，然后调用Align*方法
 */
class Icp2d {
   public:
    using Point2d = pcl::PointXY;
    using Cloud2d = pcl::PointCloud<Point2d>;
    Icp2d() {}

    /// 设置目标的Scan
    void SetTarget(Scan2d::Ptr target) {
        target_scan_ = target;
        BuildTargetKdTree();
    }

    /// 设置被配准的Scan
    void SetSource(Scan2d::Ptr source) { source_scan_ = source; }

    /// 使用高斯牛顿法进行配准
    bool AlignGaussNewton(SE2& init_pose);

    /// 使用高斯牛顿法进行配准, Point-to-Plane
    bool AlignGaussNewtonPoint2Plane(SE2& init_pose);

    //作业：实现基于优化器的点到点 点到线ICP
    bool AlignG2OPoint2Point(SE2& init_pose);
    
    bool AlignG2OPoint2Plane(SE2& init_pose);

   private:
    // 建立目标点云的Kdtree
    void BuildTargetKdTree();

    /*pcl::search::KdTree 类是PCL库中的一部分，主要用于在点云数据中进行高效的近邻搜索。
    下面是一些在 pcl::search::KdTree 类中常用的成员函数：
    setInputCloud (const PointCloudConstPtr &cloud): 这个函数将一个点云设为搜索对象。传入的点云数据会被用来构造k-d树。
    nearestKSearch (const PointT &point, int k, std::vector<int> &k_indices, std::vector<float> &k_sqr_distances): 这个函数会在点云中搜索给定点的k个最近邻。返回的结果包括最近邻点的索引（k_indices）和它们到给定点的平方距离（k_sqr_distances）。
    radiusSearch (const PointT &point, double radius, std::vector<int> &k_indices, std::vector<float> &k_sqr_distances, unsigned int max_nn = 0): 这个函数会在给定点的某个半径范围内搜索点云中的点。返回的结果与 nearestKSearch 相同。
    请注意，上述函数中 PointT 是一个模板参数，它代表点云中点的数据类型。*/
    pcl::search::KdTree<Point2d> kdtree_;
    
    Cloud2d::Ptr target_cloud_;  // PCL 形式的target cloud

    Scan2d::Ptr target_scan_ = nullptr;
    Scan2d::Ptr source_scan_ = nullptr;
};

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_ICP_2D_H
