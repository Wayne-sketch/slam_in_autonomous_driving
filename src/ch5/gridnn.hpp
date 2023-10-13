//
// Created by xiang on 2021/8/25.
//

#ifndef SLAM_IN_AUTO_DRIVING_GRID2D_HPP
#define SLAM_IN_AUTO_DRIVING_GRID2D_HPP

#include "common/eigen_types.h"
#include "common/math_utils.h"
#include "common/point_types.h"

#include <glog/logging.h>
#include <execution>
#include <map>

namespace sad {

/**
 * 栅格法最近邻
 * @tparam dim 模板参数，使用2D或3D栅格
 */
template <int dim>
class GridNN {
   public:
    using KeyType = Eigen::Matrix<int, dim, 1>;
    using PtType = Eigen::Matrix<float, dim, 1>;

    enum class NearbyType {
        CENTER,  // 只考虑中心
        // for 2D
        NEARBY4,  // 上下左右
        NEARBY8,  // 上下左右+四角

        // for 3D
        NEARBY6,  // 上下左右前后
        // for 3D
        NEARBY14,
    };

    /**
     * 构造函数
     * @param resolution 分辨率
     * @param nearby_type 近邻判定方法
     */
    explicit GridNN(float resolution = 0.1, NearbyType nearby_type = NearbyType::NEARBY4)
        : resolution_(resolution), nearby_type_(nearby_type) {
        inv_resolution_ = 1.0 / resolution_;

        // check dim and nearby
        if (dim == 2 && nearby_type_ == NearbyType::NEARBY6) {
            LOG(INFO) << "2D grid does not support nearby6, using nearby4 instead.";
            nearby_type_ = NearbyType::NEARBY4;
        } else if (dim == 3 && (nearby_type_ != NearbyType::NEARBY6 && nearby_type_ != NearbyType::CENTER)) {
            LOG(INFO) << "3D grid does not support nearby4/8, using nearby6 instead.";
            nearby_type_ = NearbyType::NEARBY6;
        }

        GenerateNearbyGrids();
    }

    /// 设置点云，建立栅格
    bool SetPointCloud(CloudPtr cloud);

    /// 获取最近邻
    bool GetClosestPoint(const PointType& pt, PointType& closest_pt, size_t& idx);

    /// 对比两个点云
    bool GetClosestPointForCloud(CloudPtr ref, CloudPtr query, std::vector<std::pair<size_t, size_t>>& matches);
    bool GetClosestPointForCloudMT(CloudPtr ref, CloudPtr query, std::vector<std::pair<size_t, size_t>>& matches);

   private:
    /// 根据最近邻的类型，生成附近网格
    void GenerateNearbyGrids();

    /// 空间坐标转到grid
    KeyType Pos2Grid(const PtType& pt);

    float resolution_ = 0.1;       // 分辨率
    float inv_resolution_ = 10.0;  // 分辨率倒数

    NearbyType nearby_type_ = NearbyType::NEARBY4;
    std::unordered_map<KeyType, std::vector<size_t>, hash_vec<dim>> grids_;  //  栅格数据
    CloudPtr cloud_;

    std::vector<KeyType> nearby_grids_;  // 附近的栅格
};

// 实现
template <int dim>
bool GridNN<dim>::SetPointCloud(CloudPtr cloud) {
    std::vector<size_t> index(cloud->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });

    std::for_each(index.begin(), index.end(), [&cloud, this](const size_t& idx) {
        auto pt = cloud->points[idx];
        auto key = Pos2Grid(ToEigen<float, dim>(pt));
        if (grids_.find(key) == grids_.end()) {
            grids_.insert({key, {idx}});
        } else {
            grids_[key].emplace_back(idx);
        }
    });

    cloud_ = cloud;
    LOG(INFO) << "grids: " << grids_.size();
    return true;
}

//传入点云坐标，传出所属栅格
template <int dim>
Eigen::Matrix<int, dim, 1> GridNN<dim>::Pos2Grid(const Eigen::Matrix<float, dim, 1>& pt) {
    return (pt * inv_resolution_).template cast<int>();
}

template <>
void GridNN<2>::GenerateNearbyGrids() {
    if (nearby_type_ == NearbyType::CENTER) {
        nearby_grids_.emplace_back(KeyType::Zero());
    } else if (nearby_type_ == NearbyType::NEARBY4) {
        nearby_grids_ = {Vec2i(0, 0), Vec2i(-1, 0), Vec2i(1, 0), Vec2i(0, 1), Vec2i(0, -1)};
    } else if (nearby_type_ == NearbyType::NEARBY8) {
        nearby_grids_ = {
            Vec2i(0, 0),   Vec2i(-1, 0), Vec2i(1, 0),  Vec2i(0, 1), Vec2i(0, -1),
            Vec2i(-1, -1), Vec2i(-1, 1), Vec2i(1, -1), Vec2i(1, 1),
        };
    }
}

template <>
void GridNN<3>::GenerateNearbyGrids() {
    if (nearby_type_ == NearbyType::CENTER) {
        nearby_grids_.emplace_back(KeyType::Zero());
    } else if (nearby_type_ == NearbyType::NEARBY6) {
        nearby_grids_ = {KeyType(0, 0, 0),  KeyType(-1, 0, 0), KeyType(1, 0, 0), KeyType(0, 1, 0),
                         KeyType(0, -1, 0), KeyType(0, 0, -1), KeyType(0, 0, 1)};
    }else if (nearby_type_ == NearbyType::NEARBY14) {
        nearby_grids_ = {KeyType(0, 0, 0),  KeyType(-1, 0, 0), KeyType(1, 0, 0), KeyType(0, 1, 0),
                         KeyType(0, -1, 0), KeyType(0, 0, -1), KeyType(0, 0, 1), KeyType(-1, -1, -1),
                         KeyType(-1, -1, 1), KeyType(1, -1, -1), KeyType(1, -1, 1), KeyType(-1, 1, -1),
                         KeyType(-1, 1, 1), KeyType(1, 1, -1), KeyType(1, 1, 1)};
    }
}

template <int dim>
bool GridNN<dim>::GetClosestPoint(const PointType& pt, PointType& closest_pt, size_t& idx) {
    // 在pt栅格周边寻找最近邻
    //存储待检查点云索引
    std::vector<size_t> idx_to_check;
    //计算pt所属栅格坐标
    auto key = Pos2Grid(ToEigen<float, dim>(pt));

    //遍历nearby_grids_，对于每个栅格偏移量，获取每个待检查栅格的
    std::for_each(nearby_grids_.begin(), nearby_grids_.end(), [&key, &idx_to_check, this](const KeyType& delta) {
        //待检查栅格的key
        auto dkey = key + delta;
        //找到对应栅格的所有点云索引
        auto iter = grids_.find(dkey);
        //把点云索引存入待检查点云索引
        if (iter != grids_.end()) {
            idx_to_check.insert(idx_to_check.end(), iter->second.begin(), iter->second.end());
        }
    });
    //如果idx_to_check为空，表示周围没有可检查的点，直接返回false
    if (idx_to_check.empty()) {
        return false;
    }

    // brute force nn in cloud_[idx]
    //新的点云对象指针
    CloudPtr nearby_cloud(new PointCloudType);
    //新的存储索引的容器
    std::vector<size_t> nearby_idx;
    //把idx_to_check中待检查的点的索引放入nearby_idx，把对应的点放入nearby_cloud
    for (auto& idx : idx_to_check) {
        nearby_cloud->points.template emplace_back(cloud_->points[idx]);
        nearby_idx.emplace_back(idx);
    }
    //暴力搜索，找到最近邻点在nearby_cloud中的索引，并存在closest_point_idx中
    size_t closest_point_idx = bfnn_point(nearby_cloud, ToVec3f(pt));
    //从nearby_idx中获取真实的索引值
    idx = nearby_idx.at(closest_point_idx);
    //获取最近邻点的坐标
    closest_pt = cloud_->points[idx];

    return true;
}

//为query中每个点找到参考点云ref中的最近邻点，并将匹配结果存储在matches向量中
template <int dim>
bool GridNN<dim>::GetClosestPointForCloud(CloudPtr ref, CloudPtr query,
                                          std::vector<std::pair<size_t, size_t>>& matches) {
    matches.clear();
    std::vector<size_t> index(query->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });
    std::for_each(index.begin(), index.end(), [this, &matches, &query](const size_t& idx) {
        PointType cp;
        size_t cp_idx;
        if (GetClosestPoint(query->points[idx], cp, cp_idx)) {
            matches.emplace_back(cp_idx, idx);
        }
    });

    return true;
}

//并行
template <int dim>
bool GridNN<dim>::GetClosestPointForCloudMT(CloudPtr ref, CloudPtr query,
                                            std::vector<std::pair<size_t, size_t>>& matches) {
    matches.clear();
    // 与串行版本基本一样，但matches需要预先生成，匹配失败时填入非法匹配
    std::vector<size_t> index(query->size());
    std::for_each(index.begin(), index.end(), [idx = 0](size_t& i) mutable { i = idx++; });
    matches.resize(index.size());

    std::for_each(std::execution::par_unseq, index.begin(), index.end(), [this, &matches, &query](const size_t& idx) {
        PointType cp;
        size_t cp_idx;
        if (GetClosestPoint(query->points[idx], cp, cp_idx)) {
            matches[idx] = {cp_idx, idx};
        } else {
            matches[idx] = {math::kINVALID_ID, math::kINVALID_ID};
        }
    });

    return true;
}

}  // namespace sad

#endif  // SLAM_IN_AUTO_DRIVING_GRID2D_HPP
