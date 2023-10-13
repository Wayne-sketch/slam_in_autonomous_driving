#ifndef SLAM_IN_AUTO_DRIVING_NANOFLANN_TEST_H
#define SLAM_IN_AUTO_DRIVING_NANOFLANN_TEST_H

#include <iostream>
#include <vector>
#include <nanoflann.hpp>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

typedef pcl::PointXYZI PointType;
typedef pcl::PointCloud<PointType> PointCloud;

class NanoflannNearestNeighbor {
public:
    NanoflannNearestNeighbor(const PointCloud::Ptr& referencePoints)
        : referencePoints_(referencePoints) {
        buildIndex();
    }

    void match(const PointCloud::Ptr& first, const PointCloud::Ptr& second,
               std::vector<std::pair<size_t, size_t>>& matches) {
        for (size_t i = 0; i < second->size(); ++i) {
            const PointType& queryPoint = (*second)[i];

            // 使用nanoflann进行最近邻搜索
            std::vector<size_t> indices(1);
            std::vector<float> distances(1);
            nanoflann::KNNResultSet<float> resultSet(1);
            resultSet.init(&indices[0], &distances[0]);
            index_.findNeighbors(resultSet, queryPoint.data, nanoflann::SearchParams());

            // 将匹配结果存储到matches变量中
            size_t queryIndex = i;
            size_t matchIndex = indices[0];
            matches.emplace_back(queryIndex, matchIndex);
        }
    }

private:
    PointCloud::Ptr referencePoints_;
    typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<float, PointCloud>,
        PointCloud,
        3,
        size_t
    > KDTree;
    KDTree index_;

    void buildIndex() {
        index_.buildIndex();
    }
};

int main() {
    // 加载点云数据
    PointCloud::Ptr first(new PointCloud);
    PointCloud::Ptr second(new PointCloud);
    pcl::io::loadPCDFile("path/to/first.pcd", *first);
    pcl::io::loadPCDFile("path/to/second.pcd", *second);

    std::vector<std::pair<size_t, size_t>> matches;

    // 创建最近邻匹配对象
    NanoflannNearestNeighbor nn(second);

    // 执行匹配
    nn.match(first, second, matches);

    // 输出匹配结果
    for (const auto& match : matches) {
        std::cout << "First point index: " << match.first << ", Second point index: " << match.second << std::endl;
    }

    return 0;
}

#endif  // SLAM_IN_AUTO_DRIVING_NANOFLANN_TEST_H