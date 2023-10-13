//
// Created by BowenBZ on 2023/5/10.
//

#include "common/point_cloud_utils.h"
#include "common/point_types.h"

#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>

/// 点云的一些工具函数

namespace sad {

/// 体素滤波
//对点云数据进行体素化处理
//形参：点云数据 体素大小
void VoxelGrid(CloudPtr cloud, float voxel_size) {
    pcl::VoxelGrid<sad::PointType> voxel;
    //设置体素尺寸 X Y Z轴方向上的体素尺寸
    voxel.setLeafSize(voxel_size, voxel_size, voxel_size);
    //设置输入点云数据
    voxel.setInputCloud(cloud);

    CloudPtr output(new PointCloudType);
    //对输入点云进行体素滤波
    voxel.filter(*output);
    //交换指针，可以直接将滤波结果存储在原始的cloud指针指向的内存空间中，避免了额外的内存拷贝操作
    cloud->swap(*output);
}

/// 移除地面
void RemoveGround(CloudPtr cloud, float z_min) {
    CloudPtr output(new PointCloudType);
    for (const auto& pt : cloud->points) {
        if (pt.z > z_min) {
            output->points.emplace_back(pt);
        }
    }

    output->height = 1;
    output->is_dense = false;
    output->width = output->points.size();
    cloud->swap(*output);
}

/// 写点云文件
template<typename CloudType> 
void SaveCloudToFile(const std::string &filePath, CloudType &cloud) {
    cloud.height = 1;
    cloud.width = cloud.size();
    pcl::io::savePCDFileASCII(filePath, cloud);
}

template void SaveCloudToFile<PointCloudType>(const std::string &filePath, PointCloudType &cloud);

template void SaveCloudToFile<FullPointCloudType>(const std::string &filePath, FullPointCloudType &cloud);

}  // namespace sad