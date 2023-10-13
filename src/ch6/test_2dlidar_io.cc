//
// Created by xiang on 2022/3/15.
//
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <opencv2/highgui.hpp>

#include "ch6/lidar_2d_utils.h"
#include "common/io_utils.h"

DEFINE_string(bag_path, "./dataset/sad/2dmapping/test_2d_lidar.bag", "数据包路径");

/// 测试从rosbag中读取2d scan并plot的结果

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = google::INFO;
    FLAGS_colorlogtostderr = true;
    google::ParseCommandLineFlags(&argc, &argv, true);
    //可视化程序负责将雷达的距离和角度信息转化为笛卡尔坐标
    //如果我们输入机器人的位姿，这个绘制程序还可以在运动状态下演示
    sad::RosbagIO rosbag_io(fLS::FLAGS_bag_path);
    rosbag_io
        .AddScan2DHandle("/pavo_scan_bottom",
                         [](Scan2d::Ptr scan) {
                             cv::Mat image;
                             sad::Visualize2DScan(scan, SE2(), image, Vec3b(255, 0, 0));
                             cv::imshow("scan", image);
                             cv::waitKey(20);
                             return true;
                         })
        .Go();

    return 0;
}