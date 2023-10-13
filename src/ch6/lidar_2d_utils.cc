//
// Created by xiang on 2022/3/15.
//

#include "ch6/lidar_2d_utils.h"
#include <opencv2/imgproc.hpp>

namespace sad {

void Visualize2DScan(Scan2d::Ptr scan, const SE2& pose, cv::Mat& image, const Vec3b& color, int image_size,
                     float resolution, const SE2& pose_submap) {
    //首先检查输入的图像是否为空，如果为空，则创建一个指定大小的图像，并将其背景设置为白色。
    if (image.data == nullptr) {
        image = cv::Mat(image_size, image_size, CV_8UC3, cv::Vec3b(255, 255, 255));
    }
    /*遍历扫描数据的每个激光点，对于有效的激光点（在范围内），执行以下操作：
    根据激光点的角度和距离计算其在二维空间中的坐标。
    检查激光点的角度是否在指定的范围内，如果不在范围内，则跳过该点。
    将激光点的坐标转换到子图坐标系下，并根据分辨率和图像大小计算其在图像中的像素坐标。
    如果像素坐标在图像范围内，则在图像上将该像素位置设置为指定的颜色。*/
    for (size_t i = 0; i < scan->ranges.size(); ++i) {
        if (scan->ranges[i] < scan->range_min || scan->ranges[i] > scan->range_max) {
            continue;
        }
        //scan中存的角度是相对于机器人坐标系，一般机器人正前方为0度，也是x轴正方向
        //这里转换为机器人坐标系的xy坐标
        double real_angle = scan->angle_min + i * scan->angle_increment;
        double x = scan->ranges[i] * std::cos(real_angle);
        double y = scan->ranges[i] * std::sin(real_angle);

        if (real_angle < scan->angle_min + 30 * M_PI / 180.0 || real_angle > scan->angle_max - 30 * M_PI / 180.0) {
            continue;
        }
        //将扫描点坐标转换到子地图坐标系下
        Vec2d psubmap = pose_submap.inverse() * (pose * Vec2d(x, y));

        int image_x = int(psubmap[0] * resolution + image_size / 2);
        int image_y = int(psubmap[1] * resolution + image_size / 2);
        if (image_x >= 0 && image_x < image.cols && image_y >= 0 && image_y < image.rows) {
            image.at<cv::Vec3b>(image_y, image_x) = cv::Vec3b(color[0], color[1], color[2]);
        }
    }
    //在图像上绘制当前位姿的位置，即在图像上以圆形标记出当前位姿的位置。
    // 同时画出pose自身所在位置
    Vec2d pose_in_image =
        pose_submap.inverse() * (pose.translation()) * double(resolution) + Vec2d(image_size / 2, image_size / 2);
    cv::circle(image, cv::Point2f(pose_in_image[0], pose_in_image[1]), 5, cv::Scalar(color[0], color[1], color[2]), 2);
}

}  // namespace sad