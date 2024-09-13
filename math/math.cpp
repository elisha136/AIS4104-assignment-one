#include <Eigen/Dense>

#include <iostream>

#include "math/math.h"



bool math::floatEquals(double a, double b) {
    return std::abs(a - b) < 1e-6;
}

Eigen::Matrix3d math:: rotate_x(double radians) {
    Eigen::Matrix3d matrix;
    matrix
            <<
            1.0, 0.0, 0.0,
            0.0, std::cos(radians), -std::sin(radians),
            0.0, std::sin(radians), std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::  rotate_y(double radians) {
    Eigen::Matrix3d matrix;
    matrix
            <<
            std::cos(radians), 0.0, std::sin(radians),
            0.0, 1.0, 0.0,
            -std::sin(radians), 0.0, std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math:: rotate_z(double radians) {
    Eigen::Matrix3d matrix;
    matrix
            <<
            std::cos(radians), -std::sin(radians), 0.0,
            std::sin(radians), std::cos(radians), 0.0,
            0.0, 0.0, 1.0;
    return matrix;
}

Eigen::Matrix3d  math:: rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
    return rotate_z(e.x()) * rotate_y(e.y()) * rotate_x(e.z());
}

Eigen::Vector3d math:: euler_zyx_from_rotation(const Eigen::Matrix3d &r) {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    if (floatEquals(r(2, 0), -1.0)) {
        b = EIGEN_PI / 2.0;
        a = 0.0;
        c = std::atan2(r(0, 1), r(1, 1));
    } else if (floatEquals(r(2, 0), 1.0)) {
        b = -(EIGEN_PI / 2.0);
        a = 0.0;
        c = -std::atan2(r(0, 1), r(1, 1));
    } else {
        b = std::atan2(-r(2, 0), std::sqrt(r(0, 0) * r(0, 0) + r(1, 0) * r(1, 0)));
        a = std::atan2(r(1, 0), r(0, 0));
        c = std::atan2(r(2, 1), r(2, 2));
    }
    return Eigen::Vector3d{a, b, c};
}





