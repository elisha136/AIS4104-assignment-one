
#ifndef MATH_H
#define MATH_H


#include <Eigen/Dense>

namespace math
{
    constexpr double deg_to_rad = 57.2957795;
    constexpr double rad_to_deg = 0.0174532925;

    Eigen::Matrix3d rotate_x(double gamma);

    Eigen::Matrix3d rotate_y(double beta);

    Eigen::Matrix3d rotate_z(double alpha);

    bool floatEquals(double a, double b);
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
    Eigen::Matrix3d  rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
    Eigen::Vector3d euler_zyx_from_rotation(const Eigen::Matrix3d &r) ;


}
#endif //MATH_H
