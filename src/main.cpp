#include <Eigen/Dense>
#include <iostream>
#include "math/math.h"
#include <iostream>

#include <Eigen/Dense>


int main()
{
    Eigen::Vector3d e = Eigen::Vector3d{60.0, 45.0, 30.0} *math:: deg_to_rad;
    Eigen::Matrix3d r =math:: rotation_matrix_from_euler_zyx(e);
    Eigen::Vector3d ea = math::euler_zyx_from_rotation(r);
    std::cout << " E: " << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Ea: " << ea.transpose() *math:: rad_to_deg << std::endl;

    return 0;
}