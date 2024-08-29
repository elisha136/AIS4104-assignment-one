#include <iostream>

#include <Eigen/Dense>

double deg_to_rad(double degrees)
{
    return degrees * 0.0174532925;
}

double rad_to_deg(double radians)
{
    return radians * 57.2957795;
}

Eigen::Vector3d eigen_matrix_types_example()
{
    Eigen::Matrix3d A; //assume this matrix has assigned values.
    A.setIdentity();
    Eigen::Matrix3d B; //assume this matrix has assigned values.
    B.setIdentity();
    Eigen::Vector3d b; //assume this vector has assigned values.
    b.setOnes();

    //The following demonstrates matrix-matrix multiplication and matrix-vector multiplication.
    Eigen::Matrix3d C = A * B;
    Eigen::Vector3d c = A * b;

    //The following demonstrates accessing elements of a matrix and vector.
    double A_11 = A(0, 0); // access the top-left element of the matrix.
    double c_x = c(0); // access the "x value" of c.
    double c_x2 = c.x(); // also accesses the "x value" of c; the functions .y() and .z() also exists.

    //The following demonstrates assignment of elements in a matrix and vector
    A(2, 2) = 3.0; // assign 1 to the bottom-right element of the matrix.
    c(2) = 3.0; // assign 3 to the "z value" of c.
    c << 1.0, 2.0, 3.0; //assign [1, 2, 3]^T to the vector.
}

Eigen::Matrix3d skew_symmetric (Eigen::Vector3d x){
 Eigen::Matrix3d skew_matrix;
    skew_matrix<<
        0.0, -x.z(), x.y(),
    x.z(),0.0, -x.x(),
    -x.y(),x.x(), 0;
    return skew_matrix;
}

void skew_symmetric_test()
{
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}

int main()
{
   skew_symmetric_test();
    return 0;
}
