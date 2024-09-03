#include <iostream>

#include <Eigen/Dense>

double deg_to_rad(double degrees) {
    return degrees * 0.0174532925;
}

double rad_to_deg(double radians) {
    return radians * 57.2957795;
}

Eigen::Vector3d eigen_matrix_types_example() {
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

Eigen::Matrix3d skew_symmetric(Eigen::Vector3d x) {
    Eigen::Matrix3d skew_matrix;
    skew_matrix <<
            0.0, -x.z(), x.y(),
            x.z(), 0.0, -x.x(),
            -x.y(), x.x(), 0;
    return skew_matrix;
}

void skew_symmetric_test() {
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}

Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                const Eigen::Vector3d &y,
                                                const Eigen::Vector3d &z) {
    Eigen::Matrix3d matrix;


    matrix.col(0) = x;
    matrix.col(1) = y;
    matrix.col(2) = z;


    return matrix;
}


Eigen::Matrix3d rotate_x(double gamma) {
    Eigen::Matrix3d matrix;

    double rad_value = gamma * 0.0174532925;
    matrix <<
            1, 0, 0,
            0, std::cos(rad_value), -std::sin(rad_value),
            0, std::sin(rad_value), std::cos(rad_value);
    return matrix;
}

Eigen::Matrix3d rotate_y(double beta) {
    Eigen::Matrix3d matrix;

    double rad_value = beta * 0.0174532925;

    matrix <<
            std::cos(rad_value), 0, std::sin(rad_value),
            0, 1, 0,
            -std::sin(rad_value), 0, std::cos(rad_value);
    return matrix;
}


Eigen::Matrix3d rotate_z(double alpha) {
    Eigen::Matrix3d matrix;

    double rad_value = alpha * 0.0174532925;
    matrix <<


            std::cos(rad_value), -std::sin(rad_value), 0,
            std::sin(rad_value), std::cos(rad_value), 0,
            0, 0, 1;
    return matrix;
}

Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees) {
    Eigen::Matrix3d matrix;
    double rad_value = deg_to_rad(degrees);
    double c = std::cos(rad_value);
    double s = std::sin(rad_value);

    double w1 = axis(0);
    double w2 = axis(1);
    double w3 = axis(2);

    matrix <<
            c + (w1 * w1) * (1 - c), w1 * w2 * (1 - c) - w3 * s, w1 * w3 * (1 - c) + w2 * s,
            w1 * w2 * (1 - c) + w3 * s, c + (w2 * w2) * (1 - c), w2 * w3 * (1 - c) - w1 * s,
            w1 * w3 * (1 - c) - w2 * s, w2 * w3 * (1 - c) + w1 * s, c + (w3 * w3) * (1 - c);


    return matrix;
}

Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
    //
    double Z = e(0);
    double Y = e(1);
    double X = e(2);


    Eigen::Matrix3d matrix;


    //Euler ZYX
    matrix = rotate_z(Z) * rotate_y(Y) * rotate_x(X);

    return matrix;
}

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p) {
    Eigen::Matrix4d matrix;


    matrix <<
            r(0, 0), r(0, 1), r(0, 2), p(0),
            r(1, 0), r(1, 1), r(1, 2), p(1),
            r(2, 0), r(2, 1), r(2, 2), p(2),
            0.0, 0.0, 0.0, 1.0;

    return matrix;
}

// Function to transform vector from {a} to {w} coordinates
void transform_vector() {
    // Given Euler angles e = [60, 45, 0] in degrees
    Eigen::Vector3d euler_angles(60, 45, 0);

    // Create the rotation matrix from these Euler angles
    Eigen::Matrix3d matrix = rotation_matrix_from_euler_zyx(euler_angles);

    // Given vector va = [2.5, 3.0, -10] in {a} coordinates
    Eigen::Vector3d va(2.5, 3.0, -10);

    // Rotate va to get vw_rotated in {w} coordinates
    Eigen::Vector3d vw_rotated = matrix * va;

    // Translate the rotated vector by the position of frame {a} (which is 10 units along the z-axis in {w} coordinates)
    Eigen::Vector3d translation(0, 0, 10);
    Eigen::Vector3d vw = vw_rotated + translation;

    // Print the transformed vector vw in {w} coordinates
    std::cout << "Transformed vector in {w} coordinates: \n" << vw << std::endl;
}


void rotation_matrix_test() {
    Eigen::Matrix3d rot =
            rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
    Eigen::Matrix3d rot_aa =
            rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0, 0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa =
            rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
                                            Eigen::Vector3d{-0.5, -0.5, 0.707107},
                                            Eigen::Vector3d{0.707107, -0.707107, 0.0});
    std::cout << "Rotation matrix from Euler: " << std::endl;
    std::cout << rot << std::endl << std::endl;
    std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
    std::cout << rot_aa << std::endl << std::endl;
    std::cout << "Rotation matrix from frame axes: " << std::endl;
    std::cout << rot_fa << std::endl << std::endl;
}

void transformation_matrix_test() {
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
    Eigen::Vector3d v{1.0, -2.0, 3.0};
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << transformation_matrix(r, v) << std::endl;
}


int main() {
    skew_symmetric_test();
    rotation_matrix_test();
    transformation_matrix_test();
    transform_vector();
    return 0;
}
