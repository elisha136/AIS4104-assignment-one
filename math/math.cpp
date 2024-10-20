#include <iostream>
#include <Eigen/Dense>
#include "math/math.h"
#include <functional>
#include <numeric>
#include <Vector>


bool math::floatEquals(double a, double b) {
    return std::abs(a - b) < 1e-6;
}

double math::deg_to_rad(double degrees) {
    return degrees * M_PI / 180.0;
}

Eigen::Matrix3d math::rotate_x(double degrees) {
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
            1.0, 0.0, 0.0,
            0.0, std::cos(radians), -std::sin(radians),
            0.0, std::sin(radians), std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::rotate_y(double degrees) {
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
            std::cos(radians), 0.0, std::sin(radians),
            0.0, 1.0, 0.0,
            -std::sin(radians), 0.0, std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::rotate_z(double degrees) {
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
            std::cos(radians), -std::sin(radians), 0.0,
            std::sin(radians), std::cos(radians), 0.0,
            0.0, 0.0, 1.0;
    return matrix;
}

Eigen::Matrix3d math::skew_symmetric(Eigen::Vector3d v) {
    Eigen::Matrix3d skewed_vector;
    skewed_vector <<
            0, -v.z(), v.y(),
            v.z(), 0, -v.x(),
            -v.y(), v.x(), 0;
    return skewed_vector;
}

Eigen::Matrix3d math::rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y,
                                                      const Eigen::Vector3d &z) {
    Eigen::Matrix3d matrix;
    matrix <<
            x(0), y(0), z(0),
            x(1), y(1), z(1),
            x(2), y(2), z(2);
    return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees) {
    Eigen::Matrix3d matrix{};
    double radians{deg_to_rad(degrees)};

    Eigen::Matrix3d skewed_vector{skew_symmetric(axis)};
    Eigen::Matrix3d I{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d skewed_vector_squared{skewed_vector * skewed_vector};

    // Bruker Rodrigues' formula for rot(omega,thetta)
    matrix = I + std::sin(radians) * skewed_vector + (1 - std::cos(radians)) * skewed_vector_squared;

    return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e) {
    Eigen::Matrix3d I{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d R_z{};
    Eigen::Matrix3d R_y{};
    Eigen::Matrix3d R_x{};
    Eigen::Matrix3d R{};

    R_z = rotate_z(e(0));
    R_y = rotate_y(e(1));
    R_x = rotate_x(e(2));

    R = I * R_z * R_y * R_x;

    return R;
}

Eigen::Matrix3d math::rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e) {
    const Eigen::Matrix3d I{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d R_z{};
    Eigen::Matrix3d R_y{};
    Eigen::Matrix3d R_x{};
    Eigen::Matrix3d R{};

    R_y = rotate_y(e(0));
    R_z = rotate_z(e(1));
    R_x = rotate_x(e(2));

    R = I * R_y * R_z * R_x;

    return R;
}

Eigen::Matrix4d math::transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p) {
    Eigen::Matrix4d matrix = Eigen::Matrix4d::Identity(); // Initialize as identity matrix

    // Set the top-left 3x3 block to the rotation matrix 'r'
    matrix.block<3, 3>(0, 0) = r;

    // Set the top-right 3x1 block to the translation vector 'p'
    matrix.block<3, 1>(0, 3) = p;

    // Bottom row is already correctly set to [0 0 0 1] by Eigen::Identity initialization

    return matrix;
}


Eigen::Vector3d math::euler_zyx_from_rotation(const Eigen::Matrix3d &r) {
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

/* Eigen::Vector3d math::euler_zyx_from_transformation(const Eigen::Matrix4d &t) {
    Eigen::Matrix3d R {};
    R << t(0, 0), t(0, 1), t(0, 2),
        t(1, 0), t(1, 1), t(1, 2),
        t(2, 0), t(2, 1), t(2, 2);

    return euler_zyx_from_rotation(R);
}*/

Eigen::VectorXd math::twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v) {
    Eigen::VectorXd twist(6);

    // Concatenate the angular and linear velocity vectors directly
    twist << w, v;

    return twist;
}


Eigen::VectorXd math::screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h) {
    Eigen::Vector3d w{s.x(), s.y(), s.z()};
    Eigen::Vector3d v{-s.cross(q) + h * s};
    Eigen::VectorXd S(6);
    S << w, v;
    return S;
}

Eigen::MatrixXd math::adjoint_matrix(const Eigen::Matrix4d &tf) {
    // Extract the 3x3 rotation matrix R from the top-left corner of tf
    Eigen::Matrix3d R = tf.block<3, 3>(0, 0);

    // Extract the translation vector p from the top-right corner of tf
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);

    // Compute the skew-symmetric matrix of p
    Eigen::Matrix3d skew_p = skew_symmetric(p);

    // Compute the matrix pR = skew(p) * R
    Eigen::Matrix3d pR = skew_p * R;

    // Create the 6x6 adjoint matrix
    Eigen::MatrixXd adjoint(6, 6);

    // Assign the top-left 3x3 block to R
    adjoint.block<3, 3>(0, 0) = R;

    // Assign the top-right 3x3 block to zero matrix
    adjoint.block<3, 3>(0, 3) = Eigen::Matrix3d::Zero();

    // Assign the bottom-left 3x3 block to pR
    adjoint.block<3, 3>(3, 0) = pR;

    // Assign the bottom-right 3x3 block to R
    adjoint.block<3, 3>(3, 3) = R;

    return adjoint;
}


double math::cot(double x) {
    return 1.0 / std::tan(x);
}

void math::wrench_in_s_and_w() {
    const Eigen::Vector3d f_w{-30, 0, 0};
    const Eigen::Vector3d m_s{0, 0, 2};
    const Eigen::Vector3d e_ws{60, -60, 0}; // Euler angles YZX in degrees.
    Eigen::Matrix3d R_ws{rotation_matrix_from_euler_yzx(e_ws)};

    Eigen::Vector3d m_w{R_ws * m_s};
    Eigen::Vector3d f_s{R_ws.transpose() * f_w};

    std::cout << "f_w: " << f_w.transpose() << std::endl;
    std::cout << "m_w: " << m_w.transpose() << std::endl;
    std::cout << "f_s: " << f_s.transpose() << std::endl;
    std::cout << "m_s: " << m_s.transpose() << std::endl;
    std::cout << "R_ws: " << std::endl << R_ws << std::endl << std::endl;
}

Eigen::VectorXd math::wrench_w(const Eigen::Vector3d &f_w, const Eigen::Vector3d &m_s, const Eigen::Vector3d &e_ws) {
    // Compute the rotation matrix from Euler angles (YZX order) for the transformation from sensor to world frame
    Eigen::Matrix3d R_ws = rotation_matrix_from_euler_yzx(e_ws);

    // Transform the torque from the sensor frame to the world frame
    Eigen::Vector3d m_w = R_ws * m_s;

    // Construct the wrench in the world frame [torque (m_w), force (f_w)]
    Eigen::VectorXd F_w(6);
    F_w << m_w, f_w;

    return F_w;
}


Eigen::VectorXd math::wrench_s(const Eigen::Vector3d &f_w, const Eigen::Vector3d &m_s, const Eigen::Vector3d &e_ws) {
    Eigen::Matrix3d R_ws{rotation_matrix_from_euler_yzx(e_ws)};

    Eigen::Vector3d f_s{R_ws.transpose() * f_w};

    Eigen::VectorXd F_s(6);
    F_s << m_s, f_s;

    return F_s;
}

Eigen::VectorXd math::wrench_a_to_b(const Eigen::VectorXd &F_a, const Eigen::Matrix4d &tf) {
    Eigen::MatrixXd Ad_T(6, 6);
    Ad_T = {math::adjoint_matrix(tf)};

    return Ad_T.transpose() * F_a;
}

// Task 2b
Eigen::VectorXd math::wrench_f(const Eigen::VectorXd &F_a, const Eigen::VectorXd &F_b, const Eigen::MatrixXd &tf_af,
                               const Eigen::MatrixXd &tf_bf) {
    return wrench_a_to_b(F_a, tf_af) + wrench_a_to_b(F_b, tf_bf);
}

// Task 3a
Eigen::Matrix3d math::matrix_exponential(const Eigen::Vector3d &w, const double theta) {
    const Eigen::Matrix3d skew_w{skew_symmetric(w)};
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = theta * M_PI / 180.0;

    return I + std::sin(rads) * skew_w + (1.0 - std::cos(rads)) * skew_w * skew_w;
}


// Task 3b
std::pair<Eigen::Vector3d, double> math::matrix_logarithm(const Eigen::Matrix3d &r) {
    double theta = 0.0;
    Eigen::Vector3d w = Eigen::Vector3d::Zero();
    const double trR = r.trace(); // Calculate the trace of the matrix (sum of diagonal elements)
    const double epsilon = 1e-6; // Tolerance for floating-point comparisons

    if (r.isApprox(Eigen::Matrix3d::Identity(), epsilon)) {
        // Case 1: Identity matrix (no rotation)
        theta = 0;
        w = Eigen::Vector3d::Zero();
    } else if (std::fabs(trR + 1) < epsilon) {
        // Case 2: 180-degree rotation (trace is -1)
        theta = EIGEN_PI;
        if (std::fabs(1 + r(2, 2)) > epsilon) {
            w = Eigen::Vector3d(r(0, 2), r(1, 2), 1 + r(2, 2)) / sqrt(2 * (1 + r(2, 2)));
        } else if (std::fabs(1 + r(1, 1)) > epsilon) {
            w = Eigen::Vector3d(r(0, 1), 1 + r(1, 1), r(2, 1)) / sqrt(2 * (1 + r(1, 1)));
        } else {
            w = Eigen::Vector3d(1 + r(0, 0), r(1, 0), r(2, 0)) / sqrt(2 * (1 + r(0, 0)));
        }
    } else {
        // Case 3: General case (rotation by arbitrary angle)
        theta = acos((trR - 1) / 2);
        Eigen::Matrix3d skew_w = (r - r.transpose()) / (2 * sin(theta));
        w = {skew_w(2, 1), skew_w(0, 2), skew_w(1, 0)};
    }

    // Return axis w and angle theta (in degrees)
    return std::make_pair(w, theta * 180.0 / EIGEN_PI);
}

//Task 3c. se(3) -> SE(3)
Eigen::Matrix4d math::matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, const double theta) {
    const Eigen::Matrix3d skew_w{skew_symmetric(w)};
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = theta * M_PI / 180;

    Eigen::Matrix3d R = matrix_exponential(w, theta);

    Eigen::Vector3d p = (I * rads + (1 - cos(rads)) * skew_w + (rads - sin(rads)) * skew_w * skew_w) * v;

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity(); // Initialize as identity matrix
    T.block<3, 3>(0, 0) = R; // Set the top-left 3x3 block to the rotation matrix R
    T.block<3, 1>(0, 3) = p; // Set the top-right 3x1 block to the translation vector p
    return T;
}

// Task 3d. SE(3) -> se(3)
Eigen::Matrix3d math::G(const Eigen::Vector3d &w, const double &theta) {
    Eigen::Matrix3d skew_w{skew_symmetric(w)};
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = theta * M_PI / 180;

    return I * rads + (1 - cos(rads)) * skew_w + (rads - sin(rads)) * skew_w * skew_w;
}

Eigen::Matrix3d math::G_inverse(const Eigen::Vector3d &w, const double &degrees) {
    Eigen::Matrix3d skew_w{skew_symmetric(w)};
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = degrees * M_PI / 180;

    return I / rads - skew_w / 2 + (1 / rads - cot(rads / 2) / 2) * skew_w * skew_w;
}


std::pair<Eigen::VectorXd, double> math::matrix_logarithm(const Eigen::Matrix4d &t) {
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double degrees;

    Eigen::Matrix3d R;
    R << t(0, 0), t(0, 1), t(0, 2),
            t(1, 0), t(1, 1), t(1, 2),
            t(2, 0), t(2, 1), t(2, 2);
    Eigen::Vector3d p{t(0, 3), t(1, 3), t(2, 3)};

    if (R == Eigen::Matrix3d::Identity()) {
        w = {0, 0, 0};
        v = p.normalized();
        degrees = sqrt(p(0) * p(0) + p(1) * p(1) + p(2) * p(2));
    } else {
        auto [fst, scd] = math::matrix_logarithm(R);
        w = fst;
        degrees = scd;

        const Eigen::Matrix3d G_inv = math::G_inverse(w, degrees);
        v = G_inv * p;
    }
    Eigen::VectorXd S(6);
    S << w(0), w(1), w(2), v(0), v(1), v(2);

    return std::make_pair(S, degrees);
}

void math::print_pose(const std::string &label, const Eigen::Matrix4d &tf) {
    Eigen::Matrix3d R;
    R << tf(0, 0), tf(0, 1), tf(0, 2),
            tf(1, 0), tf(1, 1), tf(1, 2),
            tf(2, 0), tf(2, 1), tf(2, 2);
    Eigen::Vector3d p{tf(0, 3), tf(1, 3), tf(2, 3)};

    Eigen::Vector3d e_zyx = euler_zyx_from_rotation(R);

    std::cout << "Label: " << label << std::endl;
    std::cout << "Euler ZYX angles: " << e_zyx.transpose() * 180 / EIGEN_PI << std::endl;
    std::cout << "Linear position: " << p.transpose() << std::endl;
    std::cout << " " << std::endl;
}

// Task 4b
Eigen::Matrix4d math::planar_3r_fk_transform(const std::vector<double> &joint_positions) {
    constexpr double L1 = 10.0, L2 = 10.0, L3 = 10.0;

    // Making transformation matrix for each succeeding frame. Every joint rotates around z-axis.
    const Eigen::Matrix4d T01 = transformation_matrix(rotate_z(joint_positions[0]), Eigen::Vector3d(0, 0, 0));
    const Eigen::Matrix4d T12 = transformation_matrix(rotate_z(joint_positions[1]), Eigen::Vector3d(L1, 0, 0));
    const Eigen::Matrix4d T23 = transformation_matrix(rotate_z(joint_positions[2]), Eigen::Vector3d(L2, 0, 0));
    const Eigen::Matrix4d T34 = transformation_matrix(rotate_z(0), Eigen::Vector3d(L3, 0, 0));

    Eigen::Matrix4d T04 = T01 * T12 * T23 * T34;

    return T04;
}

void math::test_planar_3r_fk_transform(const std::string &label, const std::vector<double> &joint_positions) {
    const Eigen::Matrix4d T{math::planar_3r_fk_transform(joint_positions)};

    print_pose(label, T);
}

// Task 4c
Eigen::Matrix4d math::planar_3r_fk_screw(const std::vector<double> &joint_positions) {
    constexpr double L1 = 10, L2 = 10, L3 = 10;

    Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
    M(0, 3) = L1 + L2 + L3;

    Eigen::Vector3d w1, w2, w3;
    w1 << 0, 0, 1;
    w2 << 0, 0, 1;
    w3 << 0, 0, 1;

    Eigen::Vector3d v1, v2, v3;
    v1 << 0, 0, 0;
    v2 << 0, -L1, 0;
    v3 << 0, -L1 - L2, 0;

    const Eigen::Matrix4d e1 = matrix_exponential(w1, v1, joint_positions[0]);
    const Eigen::Matrix4d e2 = matrix_exponential(w2, v2, joint_positions[1]);
    const Eigen::Matrix4d e3 = matrix_exponential(w3, v3, joint_positions[2]);

    const Eigen::Matrix4d T04 = e1 * e2 * e3 * M;

    return T04;
}

void math::test_planar_3r_fk_screw(const std::string &label, const std::vector<double> &joint_positions) {
    const Eigen::Matrix4d T{planar_3r_fk_screw(joint_positions)};

    print_pose(label, T);
}


// TASK 5
// Task 5a
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd> > math::ur3e_space_chain() {
    //OK
    constexpr double h1{0.15185}, l1{0.24355}, l2{0.2132}, h2{0.08535},
            y1{0.13105}, y2{0.0921};

    Eigen::Matrix4d M; //ok
    M << 1, 0, 0, -l1 - l2,
            0, 0, -1, -y1 - y2,
            0, 1, 0, h1 - h2,
            0, 0, 0, 1;

    Eigen::Vector3d w0, w1, w2, w3, w4, w5;
    w0 << 0, 0, 1; //(base) //ok
    w1 << 0, -1, 0; //ok
    w2 << 0, -1, 0; //ok
    w3 << 0, -1, 0; //ok
    w4 << 0, 0, -1; //ok
    w5 << 0, -1, 0; //ok

    Eigen::Vector3d v0, v1, v2, v3, v4, v5;
    v0 << 0, 0, 0;
    v1 << h1, 0, 0;
    v2 << h1, 0, l1;
    v3 << h1, 0, l1 + l2;
    v4 << y1, -l1 - l2, 0;
    v5 << h1 - h2, 0, l1 + l2;

    std::vector<Eigen::VectorXd> screws{
        twist(w0, v0),
        twist(w1, v1),
        twist(w2, v2),
        twist(w3, v3),
        twist(w4, v4),
        twist(w5, v5)

    };

    return std::make_pair(M, screws);
}


Eigen::Matrix4d math::ur3e_space_fk(const Eigen::VectorXd &joint_positions) {
    auto [m, space_screws] = math::ur3e_body_chain();
    Eigen::Matrix4d t06 = Eigen::Matrix4d::Identity();
    for (int i = 0; i < joint_positions.size(); i++) {
        Eigen::VectorXd screw = space_screws[i];
        Eigen::Vector3d w {screw(0), screw(1), screw(2)};
        Eigen::Vector3d v {screw(3), screw(4), screw(5)};
        t06 *= math::matrix_exponential(w, v, joint_positions[i]);

    }


    return t06 * m;
}


Eigen::Matrix4d math::ur3e_body_fk(const Eigen::VectorXd &joint_positions) {
    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> par =  math::ur3e_body_chain();// auto [m, body_screws] = math::ur3e_body_chain();
    Eigen::Matrix4d t06 = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d m = par.first;
    std::vector<Eigen::VectorXd> body_screws = par.second;
    for (int i = 0; i < joint_positions.size(); i++) {
        Eigen::VectorXd screw = body_screws[i];
        Eigen::Vector3d w {screw(0), screw(1), screw(2)};
        Eigen::Vector3d v {screw(3), screw(4), screw(5)};
        t06 *= math::matrix_exponential(w, v, joint_positions[i]);
    }
    return t06 * m;
}



//task 1d ur3e body chain



std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd> > math::ur3e_body_chain() {
    // Get the space frame data (M and screw axes)
    auto [M, space_screws] = ur3e_space_chain();

    // Compute the inverse of M (for transforming from space frame to body frame)
    Eigen::Matrix4d M_inv = M.inverse();

    // Compute the adjoint transformation of M_inv
    Eigen::MatrixXd adj_M_inv = adjoint_matrix(M_inv);

    // Transform the space screw axes into the body frame
    std::vector<Eigen::VectorXd> body_screws;
    for (const auto &screw: space_screws) {
        Eigen::VectorXd body_screw = adj_M_inv * screw; // Apply adjoint transformation
        body_screws.push_back(body_screw);
    }

    // Return the matrix M and the body-frame screw axes
    return std::make_pair(M, body_screws);
}



/*void math::test_ur3e_fk_screw(const std::string &label, const std::vector<double> &joint_positions) {
   const Eigen::Matrix4d T = math::ur3e_fk_screw(joint_positions);

  print_pose(label, T);
}*/

Eigen::Matrix4d math::DH_transformation_matrix(const double &joint_angle, const double &alpha, const double &a,
                                         const double &d) {
    Eigen::Matrix4d T_mn;
    double joint_rads = joint_angle * M_PI / 180.0;
    double alpha_rads = alpha * M_PI / 180.0;
    T_mn << cos(joint_rads), -sin(joint_rads) * cos(alpha_rads), sin(joint_rads) * sin(alpha_rads), a * cos(joint_rads),
            //ok
            sin(joint_rads), cos(joint_rads) * cos(alpha_rads), -cos(joint_rads) * sin(alpha_rads), a * sin(joint_rads),
            0.0, sin(alpha_rads), cos(alpha_rads), d,
            0.0, 0.0, 0.0, 1.0;
    return T_mn;
}

Eigen::Matrix4d math::ur3e_fk_transform(const std::vector<double> &joint_positions) {
    //YES!
    constexpr double h0{0.15185}, h1{0.24355}, h2{0.2132}, h3{0.08535},
            y0{0.13105}, y1{0.0921};

    const Eigen::Matrix4d T01 = DH_transformation_matrix(joint_positions[0], 90, 0, h0);
    const Eigen::Matrix4d T12 = DH_transformation_matrix(joint_positions[1], 0, -h1, 0);
    const Eigen::Matrix4d T23 = DH_transformation_matrix(joint_positions[2], 0, -h2, 0);
    const Eigen::Matrix4d T34 = DH_transformation_matrix(joint_positions[3], 90, 0, y0);
    const Eigen::Matrix4d T45 = DH_transformation_matrix(joint_positions[4], -90, 0, h3);
    const Eigen::Matrix4d T56 = DH_transformation_matrix(joint_positions[5], 0, 0, y1);


    Eigen::Matrix4d T06 = T01 * T12 * T23 * T34 * T45 * T56;

    return T06;
}

void math::test_ur3e_fk_transform(const std::string &label, const std::vector<double> &joint_positions) {
    const Eigen::Matrix4d T{ur3e_fk_transform(joint_positions)};

    print_pose(label, T);
}


Eigen::VectorXd math::std_vector_to_eigen1(const std::vector<double> &v) {
    Eigen::Map<const Eigen::VectorXd> eigenvec(v.data(), v.size());
    return eigenvec;
}


//Assignment 3
//TASK 1a
Eigen::VectorXd math::std_vector_to_eigen(const std::vector<double> &v) {
    Eigen::VectorXd r(v.size());
    for (int i = 0; i < v.size(); ++i)
        r(i) = v[i];
    return r;
}


//Task 1b
bool math::is_average_below_eps(const std::vector<double> &values, double eps, uint8_t n_values) {
    if (values.size() < n_values)
        return false;
    const double sum = std::accumulate(values.end(), values.end(), 0.0);
    return (sum / n_values) < eps;
}

//1c


//Task 2 Newton_Raphson and Gradient_Descent methods
// Newton-Raphson method
std::pair<uint32_t, double> math::newton_raphson_root_find(const std::function<double(double)> &f, double x_0,
                                                           double dx_0, double eps,
                                                           uint32_t max_iters) {
    double x = x_0;
    uint32_t iter = 0;
    while (iter < max_iters) {
        double fx = f(x);
        double dfx = (f(x + dx_0) - f(x)) / dx_0;
        if (std::abs(dfx) < eps) {
            std::cout << "Derivative is too small" << std::endl;
            break;
        }
        double x_next = x - fx / dfx;
        if (std::abs(x_next - x) < eps) {
            x = x_next;
            break;
        }
        x = x_next;
        iter++;
    }
    return std::make_pair(iter, x);
}

// Gradient Descent method
std::pair<uint32_t, double> math::gradient_descent_root_find(const std::function<double(double)> &f, double x_0,
                                                             double gamma, double dx_0, double eps,
                                                             uint32_t max_iters) {
    double x = x_0;
    uint32_t iter = 0;
    while (iter < max_iters) {
        double fx = f(x);
        double dfx = (f(x + dx_0) - fx) / dx_0;
        if (std::abs(dfx) < eps) {
            break;
        }
        double x_next = x - gamma * dfx;
        if (std::abs(x_next - x) < eps) {
            x = x_next;
            break;
        }
        x = x_next;
        iter++;
    }
    return std::make_pair(iter, x);
}

// Test Newton-Raphson root finding
void math::test_newton_raphson_root_find(const std::function<double(double)> &f, double x0) {
    auto [iterations, x_hat] = newton_raphson_root_find(f, x0, 0.5, 1E-07, 1000);
    std::cout << "NR root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) <<
            std::endl;
}

// Test Gradient Descent root finding
void math::test_gradient_descent_root_find(const std::function<double(double)> &f, double x0) {
    auto [iterations, x_hat] = gradient_descent_root_find(f, x0, 0.5, 1E-07, 1000);
    std::cout << "GD root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) <<
            std::endl;
}

// Root-finding tests
void math::test_root_find() {
    std::cout << "Root finding tests" << std::endl;
    auto f1 = [](double x) {
        return (x - 3.0) * (x - 3.0) - 1.0;
    };
    math::test_newton_raphson_root_find(f1, -20);
    math::test_gradient_descent_root_find(f1, -20);
}

//task 3 space jacobian
Eigen::MatrixXd math::ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions) {
    // Get the number of joints (6)
    int num_joints = current_joint_positions.size();

    // Get the space chain (M and screws)
    auto [M, screws] = math::ur3e_space_chain();
    Eigen::MatrixXd J_s(6, num_joints);
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    for (int i = 0; i < num_joints; ++i) {
        if (i == 0) {
            J_s.col(i) = screws[i];
        } else {
            Eigen::MatrixXd Ad_T = Eigen::MatrixXd::Identity(6, 6);
            J_s.col(i) = Ad_T * screws[i];
        }
        Eigen::Matrix4d exp_twist = Eigen::Matrix4d::Identity();
        T = T * exp_twist;
    }

    return J_s;
}

// task 3 body jacobian

Eigen::MatrixXd math::ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions) {
    int num_joints = current_joint_positions.size();
    auto [M, body_screws] = ur3e_body_chain();
    Eigen::MatrixXd j_b(6, num_joints);
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    for (int i = num_joints - 1; i >= 0; --i) {
        if (i == num_joints - 1) {
            j_b.col(i) = body_screws[i];
        } else {
            Eigen::MatrixXd Ad_T = adjoint_matrix(T);
            j_b.col(i) = Ad_T * body_screws[i];
        }
        Eigen::VectorXd screw = body_screws[i];
        Eigen::Vector3d w {screw(0), screw(1), screw(2)};
        Eigen::Vector3d v {screw(3), screw(4), screw(5)};

        Eigen::Matrix4d exp_twist = matrix_exponential( w,v,current_joint_positions[i]);
        T = exp_twist * T;
    }

    return j_b;
}


// space and body jacobian test functions
void math::ur3e_test_jacobian(const Eigen::VectorXd &joint_positions) {
    Eigen::Matrix4d tsb = math::ur3e_body_fk(joint_positions);
    auto [m, space_screws] = math::ur3e_space_chain();
    Eigen::MatrixXd jb = math::ur3e_body_jacobian(joint_positions);
    Eigen::MatrixXd js = math::ur3e_space_jacobian(joint_positions);
    Eigen::MatrixXd ad_tsb = math::adjoint_matrix(tsb);
    Eigen::MatrixXd ad_tbs = math::adjoint_matrix(tsb.inverse());
    std::cout << "Jb: " << std::endl << jb << std::endl
            << "Ad_tbs*Js:" << std::endl << ad_tbs * js << std::endl << std::endl;
    std::cout << "Js: " << std::endl << js << std::endl
            << "Ad_tsb*Jb:" << std::endl << ad_tsb * jb << std::endl << std::endl;
    std::cout << "d Jb: " << std::endl << jb - ad_tbs * js << std::endl << std::endl;
    std::cout << "d Js: " << std::endl << js - ad_tsb * jb << std::endl << std::endl;
}

void math::ur3e_test_jacobian() {
    std::cout << "Jacobian matrix tests" << std::endl;
    ur3e_test_jacobian(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * M_PI / 180);
    ur3e_test_jacobian(
        math::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}) * M_PI / 180);
}


std::pair<size_t, Eigen::VectorXd> math::ur3e_ik_body(const Eigen::Matrix4d &t_sd,
                                                      const Eigen::VectorXd &current_joint_positions,
                                                      double gamma, double v_e, double w_e) {
    const int max_iterations = 1000;
    int iterations = 0;
    Eigen::VectorXd joint_positions = current_joint_positions;
    std::pair<Eigen::VectorXd, double> pair;
    Eigen::Matrix4d t_error;
    Eigen::Matrix3d R_error;
    Eigen::VectorXd twist(6);
    Eigen::MatrixXd J_b(6, 6);
    double w_error=1;
    double v_error=1;
    while ((v_error > v_e and w_error > w_e) or (iterations < max_iterations)) {
        t_error = math::ur3e_body_fk(joint_positions).inverse() * t_sd;
        R_error = t_error.block<3, 3>(0, 0);
        pair = math::matrix_logarithm(t_error);
        twist = pair.first * pair.second;
        w_error = twist.head(3).norm();
        v_error = twist.tail(3).norm();
        J_b = math::ur3e_body_jacobian(joint_positions);
        joint_positions += gamma * J_b.completeOrthogonalDecomposition().pseudoInverse() * twist;
        iterations++;
    }
    if (iterations == max_iterations) {
        std::cerr << "Warning: IK solver did not converge after " << max_iterations << " iterations." << std::endl;
    }
    return std::make_pair(iterations, joint_positions);
}


void math::ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd &j0) {
    std::cout << "Test from pose" << std::endl;
    Eigen::Matrix4d t_sd = math::transformation_matrix(math::rotation_matrix_from_euler_zyx(zyx), pos);
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = math::ur3e_body_fk(j_ik);
    math::print_pose(" IK pose",t_ik );
    math::print_pose( "Desired pose",t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * M_PI / 180 << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * M_PI / 180 << std::endl << std::endl;
}

void math::ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0) {
    std::cout << "Test from configuration" << std::endl;
    Eigen::Matrix4d t_sd = math::ur3e_space_fk(joint_positions);
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = math::ur3e_body_fk(j_ik);
    math::print_pose(" IK pose",t_ik );
    math::print_pose("Desired pose",t_sd );
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * M_PI / 180 << std::endl;
    std::cout << "J_d: " << joint_positions.transpose() * M_PI / 180 << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * M_PI / 180 << std::endl << std::endl;
}

void math::ur3e_ik_test() {
    Eigen::VectorXd j_t0 = math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * M_PI / 180;
    Eigen::VectorXd j_t1 = math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -89.0, 0.0, 0.0, 0.0}) * M_PI / 180;

    math::ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} * M_PI / 180, j_t0);
    math::ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} * M_PI / 180, j_t1);

    Eigen::VectorXd j_t2 = math::std_vector_to_eigen(std::vector<double>{50.0, -30.0, 20.0, 0.0, -30.0, 50.0}) * M_PI /
                           180;
    Eigen::VectorXd j_d1 = math::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}) * M_PI /
                           180;

    ur3e_ik_test_configuration(j_d1, j_t0);
    ur3e_ik_test_configuration(j_d1, j_t2);
}

