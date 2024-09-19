#include <iostream>
#include <Eigen/Dense>
#include "math/math.h"



bool math::floatEquals(double a, double b) {
    return std::abs(a - b) < 1e-6;
}

double math::deg_to_rad(double degrees)
{
    return degrees * M_PI / 180;
}

Eigen::Matrix3d math::rotate_x(double degrees)
{
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
        1, 0, 0,
        0, std::cos(radians), -std::sin(radians),
        0, std::sin(radians), std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::rotate_y(double degrees) {
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
         std::cos(radians), 0, std::sin(radians),
                         0, 1,                  0,
        -std::sin(radians), 0,  std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::rotate_z(double degrees) {
    double radians = deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix <<
        std::cos(radians), -std::sin(radians), 0,
        std::sin(radians),  std::cos(radians), 0,
                        0,                  0, 1;
    return matrix;
}

Eigen::Matrix3d math::skew_symmetric(Eigen::Vector3d v) {
    Eigen::Matrix3d skewed_vector;
    skewed_vector <<
             0, -v.z(),  v.y(),
         v.z(),      0, -v.x(),
        -v.y(),  v.x(),      0;
    return skewed_vector;
}

Eigen::Matrix3d math::rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z)
{
Eigen::Matrix3d matrix;
matrix <<
    x(0), y(0), z(0),
    x(1), y(1), z(1),
    x(2), y(2), z(2);
return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees)
{
    Eigen::Matrix3d matrix {};
    double radians {deg_to_rad(degrees)};

    Eigen::Matrix3d skewed_vector {skew_symmetric(axis)};
    Eigen::Matrix3d I {Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d skewed_vector_squared {skewed_vector * skewed_vector};

    // Bruker Rodrigues' formula for rot(omega,thetta)
    matrix = I + std::sin(radians) * skewed_vector + (1-std::cos(radians))*skewed_vector_squared;

    return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d I {Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d R_z {};
    Eigen::Matrix3d R_y {};
    Eigen::Matrix3d R_x {};
    Eigen::Matrix3d R {};

    R_z = rotate_z(e(0));
    R_y = rotate_y(e(1));
    R_x = rotate_x(e(2));

    R = I * R_z * R_y * R_x;

    return R;
}

Eigen::Matrix3d math::rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e)
{
    const Eigen::Matrix3d I {Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d R_z {};
    Eigen::Matrix3d R_y {};
    Eigen::Matrix3d R_x {};
    Eigen::Matrix3d R {};

    R_y = rotate_y(e(0));
    R_z = rotate_z(e(1));
    R_x = rotate_x(e(2));

    R = I * R_y * R_z * R_x;

    return R;
}

Eigen::Matrix4d math::transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix4d matrix;
    // implement the necessary equations and functionality.
    matrix <<
        r(0,0), r(0,1), r(0,2), p(0),
        r(1,0), r(1,1), r(1,2), p(1),
        r(2,0), r(2,1), r(2,2), p(2),
                    0,              0,            0,          1;
    return matrix;
}

Eigen::Vector3d math::euler_zyx_from_rotation(const Eigen::Matrix3d &r)
{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    if(floatEquals(r(2, 0), -1.0))
    {
        b = EIGEN_PI / 2.0;
        a = 0.0;
        c = std::atan2(r(0, 1), r(1, 1));
    }
    else if(floatEquals(r(2, 0), 1.0))
    {
        b = -(EIGEN_PI / 2.0);
        a = 0.0;
        c = -std::atan2(r(0, 1), r(1, 1));
    }
    else
    {
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
    twist << w(0), w(1), w(2), v(0), v(1), v(2);
    return twist;
}

Eigen::VectorXd math::screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h) {
    Eigen::Vector3d w {s.x(), s.y(), s.z()};
    Eigen::Vector3d v {-s.cross(w) + h*s};
    Eigen::VectorXd S(6);
    S << w(0), w(1), w(2), v(0), v(1), v(2);
    return S;
}

Eigen::MatrixXd math::adjoint_matrix(const Eigen::Matrix4d &tf) {
    Eigen::Vector3d p {tf(0, 3), tf(1, 3), tf(2, 3)};
    Eigen::Matrix3d skew_p {skew_symmetric(p)};
    Eigen::Matrix3d R;
        R << tf(0,0), tf(0,1), tf(0,2),
             tf(1,0), tf(1,1), tf(1,2),
             tf(2,0), tf(2,1), tf(2,2);
    Eigen::Matrix3d pR {skew_p*R};

    Eigen::MatrixXd adjoint(6, 6);
    adjoint <<  R(0,0),  R(0,1),  R(0,2),  0,             0,             0,
                R(1,0),  R(1,1),  R(1,2),  0,             0,             0,
                R(2,0),  R(2,1),  R(2,2),  0,             0,             0,
                pR(0,0), pR(0,1), pR(0,2), R(0,0), R(0,1), R(0,2),
                pR(1,0), pR(1,1), pR(1,2), R(1,0), R(1,1), R(1,2),
                pR(2,0), pR(2,1), pR(2,2), R(2,0), R(2,1), R(2,2);

    return adjoint;
}

double math::cot(double x){
    return 1/tan(x);
}

void math::wrench_in_s_and_w(){
    const Eigen::Vector3d f_w {-30,0,0};
    const Eigen::Vector3d m_s {0,0,2};
    const Eigen::Vector3d e_ws {60,-60,0}; // Euler angles YZX in degrees.
    Eigen::Matrix3d R_ws {rotation_matrix_from_euler_yzx(e_ws)};

    Eigen::Vector3d m_w {R_ws*m_s};
    Eigen::Vector3d f_s {R_ws.transpose()*f_w};


    //std::cout << "Still in progress" << std::endl;
    std::cout << "f_w: " << f_w.transpose() << std::endl;
    std::cout << "m_w: " << m_w.transpose() << std::endl;
    std::cout << "f_s: " << f_s.transpose() << std::endl;
    std::cout << "m_s: " << m_s.transpose() << std::endl;
    std::cout << "R_ws: " << std::endl << R_ws << std::endl << std::endl;
}

Eigen::VectorXd math::wrench_w(const Eigen::Vector3d &f_w, const Eigen::Vector3d &m_s, const Eigen::Vector3d &e_ws) {
    Eigen::Matrix3d R_ws {rotation_matrix_from_euler_yzx(e_ws)};
    Eigen::Vector3d m_w {R_ws*m_s};
    //Eigen::Vector3d f_s {R_ws.transpose()*f_w};

    Eigen::VectorXd F_w(6);
        F_w << m_w(0), m_w(1), m_w(2), f_w(0), f_w(1), f_w(2);

    return F_w;
}

Eigen::VectorXd math::wrench_s(const Eigen::Vector3d &f_w, const Eigen::Vector3d &m_s, const Eigen::Vector3d &e_ws) {
    Eigen::Matrix3d R_ws {rotation_matrix_from_euler_yzx(e_ws)};
    //Eigen::Vector3d m_w {R_ws*m_s};
    Eigen::Vector3d f_s {R_ws.transpose()*f_w};

    Eigen::VectorXd F_s(6);
    F_s << m_s(0), m_s(1), m_s(2), f_s(0), f_s(1), f_s(2);

    return F_s;
}

Eigen::VectorXd math::wrench_a_to_b(const Eigen::VectorXd &F_a, const Eigen::Matrix4d &tf) {
    Eigen::MatrixXd Ad_T(6,6);
    Ad_T = {math::adjoint_matrix(tf)};

    return Ad_T.transpose()*F_a;
}

// Task 2b
Eigen::VectorXd math::wrench_f(const Eigen::VectorXd &F_a, const Eigen::VectorXd &F_b, const Eigen::MatrixXd &tf_af, const Eigen::MatrixXd &tf_bf) {
    return wrench_a_to_b(F_a, tf_af) + wrench_a_to_b(F_b, tf_bf);
}

// Task 3a
Eigen::Matrix3d math::matrix_exponential(const Eigen::Vector3d &w, const double theta) {
    const Eigen::Matrix3d skew_w {skew_symmetric(w)};
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = theta*M_PI/180;

    return I + sin(rads)*skew_w + (1-cos(rads))*skew_w*skew_w;
}

// Task 3b
std::pair<Eigen::Vector3d, double> math::matrix_logarithm(const Eigen::Matrix3d &r) {
    double theta {};
    double w_1 {};
    double w_2 {};
    double w_3 {};
    Eigen::Vector3d w;
    const double trR = r(0,0)+r(1,1)+r(2,2);

    if(r == Eigen::Matrix3d::Identity()) {
        theta = 0;
        w = {0,0,0};
    }
    else if(trR == -1) {
        theta = EIGEN_PI;
        w_1 = r(0,2)/sqrt(2*(1+r(2,2)));
        w_2 = r(1,2)/sqrt(2*(1+r(2,2)));
        w_3 = 1 + r(2,2)/sqrt(2*(1+r(2,2)));
        w = {w_1,w_2,w_3};
    }
    else {
        theta = acos((trR-1)/2);
        Eigen::Matrix3d skew_w = (r-r.transpose())/(2*sin(theta));
        w = {skew_w(2,1), skew_w(0,2), skew_w(1,0)};
    }


    return std::make_pair(w, theta*180/EIGEN_PI);
}

//Task 3c. se(3) -> SE(3)
Eigen::Matrix4d math::matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, const double theta) {
    const Eigen::Matrix3d skew_w {skew_symmetric(w)};
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = theta*M_PI/180;

    Eigen::Matrix3d R = matrix_exponential(w, theta);

    Eigen::Vector3d p = (I*rads + (1-cos(rads))*skew_w + (rads - sin(rads))*skew_w*skew_w)*v;

    Eigen::Matrix4d T ;
    T <<    R(0,0), R(0,1), R(0,2), p(0),
            R(1,0), R(1,1), R(1,2), p(1),
            R(2,0), R(2,1), R(2,2), p(2),
                        0,             0,              0,        1;

    return T;
}

// Task 3d. SE(3) -> se(3)
Eigen::Matrix3d  math::G(const Eigen::Vector3d &w, const double &theta) {
    Eigen::Matrix3d skew_w {skew_symmetric(w)};
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = theta*M_PI/180;

    return I*rads + (1-cos(rads))*skew_w + (rads- sin(rads))*skew_w*skew_w;
}

Eigen::Matrix3d math::G_inverse(const Eigen::Vector3d &w, const double &degrees) {
    Eigen::Matrix3d skew_w {skew_symmetric(w)};
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    const double rads = degrees*M_PI/180;

    return I/rads - skew_w/2 + (1/rads-cot(rads/2)/2)*skew_w*skew_w;
}

// Working
std::pair<Eigen::VectorXd, double> math::matrix_logarithm(const Eigen::Matrix4d &t) {
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double degrees;

    Eigen::Matrix3d R;
    R << t(0,0), t(0,1), t(0,2),
        t(1,0), t(1,1), t(1,2),
        t(2,0), t(2,1), t(2,2);
    Eigen::Vector3d p {t(0,3), t(1,3), t(2,3)};

    if(R == Eigen::Matrix3d::Identity()) {
        w = {0,0,0};
        v = p.normalized();
        degrees = sqrt(p(0)*p(0) + p(1)*p(1) + p(2)*p(2));
    }

    else {
        auto [fst, scd] = math::matrix_logarithm(R);
        w = fst;
        degrees = scd;

        const Eigen::Matrix3d G_inv = math::G_inverse(w,degrees);
        v = G_inv*p;
    }
    Eigen::VectorXd S(6);
    S << w(0), w(1), w(2), v(0), v(1), v(2);

    return std::make_pair(S, degrees);
}

void math::print_pose(const std::string &label, const Eigen::Matrix4d &tf){
    Eigen::Matrix3d R;
    R << tf(0,0), tf(0,1), tf(0,2),
         tf(1,0), tf(1,1), tf(1,2),
         tf(2,0), tf(2,1), tf(2,2);
    Eigen::Vector3d p {tf(0,3), tf(1,3), tf(2,3)};

    Eigen::Vector3d e_zyx = euler_zyx_from_rotation(R);

    std::cout << "Label: " << label << std::endl;
    std::cout << "Euler ZYX angles: " << e_zyx.transpose()*180/EIGEN_PI << std::endl;
    std::cout << "Linear position: " << p.transpose() << std::endl;
    std::cout << " " << std::endl;
}

// Task 4b
Eigen::Matrix4d math::planar_3r_fk_transform(const std::vector<double> &joint_positions) {
    constexpr double L1 = 10, L2 = 10, L3 = 10;

    // Making transformation matrix for each succeeding frame. Every joint rotates around z-axis.
    const Eigen::Matrix4d T01 = transformation_matrix(rotate_z(joint_positions[0]), Eigen::Vector3d(0, 0, 0));
    const Eigen::Matrix4d T12 = transformation_matrix(rotate_z(joint_positions[1]), Eigen::Vector3d(L1, 0, 0));
    const Eigen::Matrix4d T23 = transformation_matrix(rotate_z(joint_positions[2]), Eigen::Vector3d(L2, 0, 0));
    const Eigen::Matrix4d T34 = transformation_matrix(rotate_z(0), Eigen::Vector3d(L3, 0, 0));

    Eigen::Matrix4d T04 = T01 * T12 * T23 * T34;

    return T04;
}

void math::test_planar_3r_fk_transform(const std::string &label, const std::vector<double> &joint_positions) {
    const Eigen::Matrix4d T {math::planar_3r_fk_transform(joint_positions)};

    print_pose(label, T);
}

// Task 4c
Eigen::Matrix4d math::planar_3r_fk_screw(const std::vector<double> &joint_positions) {
    constexpr double L1 = 10, L2 = 10, L3 = 10;

    Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
    M(0,3) = L1+L2+L3;

    Eigen::Vector3d w1, w2, w3;
    w1 << 0,0,1;
    w2 << 0,0,1;
    w3 << 0,0,1;

    Eigen::Vector3d v1, v2, v3;
    v1 << 0,0,0;
    v2 << 0,-L1,0;
    v3 << 0,-L1-L2,0;

    const Eigen::Matrix4d e1 = matrix_exponential(w1,v1,joint_positions[0]);
    const Eigen::Matrix4d e2 = matrix_exponential(w2,v2,joint_positions[1]);
    const Eigen::Matrix4d e3 = matrix_exponential(w3,v3,joint_positions[2]);

    const Eigen::Matrix4d T04 = e1 * e2 * e3 * M;

    return T04;
}

void math::test_planar_3r_fk_screw(const std::string &label, const std::vector<double> &joint_positions) {
    const Eigen::Matrix4d T {planar_3r_fk_screw(joint_positions)};

    print_pose(label, T);
}


// TASK 5
// Task 5a
Eigen::Matrix4d math::ur3e_fk_screw(const std::vector<double> &joint_positions) {  //OK
    constexpr double h1 {0.15185}, l1 {0.24355}, l2 {0.2132}, h2 {0.08535},
    y1{0.13105}, y2{0.0921};

    Eigen::Matrix4d M; //ok
    M <<  1,  0, 0, -l1-l2,
          0,  0,-1, -y1-y2,
          0,  1, 0, h1-h2,
          0,  0, 0, 1;

    Eigen::Vector3d w0, w1, w2, w3, w4, w5;
    w0 << 0,0,1; //(base) //ok
    w1 << 0,-1,0; //ok
    w2 << 0,-1,0; //ok
    w3 << 0,-1,0; //ok
    w4 << 0,0,-1;  //ok
    w5 << 0,-1,0; //ok

    Eigen::Vector3d v0, v1, v2, v3, v4, v5;
    v0 << 0,0,0; //ok
    v1 << h1,0,0;  //ok
    v2 << h1,0,l1; //ok
    v3 << h1,0,l1+l2; //ok
    v4 << y1,-l1-l2,0; //ok
    v5 << h1-h2,0,l1+l2; //ok


    const Eigen::Matrix4d e0 = matrix_exponential(w0,v0,joint_positions[0]);
    const Eigen::Matrix4d e1 = matrix_exponential(w1,v1,joint_positions[1]);
    const Eigen::Matrix4d e2 = matrix_exponential(w2,v2,joint_positions[2]);
    const Eigen::Matrix4d e3 = matrix_exponential(w3,v3,joint_positions[3]);
    const Eigen::Matrix4d e4 = matrix_exponential(w4,v4,joint_positions[4]);
    const Eigen::Matrix4d e5 = matrix_exponential(w5,v5,joint_positions[5]);

    const Eigen::Matrix4d T = e0 * e1 * e2 * e3 * e4 * e5 * M;

    return T;
}


void math::test_ur3e_fk_screw(const std::string &label, const std::vector<double> &joint_positions) {
    const Eigen::Matrix4d T {ur3e_fk_screw(joint_positions)};

    print_pose(label, T);
}

Eigen::Matrix4d DH_transformation_matrix(const double &joint_angle, const double &alpha, const double &a, const double &d) {
    Eigen::Matrix4d T_mn;
    double joint_rads = joint_angle*M_PI/180;
    double alpha_rads = alpha*M_PI/180;
    T_mn << cos(joint_rads), -sin(joint_rads)*cos(alpha_rads),  sin(joint_rads)*sin(alpha_rads), a*cos(joint_rads),  //ok
            sin(joint_rads),  cos(joint_rads)*cos(alpha_rads), -cos(joint_rads)*sin(alpha_rads), a*sin(joint_rads),
                          0,                  sin(alpha_rads),                  cos(alpha_rads),                 d,
                          0,                                0,                                0,                 1;
    return T_mn;
}

Eigen::Matrix4d math::ur3e_fk_transform(const std::vector<double> &joint_positions) { //YES!
    constexpr double h0 {0.15185}, h1 {0.24355}, h2 {0.2132}, h3 {0.08535},
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
    const Eigen::Matrix4d T {ur3e_fk_transform(joint_positions)};

    print_pose(label, T);
}