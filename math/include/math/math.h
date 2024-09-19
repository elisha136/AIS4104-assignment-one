
#ifndef MATH_H
#define MATH_H
#include <Eigen/Dense>

namespace math {
    bool floatEquals(double a, double b);
    double deg_to_rad(double degrees);
    Eigen::Matrix3d rotate_x(double degrees);
    Eigen::Matrix3d rotate_y(double degrees);
    Eigen::Matrix3d rotate_z(double degrees);
    Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v);
    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,const Eigen::Vector3d &y,const Eigen::Vector3d &z);
    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees);
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e);
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);
    Eigen::Vector3d euler_zyx_from_rotation(const Eigen::Matrix3d &r);
    //Eigen::Vector3d euler_zyx_from_transformation(const Eigen::Matrix4d &t);
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);
    Eigen::VectorXd wrench_w(const Eigen::Vector3d &f_w, const Eigen::Vector3d &m_s, const Eigen::Vector3d &e_ws);
    Eigen::VectorXd wrench_s(const Eigen::Vector3d &f_w, const Eigen::Vector3d &m_s, const Eigen::Vector3d &e_ws);
    Eigen::VectorXd wrench_a_to_b(const Eigen::VectorXd &F_a, const Eigen::Matrix4d &tf);
    Eigen::VectorXd wrench_f(const Eigen::VectorXd &F_a, const Eigen::VectorXd &F_b, const Eigen::MatrixXd &tf_af, const Eigen::MatrixXd &tf_bf);
    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta); // so(3) -> SO(3)
    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r);  // SO(3) - so(3)
    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta); // se(3) -> SE(3)
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t); // SE(3) -> se(3)
    Eigen::Matrix3d  G(const Eigen::Vector3d &w, const double &theta);
    Eigen::Matrix3d  G_inverse(const Eigen::Vector3d &w, const  double &degrees);
    Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions);
    Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions);
    Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions);
    Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions);


    double cot(double x);
    void wrench_in_s_and_w();
    void print_pose(const std::string &label, const Eigen::Matrix4d &tf);
    void test_planar_3r_fk_transform(const std::string &label, const std::vector<double> &joint_positions);
    void test_planar_3r_fk_screw(const std::string &label, const std::vector<double> &joint_positions);
    void test_ur3e_fk_screw(const std::string &label, const std::vector<double> &joint_positions);
    void test_ur3e_fk_transform(const std::string &label, const std::vector<double> &joint_positions);
    constexpr double c_rad_to_deg{57.2957795};
    constexpr double c_deg_to_rad{0.01745329251};
}

#endif //MATH_H