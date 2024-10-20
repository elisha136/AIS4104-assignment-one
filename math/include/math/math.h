
#ifndef MATH_H
#define MATH_H
#include <Eigen/Dense>
#include <functional>

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
    Eigen::Matrix4d DH_transformation_matrix(const double &joint_angle, const double &alpha, const double &a,
                                         const double &d);
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
    //Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions);
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
    Eigen::VectorXd std_vector_to_eigen1(const std::vector<double> &v);
    Eigen::VectorXd std_vector_to_eigen(const std::vector<double> &v);
    //Assignment 3
    bool is_average_below_eps(const std::vector<double> &values, double eps = 10e-7, uint8_t n_values = 5u);
    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd> > ur3e_space_chain();
    Eigen::Matrix4d ur3e_space_fk(const Eigen::VectorXd &joint_positions);
    Eigen::Matrix4d ur3e_body_fk(const Eigen::VectorXd &joint_positions);
    std::pair<uint32_t, double> newton_raphson_root_find(const std::function<double(double)> &f, double x_0,
                                                     double dx_0 = 0.5, double eps = 1e-7, uint32_t max_iters = 1000);
    std::pair<uint32_t, double> gradient_descent_root_find(const std::function<double(double)> &f, double x_0,
                                                       double gamma = 0.1, double dx_0 = 0.5, double eps = 1e-7,
                                                       uint32_t max_iters = 1000);
    void test_newton_raphson_root_find(const std::function<double(double)> &f, double x0);
    void  test_gradient_descent_root_find(const std::function<double(double)> &f, double x0);
    void test_root_find();
    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_body_chain();
    Eigen::MatrixXd  ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions);

    Eigen::MatrixXd ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions);
    void ur3e_test_jacobian(const Eigen::VectorXd &joint_positions);
    void ur3e_test_jacobian();
  std::pair<size_t, Eigen::VectorXd> ur3e_ik_body(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &current_joint_positions,
                                           double gamma = 1e-2, double v_e = 4e-3, double w_e = 4e-3);

    void ur3e_ik_test();
    void ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0);

    void ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd &j0);


    //
}

#endif //MATH_H