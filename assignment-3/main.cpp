

#include <iostream>
#include <Eigen/Dense>
#include "math/math.h"
#include <functional>

int main() {
    math::test_root_find();
    math::ur3e_test_jacobian();
    math::ur3e_ik_test();

    return 0;
}
