#include "parameter_validator.h"
#include "curve_generator.h"
#include <iostream>

int main() {
    std::cout << "Testing Montgomery curve parameter validation..." << std::endl;

    // 生成一个曲线用于测试
    MontgomeryCurve curve;
    init_montgomery_curve(&curve);

    if (!generate_montgomery_curve(&curve, 112)) {
        std::cout << "Failed to generate test curve!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    std::cout << "Generated curve for validation test:" << std::endl;
    std::cout << "Prime p: ";
    mpz_out_str(stdout, 10, curve.p);
    std::cout << std::endl;

    // 验证曲线参数
    if (validate_montgomery_curve(&curve)) {
        std::cout << "Curve validation passed!" << std::endl;
    } else {
        std::cout << "Curve validation failed!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    // 测试无效参数的情况
    std::cout << "\nTesting invalid curve parameters..." << std::endl;
    mpz_t invalid_A;
    mpz_init_set_ui(invalid_A, 2); // A = 2 是无效的

    if (!validate_curve_parameters(invalid_A, curve.B, curve.p)) {
        std::cout << "Correctly detected invalid parameter A = 2" << std::endl;
    } else {
        std::cout << "Failed to detect invalid parameter A = 2" << std::endl;
        mpz_clear(invalid_A);
        clear_montgomery_curve(&curve);
        return 1;
    }

    mpz_clear(invalid_A);
    clear_montgomery_curve(&curve);
    std::cout << "Validation test passed!" << std::endl;
    return 0;
}
