// test_validator.cpp
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

    // 测试无效的B参数
    mpz_t invalid_B;
    mpz_init(invalid_B);

    if (!validate_curve_parameters(curve.A, invalid_B, curve.p)) {
        std::cout << "Correctly detected invalid parameter B = 0" << std::endl;
    } else {
        std::cout << "Failed to detect invalid parameter B = 0" << std::endl;
        mpz_clear(invalid_B);
        clear_montgomery_curve(&curve);
        return 1;
    }

    mpz_clear(invalid_B);

    // 测试无效的素数域
    mpz_t invalid_p;
    mpz_init_set_ui(invalid_p, 2); // p = 2 太小

    if (!validate_prime_field(invalid_p)) {
        std::cout << "Correctly detected invalid prime field p = 2" << std::endl;
    } else {
        std::cout << "Failed to detect invalid prime field p = 2" << std::endl;
        mpz_clear(invalid_p);
        clear_montgomery_curve(&curve);
        return 1;
    }

    mpz_clear(invalid_p);

    // 测试基点验证
    std::cout << "\nTesting base point validation..." << std::endl;
    if (validate_base_point(curve.x, curve.y, curve.A, curve.B, curve.p)) {
        std::cout << "Base point validation passed!" << std::endl;
    } else {
        std::cout << "Base point validation failed!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    // 测试无效基点
    mpz_t invalid_x, invalid_y;
    mpz_init_set_ui(invalid_x, 1);
    mpz_init_set_ui(invalid_y, 1);

    if (!validate_base_point(invalid_x, invalid_y, curve.A, curve.B, curve.p)) {
        std::cout << "Correctly detected invalid base point (1,1)" << std::endl;
    } else {
        std::cout << "Failed to detect invalid base point (1,1)" << std::endl;
        mpz_clear(invalid_x);
        mpz_clear(invalid_y);
        clear_montgomery_curve(&curve);
        return 1;
    }

    mpz_clear(invalid_x);
    mpz_clear(invalid_y);

    clear_montgomery_curve(&curve);
    std::cout << "All validation tests passed!" << std::endl;
    return 0;
}
