// test_generator.cpp
#include "curve_generator.h"
#include "parameter_validator.h"
#include "montgomery_curve.h"
#include <iostream>

int main() {
    std::cout << "Testing base point generation for Montgomery curve..." << std::endl;

    // 生成一个曲线用于测试
    MontgomeryCurve curve;
    init_montgomery_curve(&curve);

    // 测试生成128位安全等级的曲线
    if (generate_montgomery_curve(&curve, 128)) {
        std::cout << "Montgomery curve generated successfully!" << std::endl;
        std::cout << "Security level: " << curve.security_level << std::endl;

        std::cout << "Prime p: ";
        mpz_out_str(stdout, 10, curve.p);
        std::cout << std::endl;

        std::cout << "Parameter A: ";
        mpz_out_str(stdout, 10, curve.A);
        std::cout << std::endl;

        std::cout << "Parameter B: ";
        mpz_out_str(stdout, 10, curve.B);
        std::cout << std::endl;

        std::cout << "Base point x: ";
        mpz_out_str(stdout, 10, curve.x);
        std::cout << std::endl;

        std::cout << "Base point y: ";
        mpz_out_str(stdout, 10, curve.y);
        std::cout << std::endl;

        std::cout << "Base point order n: ";
        mpz_out_str(stdout, 10, curve.n);
        std::cout << std::endl;

        // 验证基点是否在曲线上
        if (mpz_sgn(curve.x) != 0 || mpz_sgn(curve.y) != 0) {
            MontgomeryPoint base_point;
            init_montgomery_point(&base_point);
            mpz_set(base_point.x, curve.x);
            mpz_set(base_point.y, curve.y);

            // 手动验证点是否在曲线上
            mpz_t left, right, temp1, temp2;
            mpz_init(left);
            mpz_init(right);
            mpz_init(temp1);
            mpz_init(temp2);

            // 计算左边: By^2
            mpz_mul(temp1, curve.y, curve.y);      // y^2
            mpz_mod(temp1, temp1, curve.p);
            mpz_mul(left, curve.B, temp1);         // B*y^2
            mpz_mod(left, left, curve.p);

            // 计算右边: x^3 + Ax^2 + x
            mpz_mul(temp1, curve.x, curve.x);      // x^2
            mpz_mod(temp1, temp1, curve.p);

            mpz_mul(temp2, curve.A, temp1);        // A*x^2
            mpz_mod(temp2, temp2, curve.p);

            mpz_mul(temp1, temp1, curve.x);        // x^3
            mpz_mod(temp1, temp1, curve.p);

            mpz_add(right, temp1, temp2);          // x^3 + A*x^2
            mpz_mod(right, right, curve.p);

            mpz_add(temp1, right, curve.x);        // x^3 + A*x^2 + x
            mpz_mod(temp1, temp1, curve.p);

            std::cout << "Left side (By^2): ";
            mpz_out_str(stdout, 10, left);
            std::cout << std::endl;

            std::cout << "Right side (x^3 + Ax^2 + x): ";
            mpz_out_str(stdout, 10, temp1);
            std::cout << std::endl;

            if (mpz_cmp(left, temp1) == 0) {
                std::cout << "SUCCESS: Base point is on the curve!" << std::endl;
            } else {
                std::cout << "ERROR: Base point is NOT on the curve!" << std::endl;
                mpz_clear(left);
                mpz_clear(right);
                mpz_clear(temp1);
                mpz_clear(temp2);
                clear_montgomery_point(&base_point);
                clear_montgomery_curve(&curve);
                return 1;
            }

            mpz_clear(left);
            mpz_clear(right);
            mpz_clear(temp1);
            mpz_clear(temp2);

            // 使用现有的函数验证
            if (is_point_on_curve(&base_point, curve.A, curve.B, curve.p)) {
                std::cout << "SUCCESS: Base point validated with is_point_on_curve function!" << std::endl;
            } else {
                std::cout << "ERROR: Base point not validated with is_point_on_curve function!" << std::endl;
                clear_montgomery_point(&base_point);
                clear_montgomery_curve(&curve);
                return 1;
            }

            // 验证基点是否为生成元
            if (verify_generator(&curve)) {
                std::cout << "SUCCESS: Base point is verified as a generator!" << std::endl;
            } else {
                std::cout << "ERROR: Base point is NOT a generator!" << std::endl;
                clear_montgomery_point(&base_point);
                clear_montgomery_curve(&curve);
                return 1;
            }

            // 验证阶数是否为素数
            if (mpz_probab_prime_p(curve.n, 25)) {
                std::cout << "SUCCESS: Base point order is prime!" << std::endl;
            } else {
                std::cout << "WARNING: Base point order is not prime!" << std::endl;
            }

            clear_montgomery_point(&base_point);
        } else {
            std::cout << "ERROR: Base point coordinates are both zero!" << std::endl;
            clear_montgomery_curve(&curve);
            return 1;
        }
    } else {
        std::cout << "Failed to generate Montgomery curve with valid base point!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    clear_montgomery_curve(&curve);
    std::cout << "Base point generation test completed successfully!" << std::endl;

    // 额外测试不同安全级别的曲线生成
    std::cout << "\nTesting different security levels..." << std::endl;
    int security_levels[] = {112, 128, 192, 256};
    int num_levels = sizeof(security_levels) / sizeof(security_levels[0]);

    for (int i = 0; i < num_levels; i++) {
        MontgomeryCurve test_curve;
        init_montgomery_curve(&test_curve);

        std::cout << "Testing security level " << security_levels[i] << "..." << std::endl;
        if (generate_montgomery_curve(&test_curve, security_levels[i])) {
            std::cout << "SUCCESS: Curve generated for security level " << security_levels[i] << std::endl;

            // 验证曲线参数
            if (validate_montgomery_curve(&test_curve)) {
                std::cout << "SUCCESS: Curve parameters validated for security level " << security_levels[i] << std::endl;
            } else {
                std::cout << "ERROR: Curve parameters validation failed for security level " << security_levels[i] << std::endl;
                clear_montgomery_curve(&test_curve);
                return 1;
            }

            // 验证基点是否为生成元
            if (verify_generator(&test_curve)) {
                std::cout << "SUCCESS: Base point is verified as a generator for security level " << security_levels[i] << std::endl;
            } else {
                std::cout << "ERROR: Base point is NOT a generator for security level " << security_levels[i] << std::endl;
                clear_montgomery_curve(&test_curve);
                return 1;
            }
        } else {
            std::cout << "ERROR: Failed to generate curve for security level " << security_levels[i] << std::endl;
            clear_montgomery_curve(&test_curve);
            return 1;
        }

        clear_montgomery_curve(&test_curve);
    }

    std::cout << "All security level tests passed!" << std::endl;
    return 0;
}
