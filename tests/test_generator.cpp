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

            clear_montgomery_point(&base_point);
        } else {
            std::cout << "WARNING: Base point was not found!" << std::endl;
        }
    } else {
        std::cout << "Failed to generate Montgomery curve!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    clear_montgomery_curve(&curve);
    std::cout << "Base point generation test completed successfully!" << std::endl;
    return 0;
}
