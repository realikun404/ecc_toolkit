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

            if (is_point_on_curve(&base_point, curve.A, curve.B, curve.p)) {
                std::cout << "SUCCESS: Base point is on the curve!" << std::endl;
            } else {
                std::cout << "ERROR: Base point is NOT on the curve!" << std::endl;
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
