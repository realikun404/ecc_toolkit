#include "curve_generator.h"
#include <iostream>

int main() {
    std::cout << "Testing Montgomery curve parameter generation..." << std::endl;

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
    } else {
        std::cout << "Failed to generate Montgomery curve!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    clear_montgomery_curve(&curve);
    std::cout << "Test passed!" << std::endl;
    return 0;
}
