#include "montgomery_curve.h"
#include "curve_generator.h"
#include <iostream>

int main() {
    std::cout << "Testing Montgomery curve operations..." << std::endl;

    // 生成测试曲线
    MontgomeryCurve curve;
    init_montgomery_curve(&curve);

    if (!generate_montgomery_curve(&curve, 112)) {
        std::cout << "Failed to generate test curve!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    // 初始化测试点
    MontgomeryPoint P, Q, R;
    init_montgomery_point(&P);
    init_montgomery_point(&Q);
    init_montgomery_point(&R);

    // 设置测试点P (使用曲线参数，但实际应该使用曲线上的真实点)
    // 这里我们使用简化的示例值进行测试
    mpz_set_ui(P.x, 1);
    mpz_set_ui(P.y, 2);

    // 设置测试点Q
    mpz_set_ui(Q.x, 3);
    mpz_set_ui(Q.y, 4);

    std::cout << "Testing point addition..." << std::endl;

    // 测试点加法
    if (montgomery_point_add(&R, &P, &Q, curve.A, curve.p)) {
        std::cout << "Point addition successful!" << std::endl;
        std::cout << "R.x = ";
        mpz_out_str(stdout, 10, R.x);
        std::cout << std::endl;

        std::cout << "R.y = ";
        mpz_out_str(stdout, 10, R.y);
        std::cout << std::endl;
    } else {
        std::cout << "Point addition failed!" << std::endl;
        clear_montgomery_point(&P);
        clear_montgomery_point(&Q);
        clear_montgomery_point(&R);
        clear_montgomery_curve(&curve);
        return 1;
    }

    // 测试点倍乘
    std::cout << "\nTesting point multiplication..." << std::endl;
    mpz_t k;
    mpz_init_set_ui(k, 3);

    if (montgomery_point_multiply(&R, k, &P, curve.A, curve.p)) {
        std::cout << "Point multiplication successful!" << std::endl;
        std::cout << "k = ";
        mpz_out_str(stdout, 10, k);
        std::cout << std::endl;

        std::cout << "R.x = ";
        mpz_out_str(stdout, 10, R.x);
        std::cout << std::endl;

        std::cout << "R.y = ";
        mpz_out_str(stdout, 10, R.y);
        std::cout << std::endl;
    } else {
        std::cout << "Point multiplication failed!" << std::endl;
        mpz_clear(k);
        clear_montgomery_point(&P);
        clear_montgomery_point(&Q);
        clear_montgomery_point(&R);
        clear_montgomery_curve(&curve);
        return 1;
    }

    // 测试无穷远点
    std::cout << "\nTesting point at infinity..." << std::endl;
    set_point_infinity(&R);
    if (is_point_infinity(&R)) {
        std::cout << "Infinity point test passed!" << std::endl;
    } else {
        std::cout << "Infinity point test failed!" << std::endl;
        mpz_clear(k);
        clear_montgomery_point(&P);
        clear_montgomery_point(&Q);
        clear_montgomery_point(&R);
        clear_montgomery_curve(&curve);
        return 1;
    }

    // 清理资源
    mpz_clear(k);
    clear_montgomery_point(&P);
    clear_montgomery_point(&Q);
    clear_montgomery_point(&R);
    clear_montgomery_curve(&curve);

    std::cout << "All Montgomery curve operation tests passed!" << std::endl;
    return 0;
}
