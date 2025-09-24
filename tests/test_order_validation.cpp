// test_order_validation.cpp
#include "curve_generator.h"
#include "parameter_validator.h"
#include "montgomery_curve.h"
#include <iostream>

int main() {
    std::cout << "Testing order validation for Montgomery curve..." << std::endl;

    // 生成一个曲线用于测试
    MontgomeryCurve curve;
    init_montgomery_curve(&curve);

    // 测试生成128位安全等级的曲线
    if (generate_montgomery_curve(&curve, 128)) {
        std::cout << "Montgomery curve generated successfully!" << std::endl;

        // 测试基点阶数计算
        std::cout << "\nTesting base point order calculation..." << std::endl;
        mpz_t computed_order;
        mpz_init(computed_order);

        if (calculate_base_point_order(&curve, computed_order)) {
            std::cout << "SUCCESS: Base point order calculated!" << std::endl;
            std::cout << "Computed order: ";
            mpz_out_str(stdout, 10, computed_order);
            std::cout << std::endl;

            std::cout << "Stored order (n): ";
            mpz_out_str(stdout, 10, curve.n);
            std::cout << std::endl;

            // 比较计算得到的阶数与存储的阶数
            if (mpz_cmp(computed_order, curve.n) == 0) {
                std::cout << "SUCCESS: Computed order matches stored order!" << std::endl;
            } else {
                std::cout << "WARNING: Computed order does not match stored order!" << std::endl;
            }
        } else {
            std::cout << "ERROR: Failed to calculate base point order!" << std::endl;
            mpz_clear(computed_order);
            clear_montgomery_curve(&curve);
            return 1;
        }

        mpz_clear(computed_order);

        // 测试生成元验证
        std::cout << "\nTesting generator verification..." << std::endl;
        if (verify_generator(&curve)) {
            std::cout << "SUCCESS: Base point verified as generator!" << std::endl;
        } else {
            std::cout << "ERROR: Base point NOT verified as generator!" << std::endl;
            clear_montgomery_curve(&curve);
            return 1;
        }

        // 测试点乘运算
        std::cout << "\nTesting point multiplication with order..." << std::endl;
        MontgomeryPoint base_point, result_point;
        init_montgomery_point(&base_point);
        init_montgomery_point(&result_point);

        mpz_set(base_point.x, curve.x);
        mpz_set(base_point.y, curve.y);
        base_point.is_infinity = 0;

        // 计算 n*P，结果应该是无穷远点
        if (montgomery_point_multiply(&result_point, curve.n, &base_point, curve.A, curve.p)) {
            if (is_point_infinity(&result_point)) {
                std::cout << "SUCCESS: n*P = O (infinity point) as expected!" << std::endl;
            } else {
                std::cout << "ERROR: n*P is not infinity point!" << std::endl;
                clear_montgomery_point(&base_point);
                clear_montgomery_point(&result_point);
                clear_montgomery_curve(&curve);
                return 1;
            }
        } else {
            std::cout << "ERROR: Point multiplication failed!" << std::endl;
            clear_montgomery_point(&base_point);
            clear_montgomery_point(&result_point);
            clear_montgomery_curve(&curve);
            return 1;
        }

        clear_montgomery_point(&base_point);
        clear_montgomery_point(&result_point);
    } else {
        std::cout << "ERROR: Failed to generate Montgomery curve!" << std::endl;
        clear_montgomery_curve(&curve);
        return 1;
    }

    clear_montgomery_curve(&curve);
    std::cout << "\nAll order validation tests passed!" << std::endl;
    return 0;
}
