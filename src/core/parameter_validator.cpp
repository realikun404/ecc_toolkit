#include "parameter_validator.h"
#include <gmp.h>
#include <iostream>

// 验证Montgomery曲线参数的合法性
bool validate_montgomery_curve(const MontgomeryCurve* curve) {
    // 验证素数域p
    if (!validate_prime_field(curve->p)) {
        std::cerr << "Invalid prime field" << std::endl;
        return false;
    }

    // 验证曲线参数A和B
    if (!validate_curve_parameters(curve->A, curve->B, curve->p)) {
        std::cerr << "Invalid curve parameters" << std::endl;
        return false;
    }

    // 验证基点是否在曲线上（如果基点已设置）
    if (mpz_sgn(curve->x) != 0 || mpz_sgn(curve->y) != 0) {
        if (!validate_base_point(curve->x, curve->y, curve->A, curve->B, curve->p)) {
            std::cerr << "Base point is not on the curve" << std::endl;
            return false;
        }
    }

    return true;
}

// 验证素数域p是否为有效素数
bool validate_prime_field(const mpz_t p) {
    // 检查p是否大于3
    if (mpz_cmp_ui(p, 3) <= 0) {
        return false;
    }

    // 检查p是否为素数
    return mpz_probab_prime_p(p, 25) > 0; // 25轮Miller-Rabin测试
}

// 验证曲线参数A和B是否满足Montgomery曲线条件
bool validate_curve_parameters(const mpz_t A, const mpz_t B, const mpz_t p) {
    // B不能为0 (mod p)
    if (mpz_divisible_p(B, p)) {
        return false;
    }

    // A不能等于±2 (mod p)
    mpz_t two, minus_two;
    mpz_init_set_ui(two, 2);
    mpz_init_set_si(minus_two, -2);

    // 计算A mod p
    mpz_t A_mod_p;
    mpz_init(A_mod_p);
    mpz_mod(A_mod_p, A, p);

    bool result = true;
    if (mpz_congruent_p(A_mod_p, two, p) || mpz_congruent_p(A_mod_p, minus_two, p)) {
        result = false;
    }

    mpz_clear(two);
    mpz_clear(minus_two);
    mpz_clear(A_mod_p);

    return result;
}

// 验证基点是否在曲线上: By^2 = x^3 + Ax^2 + x
bool validate_base_point(const mpz_t x, const mpz_t y, const mpz_t A, const mpz_t B, const mpz_t p) {
    // 计算左边: By^2
    mpz_t left, temp1, temp2;
    mpz_init(left);
    mpz_init(temp1);
    mpz_init(temp2);

    // left = B * y^2 mod p
    mpz_mul(temp1, y, y);        // y^2
    mpz_mul(left, B, temp1);     // B * y^2
    mpz_mod(left, left, p);      // (B * y^2) mod p

    // 计算右边: x^3 + Ax^2 + x
    mpz_t right;
    mpz_init(right);

    // x^2
    mpz_mul(temp1, x, x);
    mpz_mod(temp1, temp1, p);

    // Ax^2
    mpz_mul(temp2, A, temp1);
    mpz_mod(temp2, temp2, p);

    // x^3
    mpz_mul(temp1, temp1, x);
    mpz_mod(temp1, temp1, p);

    // x^3 + Ax^2
    mpz_add(right, temp1, temp2);
    mpz_mod(right, right, p);

    // x^3 + Ax^2 + x
    mpz_add(right, right, x);
    mpz_mod(right, right, p);

    // 比较左右两边是否相等
    bool result = (mpz_cmp(left, right) == 0);

    mpz_clear(left);
    mpz_clear(right);
    mpz_clear(temp1);
    mpz_clear(temp2);

    return result;
}
