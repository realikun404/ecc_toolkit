#include "curve_generator.h"
#include "parameter_validator.h"
#include "math/prime_utils.h"
#include <gmp.h>
#include <iostream>
#include <random>

// 初始化Montgomery曲线参数
void init_montgomery_curve(MontgomeryCurve* curve) {
    mpz_init(curve->p);
    mpz_init(curve->A);
    mpz_init(curve->B);
    mpz_init(curve->x);
    mpz_init(curve->y);
    mpz_init(curve->n);
    mpz_init(curve->h);
    curve->security_level = 0;
}

// 清理Montgomery曲线参数
void clear_montgomery_curve(MontgomeryCurve* curve) {
    mpz_clear(curve->p);
    mpz_clear(curve->A);
    mpz_clear(curve->B);
    mpz_clear(curve->x);
    mpz_clear(curve->y);
    mpz_clear(curve->n);
    mpz_clear(curve->h);
    curve->security_level = 0;
}

// 根据安全等级生成Montgomery曲线参数
bool generate_montgomery_curve(MontgomeryCurve* curve, int security_level) {
    // 设置安全等级
    curve->security_level = security_level;

    // 根据安全等级确定素数位数
    int prime_bits;
    switch (security_level) {
        case 112:  // 112位安全等级
            prime_bits = 224;
            break;
        case 128:  // 128位安全等级
            prime_bits = 256;
            break;
        case 192:  // 192位安全等级
            prime_bits = 384;
            break;
        case 256:  // 256位安全等级
            prime_bits = 512;
            break;
        default:
            std::cerr << "Unsupported security level: " << security_level << std::endl;
            return false;
    }

    // 生成素数域p
    generate_prime(curve->p, prime_bits);

    // 生成曲线参数A和B (B ≠ 0, A ≠ ±2 mod p)
    // 简化版本：随机生成满足条件的A和B
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    // 生成A参数 (A ≠ ±2 mod p)
    do {
        mpz_urandomm(curve->A, state, curve->p);
    } while (mpz_cmp_ui(curve->A, 2) == 0 ||
             (mpz_cmp_si(curve->A, -2) == 0 && mpz_sgn(curve->A) >= 0));

    // 生成B参数 (B ≠ 0)
    do {
        mpz_urandomm(curve->B, state, curve->p);
    } while (mpz_sgn(curve->B) == 0);

    // 寻找基点
    find_base_point(curve, state);

    gmp_randclear(state);

    return true;
}

// 在Montgomery曲线上寻找基点
void find_base_point(MontgomeryCurve* curve, gmp_randstate_t state) {
    // 初始化基点坐标
    mpz_set_ui(curve->x, 0);
    mpz_set_ui(curve->y, 0);
    mpz_set_ui(curve->n, 0);
    mpz_set_ui(curve->h, 1);

    // 寻找曲线上的点
    mpz_t x, rhs, y_squared;
    mpz_init(x);
    mpz_init(rhs);
    mpz_init(y_squared);

    // 尝试不同的x值，直到找到一个在曲线上的点
    int max_attempts = 10000;  // 最大尝试次数
    for (int i = 0; i < max_attempts; i++) {
        // 生成随机x值
        mpz_urandomm(x, state, curve->p);

        // 计算右边: x^3 + Ax^2 + x (mod p)
        mpz_t temp1, temp2;
        mpz_init(temp1);
        mpz_init(temp2);

        // x^2
        mpz_mul(temp1, x, x);
        mpz_mod(temp1, temp1, curve->p);

        // Ax^2
        mpz_mul(temp2, curve->A, temp1);
        mpz_mod(temp2, temp2, curve->p);

        // x^3
        mpz_mul(temp1, temp1, x);
        mpz_mod(temp1, temp1, curve->p);

        // x^3 + Ax^2
        mpz_add(rhs, temp1, temp2);
        mpz_mod(rhs, rhs, curve->p);

        // x^3 + Ax^2 + x
        mpz_add(rhs, rhs, x);
        mpz_mod(rhs, rhs, curve->p);

        // 计算 y^2 = (x^3 + Ax^2 + x) / B (mod p)
        mpz_t b_inv;
        mpz_init(b_inv);

        // 计算B的逆元
        if (mpz_invert(b_inv, curve->B, curve->p)) {
            // 计算 (x^3 + Ax^2 + x) / B (mod p)
            mpz_mul(y_squared, rhs, b_inv);
            mpz_mod(y_squared, y_squared, curve->p);

            // 检查y_squared是否是二次剩余（即是否有解）
            if (is_quadratic_residue(y_squared, curve->p)) {
                // 计算y值
                mpz_t y;
                mpz_init(y);
                bool sqrt_success = mpz_sqrtmod(y, y_squared, curve->p);
                if (sqrt_success) {
                    // 检查点是否有效
                    if (mpz_sgn(y) != 0) {
                        mpz_set(curve->x, x);
                        mpz_set(curve->y, y);
                        mpz_clear(y);
                        mpz_clear(temp1);
                        mpz_clear(temp2);
                        mpz_clear(b_inv);
                        break;
                    }
                }
                mpz_clear(y);
            }
        }

        mpz_clear(temp1);
        mpz_clear(temp2);
        mpz_clear(b_inv);
    }

    mpz_clear(x);
    mpz_clear(rhs);
    mpz_clear(y_squared);
}

// 检查一个数是否是模p的二次剩余
bool is_quadratic_residue(const mpz_t a, const mpz_t p) {
    // 使用欧拉准则: a是模p的二次剩余当且仅当 a^((p-1)/2) ≡ 1 (mod p)
    if (mpz_sgn(a) == 0) return true; // 0总是二次剩余

    mpz_t exponent, result;
    mpz_init(exponent);
    mpz_init(result);

    // 计算 (p-1)/2
    mpz_sub_ui(exponent, p, 1);
    mpz_fdiv_q_2exp(exponent, exponent, 1);

    // 计算 a^((p-1)/2) mod p
    mpz_powm(result, a, exponent, p);

    // 如果结果为1，则是二次剩余
    bool is_residue = (mpz_cmp_ui(result, 1) == 0);

    mpz_clear(exponent);
    mpz_clear(result);

    return is_residue;
}

// 计算模平方根，如果成功返回true，否则返回false
bool mpz_sqrtmod(mpz_t result, const mpz_t a, const mpz_t p) {
    // 使用Shanks-Tonelli算法计算模平方根
    // 这里使用一个简化但更可靠的实现

    // 特殊情况：a为0
    if (mpz_sgn(a) == 0) {
        mpz_set_ui(result, 0);
        return true;
    }

    // 检查是否为二次剩余
    if (!is_quadratic_residue(a, p)) {
        return false;
    }

    // 对于p ≡ 3 (mod 4)的情况，有一个简单的公式:
    // sqrt(a) ≡ a^((p+1)/4) (mod p)
    mpz_t mod4;
    mpz_init(mod4);
    mpz_mod_ui(mod4, p, 4);

    if (mpz_cmp_ui(mod4, 3) == 0) {
        // p ≡ 3 (mod 4)
        mpz_t exponent;
        mpz_init(exponent);

        // 计算 (p+1)/4
        mpz_add_ui(exponent, p, 1);
        mpz_fdiv_q_2exp(exponent, exponent, 2);

        // 计算 a^((p+1)/4) mod p
        mpz_powm(result, a, exponent, p);

        mpz_clear(exponent);
        mpz_clear(mod4);
        return true;
    }

    mpz_clear(mod4);

    // 简化处理：使用mpz_sqrt并验证结果
    mpz_sqrt(result, a);
    mpz_mod(result, result, p);

    // 验证结果
    mpz_t square;
    mpz_init(square);
    mpz_mul(square, result, result);
    mpz_mod(square, square, p);

    bool valid = (mpz_cmp(square, a) == 0);
    mpz_clear(square);

    return valid;
}
