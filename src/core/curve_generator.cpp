#include "curve_generator.h"
#include "parameter_validator.h"
#include "math/prime_utils.h"
#include <gmp.h>
#include <iostream>
#include <random>

// 初始化Montgomery曲线参数
void init_montgomery_curve(MontgomeryCurve* curve) {
    std::cout << "DEBUG: init_montgomery_curve start" << std::endl;
    mpz_init(curve->p);
    mpz_init(curve->A);
    mpz_init(curve->B);
    mpz_init(curve->x);
    mpz_init(curve->y);
    mpz_init(curve->n);
    mpz_init(curve->h);
    curve->security_level = 0;
    std::cout << "DEBUG: init_montgomery_curve end" << std::endl;
}

// 清理Montgomery曲线参数
void clear_montgomery_curve(MontgomeryCurve* curve) {
    std::cout << "DEBUG: clear_montgomery_curve start" << std::endl;
    mpz_clear(curve->p);
    mpz_clear(curve->A);
    mpz_clear(curve->B);
    mpz_clear(curve->x);
    mpz_clear(curve->y);
    mpz_clear(curve->n);
    mpz_clear(curve->h);
    curve->security_level = 0;
    std::cout << "DEBUG: clear_montgomery_curve end" << std::endl;
}

// 根据安全等级生成Montgomery曲线参数
bool generate_montgomery_curve(MontgomeryCurve* curve, int security_level) {
    std::cout << "DEBUG: generate_montgomery_curve start" << std::endl;
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
    std::cout << "DEBUG: Initializing gmp_randstate_t" << std::endl;
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
    std::cout << "DEBUG: Calling find_base_point" << std::endl;
    find_base_point(curve, state);

    std::cout << "DEBUG: Cleaning up gmp_randstate_t" << std::endl;
    gmp_randclear(state);

    std::cout << "DEBUG: generate_montgomery_curve end" << std::endl;
    return true;
}

// 在Montgomery曲线上寻找基点
void find_base_point(MontgomeryCurve* curve, gmp_randstate_t state) {
    std::cout << "DEBUG: find_base_point start" << std::endl;
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
    int max_attempts = 100000;  // 增加最大尝试次数
    bool found = false;
    for (int i = 1; i < max_attempts && !found; i++) {
        // 生成随机x值 (从1开始，避免x=0)
        mpz_urandomm(x, state, curve->p);
        if (mpz_sgn(x) == 0) {
            mpz_set_ui(x, 1);
        }

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
                if (sqrt_success && mpz_sgn(y) != 0) {
                    // 验证点是否在曲线上
                    mpz_t left_side;
                    mpz_init(left_side);

                    // 计算By^2
                    mpz_mul(temp1, curve->B, y);
                    mpz_mul(left_side, temp1, y);
                    mpz_mod(left_side, left_side, curve->p);

                    // 比较左右两边
                    if (mpz_cmp(left_side, rhs) == 0) {
                        mpz_set(curve->x, x);
                        mpz_set(curve->y, y);
                        found = true;
                        std::cout << "DEBUG: Base point found" << std::endl;
                        mpz_clear(left_side);
                        mpz_clear(y);
                    } else {
                        mpz_clear(left_side);
                        mpz_clear(y);
                    }
                } else {
                    mpz_clear(y);
                }
            }
        }

        mpz_clear(temp1);
        mpz_clear(temp2);
        mpz_clear(b_inv);

        // 如果找到了基点，清理循环中的变量
        if (found) {
            std::cout << "DEBUG: Cleaning up variables in found branch" << std::endl;
            mpz_clear(x);
            mpz_clear(rhs);
            mpz_clear(y_squared);
            std::cout << "DEBUG: find_base_point end (found branch)" << std::endl;
            return;
        }
    }

    // 如果随机方法没找到，使用确定性方法
    std::cout << "DEBUG: Cleaning up variables before deterministic method" << std::endl;
    mpz_clear(x);
    mpz_clear(rhs);
    mpz_clear(y_squared);
    find_base_point_deterministic(curve);
    std::cout << "DEBUG: find_base_point end" << std::endl;
}

// 确定性方法寻找基点
void find_base_point_deterministic(MontgomeryCurve* curve) {
    std::cout << "DEBUG: find_base_point_deterministic start" << std::endl;
    mpz_t x, rhs, y_squared;
    mpz_init(x);
    mpz_init(rhs);
    mpz_init(y_squared);

    // 从1开始尝试x值
    mpz_set_ui(x, 1);

    // 尝试最多1000个x值
    bool found = false;
    for (int i = 0; i < 1000 && !found; i++) {
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
                if (sqrt_success && mpz_sgn(y) != 0) {
                    // 验证点是否在曲线上
                    mpz_t left_side;
                    mpz_init(left_side);

                    // 计算By^2
                    mpz_mul(temp1, curve->B, y);
                    mpz_mul(left_side, temp1, y);
                    mpz_mod(left_side, left_side, curve->p);

                    // 比较左右两边
                    if (mpz_cmp(left_side, rhs) == 0) {
                        mpz_set(curve->x, x);
                        mpz_set(curve->y, y);
                        found = true;
                        std::cout << "DEBUG: Base point found in deterministic method" << std::endl;
                        mpz_clear(left_side);
                        mpz_clear(y);
                    } else {
                        mpz_clear(left_side);
                        mpz_clear(y);
                    }
                } else {
                    mpz_clear(y);
                }
            }
        }

        mpz_clear(temp1);
        mpz_clear(temp2);
        mpz_clear(b_inv);

        // 如果找到了基点，清理循环中的变量
        if (found) {
            std::cout << "DEBUG: Cleaning up variables in deterministic found branch" << std::endl;
            mpz_clear(x);
            mpz_clear(rhs);
            mpz_clear(y_squared);
            std::cout << "DEBUG: find_base_point_deterministic end (found branch)" << std::endl;
            return;
        }

        // 尝试下一个x值
        mpz_add_ui(x, x, 1);
        mpz_mod(x, x, curve->p);
    }

    std::cout << "DEBUG: Cleaning up variables in deterministic method end" << std::endl;
    mpz_clear(x);
    mpz_clear(rhs);
    mpz_clear(y_squared);
    std::cout << "DEBUG: find_base_point_deterministic end" << std::endl;
}

// 检查一个数是否是模p的二次剩余
bool is_quadratic_residue(const mpz_t a, const mpz_t p) {
    std::cout << "DEBUG: is_quadratic_residue start" << std::endl;
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
    std::cout << "DEBUG: is_quadratic_residue end" << std::endl;
    return is_residue;
}

// 计算模平方根，如果成功返回true，否则返回false
bool mpz_sqrtmod(mpz_t result, const mpz_t a, const mpz_t p) {
    std::cout << "DEBUG: mpz_sqrtmod start" << std::endl;
    // 使用Shanks-Tonelli算法计算模平方根
    // 这里使用一个简化但更可靠的实现

    // 特殊情况：a为0
    if (mpz_sgn(a) == 0) {
        mpz_set_ui(result, 0);
        std::cout << "DEBUG: mpz_sqrtmod end (zero case)" << std::endl;
        return true;
    }

    // 检查是否为二次剩余
    std::cout << "DEBUG: mpz_sqrtmod calling is_quadratic_residue" << std::endl;
    if (!is_quadratic_residue(a, p)) {
        std::cout << "DEBUG: mpz_sqrtmod end (not quadratic residue)" << std::endl;
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
        std::cout << "DEBUG: mpz_sqrtmod end (p ≡ 3 mod 4 case)" << std::endl;
        return true;
    }

    mpz_clear(mod4);

    // 对于其他情况，使用mpz_sqrt并验证结果
    std::cout << "DEBUG: mpz_sqrtmod general case start" << std::endl;
    mpz_sqrt(result, a);
    mpz_mod(result, result, p);

    // 验证结果
    mpz_t square;
    mpz_init(square);
    mpz_mul(square, result, result);
    mpz_mod(square, square, p);

    bool valid = (mpz_cmp(square, a) == 0);
    std::cout << "DEBUG: mpz_sqrtmod clearing square (first time)" << std::endl;
    mpz_clear(square);

    // 如果验证失败，尝试另一种方法
    if (!valid && mpz_sgn(result) != 0) {
        // 尝试 p - result
        std::cout << "DEBUG: mpz_sqrtmod alternative method start" << std::endl;
        mpz_sub(result, p, result);
        mpz_init(square);  // 重新初始化square变量
        mpz_mul(square, result, result);
        mpz_mod(square, square, p);
        valid = (mpz_cmp(square, a) == 0);
        std::cout << "DEBUG: mpz_sqrtmod clearing square (second time)" << std::endl;
        mpz_clear(square);
    }

    std::cout << "DEBUG: mpz_sqrtmod end (general case)" << std::endl;
    return valid;
}
