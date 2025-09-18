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

    // 暂时不设置基点，因为计算在曲线上的点比较复杂
    // 在实际应用中，需要更复杂的算法来找到曲线上的点和相应的阶数
    mpz_set_ui(curve->x, 0);  // 0表示未设置
    mpz_set_ui(curve->y, 0);  // 0表示未设置
    mpz_set_ui(curve->n, 0);  // 0表示未设置
    mpz_set_ui(curve->h, 1);

    gmp_randclear(state);

    return true;
}
