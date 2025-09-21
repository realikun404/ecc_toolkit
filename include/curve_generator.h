#ifndef CURVE_GENERATOR_H
#define CURVE_GENERATOR_H

#include "types.h"
#include <gmp.h>

// 根据安全等级生成Montgomery曲线参数
bool generate_montgomery_curve(MontgomeryCurve* curve, int security_level);

// 寻找基点
void find_base_point(MontgomeryCurve* curve, gmp_randstate_t state);

// 检查是否为二次剩余
bool is_quadratic_residue(const mpz_t a, const mpz_t p);

// 计算模平方根（简化版）
void mpz_sqrtmod(mpz_t result, const mpz_t a, const mpz_t p);

#endif // CURVE_GENERATOR_H
