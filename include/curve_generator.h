#ifndef CURVE_GENERATOR_H
#define CURVE_GENERATOR_H

#include "types.h"
#include <gmp.h>

// 根据安全等级生成Montgomery曲线参数
bool generate_montgomery_curve(MontgomeryCurve* curve, int security_level);

// 辅助函数：验证基点是否在曲线上
bool is_point_on_curve_valid(const MontgomeryCurve* curve);


// 寻找基点
void find_base_point(MontgomeryCurve* curve, gmp_randstate_t state);

// 确定性方法寻找基点
void find_base_point_deterministic(MontgomeryCurve* curve);

// 检查是否为二次剩余
bool is_quadratic_residue(const mpz_t a, const mpz_t p);

// 计算模平方根，如果成功返回true，否则返回false
bool mpz_sqrtmod(mpz_t result, const mpz_t a, const mpz_t p);

// 计算曲线阶数（点数）
bool calculate_curve_order(const MontgomeryCurve* curve, mpz_t order);
// 计算基点的阶
bool calculate_base_point_order(const MontgomeryCurve* curve, mpz_t order);
// 验证基点是否为生成元
bool verify_generator(const MontgomeryCurve* curve);


#endif // CURVE_GENERATOR_H
