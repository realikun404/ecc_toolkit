#ifndef PARAMETER_VALIDATOR_H
#define PARAMETER_VALIDATOR_H

#include "types.h"

// 验证Montgomery曲线参数的合法性
bool validate_montgomery_curve(const MontgomeryCurve* curve);

// 验证素数域p是否为有效素数
bool validate_prime_field(const mpz_t p);

// 验证曲线参数A和B是否满足Montgomery曲线条件
bool validate_curve_parameters(const mpz_t A, const mpz_t B, const mpz_t p);

// 验证基点是否在曲线上
bool validate_base_point(const mpz_t x, const mpz_t y, const mpz_t A, const mpz_t B, const mpz_t p);

#endif // PARAMETER_VALIDATOR_H
