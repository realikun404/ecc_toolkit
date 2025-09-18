#ifndef MONTGOMERY_CURVE_H
#define MONTGOMERY_CURVE_H

#include "types.h"

// Montgomery曲线上的点
typedef struct {
    mpz_t x;
    mpz_t y;
    int is_infinity; // 1表示无穷远点，0表示普通点
} MontgomeryPoint;

// 初始化Montgomery点
void init_montgomery_point(MontgomeryPoint* point);

// 清理Montgomery点
void clear_montgomery_point(MontgomeryPoint* point);

// 设置点为无穷远点
void set_point_infinity(MontgomeryPoint* point);

// 检查点是否为无穷远点
int is_point_infinity(const MontgomeryPoint* point);

// Montgomery曲线点加法: R = P + Q
bool montgomery_point_add(MontgomeryPoint* R, const MontgomeryPoint* P, const MontgomeryPoint* Q,
                         const mpz_t A, const mpz_t p);

// Montgomery曲线点倍乘: R = k * P
bool montgomery_point_multiply(MontgomeryPoint* R, const mpz_t k, const MontgomeryPoint* P,
                              const mpz_t A, const mpz_t p);

// 检查点是否在Montgomery曲线上
bool is_point_on_curve(const MontgomeryPoint* point, const mpz_t A, const mpz_t B, const mpz_t p);

#endif // MONTGOMERY_CURVE_H
