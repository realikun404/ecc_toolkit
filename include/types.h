#ifndef TYPES_H
#define TYPES_H

#include <gmp.h>

// 定义Montgomery曲线参数结构
typedef struct {
    mpz_t p;  // 素数域
    mpz_t A;  // 曲线参数A
    mpz_t B;  // 曲线参数B
    mpz_t x;  // 基点x坐标
    mpz_t y;  // 基点y坐标
    mpz_t n;  // 基点阶数
    mpz_t h;  // 余因子
    int security_level; // 安全等级（以位为单位）
} MontgomeryCurve;

// 初始化Montgomery曲线参数
void init_montgomery_curve(MontgomeryCurve* curve);

// 清理Montgomery曲线参数
void clear_montgomery_curve(MontgomeryCurve* curve);

#endif // TYPES_H
