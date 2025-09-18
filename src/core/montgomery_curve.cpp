#include "montgomery_curve.h"
#include <gmp.h>
#include <iostream>

// 初始化Montgomery点
void init_montgomery_point(MontgomeryPoint* point) {
    mpz_init(point->x);
    mpz_init(point->y);
    point->is_infinity = 0;
}

// 清理Montgomery点
void clear_montgomery_point(MontgomeryPoint* point) {
    mpz_clear(point->x);
    mpz_clear(point->y);
    point->is_infinity = 0;
}

// 设置点为无穷远点
void set_point_infinity(MontgomeryPoint* point) {
    point->is_infinity = 1;
}

// 检查点是否为无穷远点
int is_point_infinity(const MontgomeryPoint* point) {
    return point->is_infinity;
}

// Montgomery曲线点加法: R = P + Q
// 使用Montgomery梯度公式
bool montgomery_point_add(MontgomeryPoint* R, const MontgomeryPoint* P, const MontgomeryPoint* Q,
                         const mpz_t A, const mpz_t p) {
    // 如果P是无穷远点，则R = Q
    if (is_point_infinity(P)) {
        mpz_set(R->x, Q->x);
        mpz_set(R->y, Q->y);
        R->is_infinity = Q->is_infinity;
        return true;
    }

    // 如果Q是无穷远点，则R = P
    if (is_point_infinity(Q)) {
        mpz_set(R->x, P->x);
        mpz_set(R->y, P->y);
        R->is_infinity = P->is_infinity;
        return true;
    }

    // 如果P和Q是相同的点，则使用点倍乘公式
    if (mpz_cmp(P->x, Q->x) == 0 && mpz_cmp(P->y, Q->y) == 0) {
        // 点倍乘情况
        // λ = (3*x^2 + 2*A*x + 1) / (2*B*y)
        mpz_t lambda, temp1, temp2, temp3, denominator;
        mpz_init(lambda);
        mpz_init(temp1);
        mpz_init(temp2);
        mpz_init(temp3);
        mpz_init(denominator);

        // 计算分子: 3*x^2 + 2*A*x + 1
        mpz_mul(temp1, P->x, P->x);              // x^2
        mpz_mul_ui(temp1, temp1, 3);             // 3*x^2
        mpz_mod(temp1, temp1, p);

        mpz_mul(temp2, A, P->x);                 // A*x
        mpz_mul_ui(temp2, temp2, 2);             // 2*A*x
        mpz_mod(temp2, temp2, p);

        mpz_add(temp3, temp1, temp2);            // 3*x^2 + 2*A*x
        mpz_mod(temp3, temp3, p);

        mpz_add_ui(temp3, temp3, 1);             // 3*x^2 + 2*A*x + 1
        mpz_mod(temp3, temp3, p);

        // 计算分母: 2*y
        mpz_mul_ui(denominator, P->y, 2);        // 2*y
        mpz_mod(denominator, denominator, p);

        // 计算λ = 分子 / 分母
        mpz_invert(temp1, denominator, p);       // 求分母的逆元
        mpz_mul(lambda, temp3, temp1);           // λ = 分子 * 逆元
        mpz_mod(lambda, lambda, p);

        // 计算新的x坐标: x3 = λ^2 - A - x1 - x2 (这里x1 = x2)
        mpz_mul(temp1, lambda, lambda);          // λ^2
        mpz_mod(temp1, temp1, p);

        mpz_sub(temp2, temp1, A);                // λ^2 - A
        mpz_mod(temp2, temp2, p);

        mpz_sub(temp3, temp2, P->x);             // λ^2 - A - x1
        mpz_mod(temp3, temp3, p);

        mpz_sub(R->x, temp3, P->x);              // λ^2 - A - x1 - x2
        mpz_mod(R->x, R->x, p);

        // 计算新的y坐标: y3 = λ*(x1 - x3) - y1
        mpz_sub(temp1, P->x, R->x);              // x1 - x3
        mpz_mod(temp1, temp1, p);

        mpz_mul(temp2, lambda, temp1);           // λ*(x1 - x3)
        mpz_mod(temp2, temp2, p);

        mpz_sub(R->y, temp2, P->y);              // λ*(x1 - x3) - y1
        mpz_mod(R->y, R->y, p);

        R->is_infinity = 0;

        mpz_clear(lambda);
        mpz_clear(temp1);
        mpz_clear(temp2);
        mpz_clear(temp3);
        mpz_clear(denominator);
        return true;
    }

    // 如果P和Q的x坐标相同但y坐标不同，则结果为无穷远点
    if (mpz_cmp(P->x, Q->x) == 0 && mpz_cmp(P->y, Q->y) != 0) {
        set_point_infinity(R);
        return true;
    }

    // 一般情况的点加法
    // λ = (y2 - y1) / (x2 - x1)
    mpz_t lambda, temp1, temp2, temp3;
    mpz_init(lambda);
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(temp3);

    // 计算λ = (y2 - y1) / (x2 - x1)
    mpz_sub(temp1, Q->y, P->y);                 // y2 - y1
    mpz_mod(temp1, temp1, p);

    mpz_sub(temp2, Q->x, P->x);                 // x2 - x1
    mpz_mod(temp2, temp2, p);

    mpz_invert(temp3, temp2, p);                // 求(x2 - x1)的逆元
    mpz_mul(lambda, temp1, temp3);              // λ = (y2 - y1) * 逆元
    mpz_mod(lambda, lambda, p);

    // 计算新的x坐标: x3 = λ^2 - A - x1 - x2
    mpz_mul(temp1, lambda, lambda);             // λ^2
    mpz_mod(temp1, temp1, p);

    mpz_sub(temp2, temp1, A);                   // λ^2 - A
    mpz_mod(temp2, temp2, p);

    mpz_sub(temp3, temp2, P->x);                // λ^2 - A - x1
    mpz_mod(temp3, temp3, p);

    mpz_sub(R->x, temp3, Q->x);                 // λ^2 - A - x1 - x2
    mpz_mod(R->x, R->x, p);

    // 计算新的y坐标: y3 = λ*(x1 - x3) - y1
    mpz_sub(temp1, P->x, R->x);                 // x1 - x3
    mpz_mod(temp1, temp1, p);

    mpz_mul(temp2, lambda, temp1);              // λ*(x1 - x3)
    mpz_mod(temp2, temp2, p);

    mpz_sub(R->y, temp2, P->y);                 // λ*(x1 - x3) - y1
    mpz_mod(R->y, R->y, p);

    R->is_infinity = 0;

    mpz_clear(lambda);
    mpz_clear(temp1);
    mpz_clear(temp2);
    mpz_clear(temp3);
    return true;
}

// Montgomery曲线点倍乘: R = k * P (使用double-and-add算法)
bool montgomery_point_multiply(MontgomeryPoint* R, const mpz_t k, const MontgomeryPoint* P,
                              const mpz_t A, const mpz_t p) {
    // 如果k为0或P是无穷远点，则结果为无穷远点
    if (mpz_sgn(k) == 0 || is_point_infinity(P)) {
        set_point_infinity(R);
        return true;
    }

    // 初始化结果为无穷远点
    set_point_infinity(R);

    // 创建临时点用于计算
    MontgomeryPoint Q;
    init_montgomery_point(&Q);
    mpz_set(Q.x, P->x);
    mpz_set(Q.y, P->y);
    Q.is_infinity = P->is_infinity;

    // 使用double-and-add算法
    mpz_t temp_k;
    mpz_init_set(temp_k, k);

    while (mpz_sgn(temp_k) > 0) {
        // 如果当前位为1，则加到结果中
        if (mpz_odd_p(temp_k)) {
            if (is_point_infinity(R)) {
                mpz_set(R->x, Q.x);
                mpz_set(R->y, Q.y);
                R->is_infinity = Q.is_infinity;
            } else {
                MontgomeryPoint temp;
                init_montgomery_point(&temp);
                montgomery_point_add(&temp, R, &Q, A, p);
                mpz_set(R->x, temp.x);
                mpz_set(R->y, temp.y);
                R->is_infinity = temp.is_infinity;
                clear_montgomery_point(&temp);
            }
        }

        // Q = 2*Q
        MontgomeryPoint temp;
        init_montgomery_point(&temp);
        montgomery_point_add(&temp, &Q, &Q, A, p);
        mpz_set(Q.x, temp.x);
        mpz_set(Q.y, temp.y);
        Q.is_infinity = temp.is_infinity;
        clear_montgomery_point(&temp);

        // k = k >> 1
        mpz_fdiv_q_2exp(temp_k, temp_k, 1);
    }

    clear_montgomery_point(&Q);
    mpz_clear(temp_k);
    return true;
}

// 检查点是否在Montgomery曲线上: By^2 = x^3 + Ax^2 + x
bool is_point_on_curve(const MontgomeryPoint* point, const mpz_t A, const mpz_t B, const mpz_t p) {
    // 无穷远点被认为在曲线上
    if (is_point_infinity(point)) {
        return true;
    }

    mpz_t left, right, temp1, temp2;
    mpz_init(left);
    mpz_init(right);
    mpz_init(temp1);
    mpz_init(temp2);

    // 计算左边: By^2
    mpz_mul(temp1, point->y, point->y);      // y^2
    mpz_mod(temp1, temp1, p);
    mpz_mul(left, B, temp1);                 // B*y^2
    mpz_mod(left, left, p);

    // 计算右边: x^3 + Ax^2 + x
    mpz_mul(temp1, point->x, point->x);      // x^2
    mpz_mod(temp1, temp1, p);

    mpz_mul(temp2, A, temp1);                // A*x^2
    mpz_mod(temp2, temp2, p);

    mpz_mul(temp1, temp1, point->x);         // x^3
    mpz_mod(temp1, temp1, p);

    mpz_add(right, temp1, temp2);            // x^3 + A*x^2
    mpz_mod(right, right, p);

    mpz_add(temp1, right, point->x);         // x^3 + A*x^2 + x
    mpz_mod(temp1, temp1, p);

    // 比较左右两边
    bool result = (mpz_cmp(left, temp1) == 0);

    mpz_clear(left);
    mpz_clear(right);
    mpz_clear(temp1);
    mpz_clear(temp2);

    return result;
}
