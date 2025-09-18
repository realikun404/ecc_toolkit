#ifndef PRIME_UTILS_H
#define PRIME_UTILS_H

#include <gmp.h>

// 生成指定位数的素数
void generate_prime(mpz_t prime, int bits);

#endif // PRIME_UTILS_H
