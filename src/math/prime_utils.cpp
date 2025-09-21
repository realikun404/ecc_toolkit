#include "math/prime_utils.h"
#include <gmp.h>
#include <ctime>

// 生成指定位数的素数
void generate_prime(mpz_t prime, int bits) {
    // 初始化随机数状态
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    // 生成指定位数的随机数
    mpz_urandomb(prime, state, bits);

    // 设置最高位和最低位为1确保是奇数且具有指定位数
    mpz_setbit(prime, bits - 1); // 设置最高位
    mpz_setbit(prime, 0);        // 设置最低位

    // 寻找下一个素数
    while (!mpz_probab_prime_p(prime, 25)) {  // 使用Miller-Rabin测试检查是否为素数
        mpz_add_ui(prime, prime, 2);          // 如果不是素数，则加2继续检查下一个奇数
    }

    // 清理随机数状态
    gmp_randclear(state);
}
