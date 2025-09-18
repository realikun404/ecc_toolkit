#include <gmp.h>
#include <iostream>

int main() {
    // 初始化大整数
    mpz_t a, b, result;
    mpz_init_set_str(a, "123456789012345678901234567890", 10);
    mpz_init_set_str(b, "987654321098765432109876543210", 10);
    mpz_init(result);

    // 执行大整数加法
    mpz_add(result, a, b);

    // 输出结果
    std::cout << "Testing GMP library installation..." << std::endl;
    std::cout << "a = ";
    mpz_out_str(stdout, 10, a);
    std::cout << std::endl;

    std::cout << "b = ";
    mpz_out_str(stdout, 10, b);
    std::cout << std::endl;

    std::cout << "a + b = ";
    mpz_out_str(stdout, 10, result);
    std::cout << std::endl;

    // 清理内存
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(result);

    std::cout << "GMP library is working correctly!" << std::endl;
    return 0;
}
