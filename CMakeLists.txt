cmake_minimum_required(VERSION 3.10)
project(ecc_toolkit)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找GMP库
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED gmp)

# 包含头文件目录
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${GMP_INCLUDE_DIRS})

# 添加编译选项
add_compile_options(${GMP_CFLAGS_OTHER})

# 创建GMP测试可执行文件
add_executable(test_gmp tests/test_gmp.cpp)

# 创建曲线生成器测试可执行文件
add_executable(test_generator tests/test_generator.cpp
    src/core/curve_generator.cpp
    src/core/parameter_validator.cpp
    src/core/montgomery_curve.cpp
    src/math/prime_utils.cpp)

# 创建参数验证器测试可执行文件
add_executable(test_validator tests/test_validator.cpp
    src/core/parameter_validator.cpp
    src/core/curve_generator.cpp
    src/math/prime_utils.cpp)

# 创建Montgomery曲线运算测试可执行文件
add_executable(test_montgomery tests/test_montgomery.cpp
    src/core/montgomery_curve.cpp
    src/core/curve_generator.cpp
    src/core/parameter_validator.cpp
    src/math/prime_utils.cpp)

# 链接GMP库
target_link_libraries(test_gmp ${GMP_LIBRARIES})
target_link_libraries(test_generator ${GMP_LIBRARIES})
target_link_libraries(test_validator ${GMP_LIBRARIES})
target_link_libraries(test_montgomery ${GMP_LIBRARIES})
