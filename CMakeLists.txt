# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(ecc_toolkit)

set(CMAKE_CXX_STANDARD 11)

# 查找GMP库
find_package(PkgConfig REQUIRED)
pkg_check_modules(GMP REQUIRED gmp)

# 包含头文件目录
include_directories(include)
include_directories(include/math)

# 添加源文件
set(SOURCES
        src/core/curve_generator.cpp
        src/core/montgomery_curve.cpp
        src/core/parameter_validator.cpp
        src/math/prime_utils.cpp
)

# 创建静态库
add_library(ecc_toolkit STATIC ${SOURCES})

# 链接GMP库
target_link_libraries(ecc_toolkit ${GMP_LIBRARIES})
target_include_directories(ecc_toolkit PRIVATE ${GMP_INCLUDE_DIRS})

# 添加测试可执行文件
add_executable(test_gmp tests/test_gmp.cpp)
target_link_libraries(test_gmp ecc_toolkit ${GMP_LIBRARIES})

add_executable(test_montgomery tests/test_montgomery.cpp)
target_link_libraries(test_montgomery ecc_toolkit ${GMP_LIBRARIES})

add_executable(test_generator tests/test_generator.cpp)
target_link_libraries(test_generator ecc_toolkit ${GMP_LIBRARIES})

add_executable(test_validator tests/test_validator.cpp)
target_link_libraries(test_validator ecc_toolkit ${GMP_LIBRARIES})

# 新增的阶验证测试
add_executable(test_order_validation tests/test_order_validation.cpp)
target_link_libraries(test_order_validation ecc_toolkit ${GMP_LIBRARIES})

# 启用测试
enable_testing()

# 添加测试
add_test(NAME test_gmp COMMAND test_gmp)
add_test(NAME test_montgomery COMMAND test_montgomery)
add_test(NAME test_generator COMMAND test_generator)
add_test(NAME test_validator COMMAND test_validator)
add_test(NAME test_order_validation COMMAND test_order_validation)
