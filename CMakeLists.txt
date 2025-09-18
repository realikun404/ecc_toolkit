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

# 收集源文件
file(GLOB_RECURSE SOURCES "src/**/*.cpp" "src/*.cpp")

# 创建GMP测试可执行文件
add_executable(test_gmp tests/test_gmp.cpp)

# 创建曲线生成器测试可执行文件
add_executable(test_generator tests/test_generator.cpp ${SOURCES})

# 链接GMP库
target_link_libraries(test_gmp ${GMP_LIBRARIES})
target_link_libraries(test_generator ${GMP_LIBRARIES})
