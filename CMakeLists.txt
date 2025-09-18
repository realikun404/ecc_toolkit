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

# 创建可执行文件用于测试GMP
add_executable(test_gmp test_gmp.cpp)

# 链接GMP库
target_link_libraries(test_gmp ${GMP_LIBRARIES})

# 其他原有的项目配置...
