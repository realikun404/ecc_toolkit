#!/bin/bash

# 测试脚本，用于多次运行测试程序以检测间歇性问题

echo "开始多次测试以检测间歇性问题..."

# 编译项目
echo "编译项目..."
cd build
make clean
make

# 检查编译是否成功
if [ $? -ne 0 ]; then
    echo "编译失败!"
    exit 1
fi

# 运行测试的次数
TEST_COUNT=50

echo "将运行 $TEST_COUNT 次测试..."

# 测试test_gmp
echo "测试 test_gmp..."
for i in $(seq 1 $TEST_COUNT); do
    echo -n "运行 test_gmp 第 $i 次... "
    ./test_gmp > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "失败!"
        exit 1
