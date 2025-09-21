#!/bin/bash

# 测试脚本，用于验证double free问题

echo "清理并重新编译项目..."
cd ~/ecc_toolkit
rm -rf build
mkdir build
cd build
cmake .. > /dev/null 2>&1
make > /dev/null 2>&1

echo "测试 curve_generator 10次运行..."
for i in {1..10}
do
    echo "运行第 $i 次测试..."
    ./test_generator > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "第 $i 次测试失败!"
        exit 1
    fi
done

echo "测试 test_montgomery 10次运行..."
for i in {1..10}
do
    echo "运行第 $i 次测试..."
    ./test_montgomery > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "第 $i 次测试失败!"
        exit 1
    fi
done

echo "测试 test_validator 10次运行..."
for i in {1..10}
do
    echo "运行第 $i 次测试..."
    ./test_validator > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "第 $i 次测试失败!"
        exit 1
    fi
done

echo "所有测试通过，double free问题已解决!"
