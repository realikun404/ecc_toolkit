#ifndef CURVE_GENERATOR_H
#define CURVE_GENERATOR_H

#include "types.h"

// 根据安全等级生成Montgomery曲线参数
bool generate_montgomery_curve(MontgomeryCurve* curve, int security_level);

#endif // CURVE_GENERATOR_H
