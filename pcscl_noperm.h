#ifndef PCSCL_NOPERM_H
#define PCSCL_NOPERM_H

//#define MEX

#ifdef MEX
#include "mex.h"
#endif
#include <vector>
#include <string>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cassert>
#include <cmath>

struct pcscl_list
{
    std::vector<bool> x;
    int idx;
};

void test();

float phi(float x);
inline float sign(float x);
inline float jacoblog(float x);
inline float cnop(float a, float b);
inline float cnop(bool a, float b);
inline float vnop(float a, float b);
std::vector<pcscl_list> pcscl(const std::vector<std::vector<float> > &y, std::vector<int>::const_iterator f_it);
std::vector<bool> decode64(const std::vector<bool>& data);

#endif // PCSCL_NOPERM_H
