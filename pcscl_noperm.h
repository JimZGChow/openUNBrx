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
#include <stdio.h>
#include <iostream>

struct pcscl_list
{
    std::vector<bool> x;
    int idx;
};

float phi(float x);
inline float sign(float x);
inline float jacoblog(float x);
inline float cnop(float a, float b);
inline float cnop(bool a, float b);
inline float vnop(float a, float b);
std::vector<std::vector<uint8_t>> polar_transform_noperm(std::vector<std::vector<uint8_t>> data);
std::vector<pcscl_list> pcscl(const std::vector<std::vector<float> >& y, std::vector<int>::const_iterator f_it);
std::vector<std::vector<uint8_t>> extract_with_filter(const std::vector<std::vector<uint8_t>>& data, const uint16_t* filtr, int K);
std::vector<std::vector<uint8_t>> pcscl_prep(int m, int L, const std::vector<uint8_t>& data, std::vector<std::vector<float>>* prob, uint8_t* info_bit_pattern);
std::vector<std::vector<uint8_t>> pcscl_prep(int m, int L, const std::vector<float>& data, std::vector<std::vector<float>>* prob, uint8_t* info_bit_pattern);
std::vector<uint8_t> crc_ok_array(uint32_t polynom, std::vector<std::vector<uint8_t>> a);
uint8_t crc_ok(uint32_t polynom, std::vector<uint8_t> a);

#endif // PCSCL_NOPERM_H
