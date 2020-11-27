#ifndef DS_H
#define DS_H

#include <vector>
#include <stdint.h>
#include <string>

std::vector<uint8_t> getVectorFromArray(uint8_t* array, int size);
std::vector<uint8_t> getVectorFromStringHex(std::string s);
std::vector<uint8_t> getVectorFromStringBin(std::string s);
std::vector<uint8_t> getRandomVector(int size);
std::string getStringHexFromVector(std::vector<uint8_t>);
std::string getStringBinFromVector(std::vector<uint8_t>);
int getErrors(std::vector<uint8_t> a, std::vector<uint8_t> b);


#endif // DS_H
