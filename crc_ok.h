#ifndef CRC_OK_H
#define CRC_OK_H

#include <vector>
#include <stdint.h>

std::vector<bool> crc_ok(uint32_t polynom, std::vector<std::vector<bool>> a);

#endif // CRC_OK_H
