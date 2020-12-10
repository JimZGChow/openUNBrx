#ifndef PREAMBLEPOINT_H
#define PREAMBLEPOINT_H

#include <vector>
#include <complex>

class PreamblePoint
{
public:
    PreamblePoint();

    int preableErrors;
    int channel;
    std::vector<std::complex<float>> data;
    uint8_t inv;
    int added = 0;
    float noise;
    int pos;

    enum {
        invBits = 1,
        invIQ = 2
    };
};

#endif // PREAMBLEPOINT_H
