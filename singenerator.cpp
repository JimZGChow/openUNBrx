#include "singenerator.h"

SinCosGenerator::SinCosGenerator(int sampleRate, float freq) {
    sr = sampleRate;
    fr = freq;

    div = fr/sr;
    delta = 0;
}

SinCosGenerator::SinCosGenerator() {
    sr = 0;
    fr = 0;

    div = 1;
    delta = 0;
}

float SinCosGenerator::nextSin() {
    if (fr == 0.0)
        return 1;
    return sin(delta);
}

float SinCosGenerator::nextCos() {
    if (fr == 0.0)
        return 1;
    return cos(delta);
}

void SinCosGenerator::nextStep() {

    delta += div * 2 * M_PI;
}

void SinCosGenerator::setFreq(float freq) {
    //delta = 0;
    fr = freq;
    div = fr/sr;
}

void SinCosGenerator::setSampleRate(int sampleRate) {
    delta = 0;
    sr = sampleRate;
    div = fr/sr;
}
