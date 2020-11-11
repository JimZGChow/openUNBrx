#include "singenerator.h"

SinCosGenerator::SinCosGenerator(int sampleRate, float freq) {
    sr = sampleRate;
    fr = freq;

    div = freq/sampleRate;
    delta = 0;
}

float SinCosGenerator::nextSin() {
    return sin(delta);
}

float SinCosGenerator::nextCos() {
    return cos(delta);
}

void SinCosGenerator::nextStep() {
    delta += div * 2 * M_PI;
}

void SinCosGenerator::setFreq(float freq) {
    delta = 0;
    fr = freq;
    div = freq/sr;
}

void SinCosGenerator::setSampleRate(int sampleRate) {
    delta = 0;
    sr = sampleRate;
    div = fr/sampleRate;
}
