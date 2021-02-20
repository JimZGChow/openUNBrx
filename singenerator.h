#ifndef SINGENERATOR_H
#define SINGENERATOR_H

#include <math.h>
#include <stdio.h>
#include <iostream>

class SinCosGenerator
{
public:
    SinCosGenerator(int sampleRate, float freq);
    SinCosGenerator();

    float nextSin();
    float nextCos();
    void nextStep();
    void setSampleRate(int sampleRate);
    void setFreq(float freq);
    void reset() {delta = 0;};

private:
    double sr;
    double fr;
    double div;
    double delta;
};

#endif // SINGENERATOR_H
