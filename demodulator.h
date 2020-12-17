#ifndef DEMODULATOR_H
#define DEMODULATOR_H

#include "mainwindow.h"
#include <QApplication>

#include <fstream>
#include <vector>
#include <fftw3.h>
#include <iostream>
#include <string.h>
#include <complex>
#include <thread>
#include <unistd.h>
#include <iterator>
#include <random>

//#include "fdacoefs.h"
#include "fdacoefs_125K_to_100.h"
#include "fdacoefs_1M_to_100.h"
#include "fdacoefs_1M_to_125K.h"
#include "singenerator.h"
#include "preamblepoint.h"
#include "decoder.h"

#define X2_CHANNELS         1
#define CHANNEL_SAMPLE_RATE 100
#define SAMPLE_RATE         125000
#define SYM_RATE            100

#define SYM_LEN             (CHANNEL_SAMPLE_RATE/SYM_RATE)
#define X2                  (X2_CHANNELS * 2)
#define CHANNELS            (SAMPLE_RATE / CHANNEL_SAMPLE_RATE * X2)

#define DATA_LEN            256
#define MAX_ERRORS          3
#define PREAMBLE_LEN        32

class Demodulator
{
public:
    Demodulator(int inSempleRate);
    ~Demodulator();

    void addIQ(void* data, int sizeInSamples);
private:
    int numOfChannels;
    int decimationK;

    std::vector<float> wideSpectorData;
    std::vector<std::complex<float>>* channeslData;
    std::vector<float> decimatedData;

    fftwf_complex* fftw_in;
    fftwf_complex* fftw_out;
    fftwf_plan fftw_p;

#ifdef USE_WINDOW
    MainWindow* window = nullptr;
    std::thread* guiTh;
#endif
    std::thread* speedTh;

    float* i1;
    float* q1;

    float i2[CHANNELS], q2[CHANNELS];
    float corr1[CHANNELS], corr2[CHANNELS];

    uint32_t corr[CHANNELS];
    uint32_t inv_corr[CHANNELS];

    std::vector<PreamblePoint*> PreamblePointWithoutFullData;
    std::vector<PreamblePoint*> PreamblePointFullData;
    SinCosGenerator scg;

    int wideSpectorDataSamples = 0;
    int decimatedDataSamples = 0;
    int channeslDataSamples = 0;
    bool invers = true;

    Decoder* dec;

    void channelize();
    void decimation();
    int findPreambles(int ch);
    void dumpToFile(std::string fileName, void* data, size_t size);
    int bitDif(uint32_t, uint32_t);
    uint64_t toFromCoding(uint64_t);

    void freqCounter();
    void guiThread(int argc, char *argv[]);
    float noise;

    std::string uint32ToSring(uint32_t);
};

#endif // DEMODULATOR_H
