#ifndef DEMODULATOR_H
#define DEMODULATOR_H

#include "mainwindow.h"
#include <QApplication>

#include <vector>
#include <fftw3.h>
#include <iostream>
#include <string.h>
#include <complex>
#include <thread>

//#include "fdacoefs.h"
#include "fdacoefs_125K_to_100.h"
#include "fdacoefs_1M_to_100.h"
#include "fdacoefs_1M_to_125K.h"

#define X2_CHANNELS         1
#define CHANNEL_SAMPLE_RATE 100
#define SAMPLE_RATE         125000
#define SYM_RATE            100

#define SYM_LEN             (CHANNEL_SAMPLE_RATE/SYM_RATE)
#define X2                  (X2_CHANNELS * 2)

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
    std::vector<float> decimatidData;

    fftwf_complex* fftw_in;
    fftwf_complex* fftw_out;
    fftwf_plan fftw_p;

    MainWindow* window = nullptr;
    std::thread* guiTh;

    float* i1;
    float* q1;

    int wideSpectorDataSamples = 0;
    int decimatidDataSamples = 0;
    int channeslDataSamples = 0;
    bool invers = true;

    void channelize();
    void decimation();
    void findPreambles(int ch);
    void dumpToFile(std::string fileName, void* data, size_t size);
    int bitDif(uint64_t, uint64_t);
    uint64_t toDifCoding(uint64_t);
    uint64_t toFromCoding(uint64_t);

    void freqCounter();
    void guiThread(int argc, char *argv[]);
};

#endif // DEMODULATOR_H
