#ifndef DEMODULATOR_H
#define DEMODULATOR_H

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
#include "fdacoefs_1M_to_125K.h"
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
#define MAX_ERRORS          2
#define PREAMBLE_LEN        32
#define MAX_NOISE_NUM       100

class OpenUNBDemodulator
{
public:
    OpenUNBDemodulator(int inSempleRate);
    ~OpenUNBDemodulator();

    void addIQ(void* data, int sizeInSamples);
    void setCallback(void (*clb_f)(uint8_t* data, size_t size));
private:
    int numOfChannels;
    int decimationK;

    std::vector<float> wideSpectorData;
    std::vector<std::complex<float>>* channeslData;
    std::vector<std::complex<float>> udpSendData;
    std::vector<float> decimatedData;

    fftwf_complex* fftw_in;
    fftwf_complex* fftw_out;
    fftwf_plan fftw_p;

    std::thread* speedTh;

    float* i1;
    float* q1;
    std::complex<float> iq1[CHANNELS];
    std::complex<float> iq2[CHANNELS];

    float i2[CHANNELS], q2[CHANNELS];
    float corr1[CHANNELS], corr2[CHANNELS];

    uint32_t corr[CHANNELS];
    uint32_t inv_corr[CHANNELS];

    std::vector<PreamblePoint*> PreamblePointWithoutFullData;
    std::vector<PreamblePoint*> PreamblePointFullData;

    int wideSpectorDataSamples = 0;
    int decimatedDataSamples = 0;
    int channeslDataSamples = 0;
    bool invers = true;

    OpenUNBDecoder* dec;

    void channelize();
    void decimation();
    int findPreambles(int ch);
    void dumpToFile(std::string fileName, void* data, size_t size);
    int bitDif(uint32_t, uint32_t);
    uint64_t toFromCoding(uint64_t);

    void freqCounter();

    unsigned int noiseNum = 0;
    unsigned int totalSamples = 0;
    unsigned int totalBatches = 0;
    float noise[CHANNELS][MAX_NOISE_NUM];
    float avg(float*, size_t);

    std::string uint32ToSring(uint32_t);
};

#endif // DEMODULATOR_H
