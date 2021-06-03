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
#include <QApplication>
#include <volk/volk.h>

//#define QT_THREAD
//#define WINDOW

#ifdef WINDOW
    #include "mainwindow.h"
#endif

#ifdef QT_THREAD
    #include <QThread>
#endif

//#include "fdacoefs.h"
#include "fdacoefs_125K_to_100.h"
#include "fdacoefs_1M_to_125K.h"
#include "preamblepoint.h"
#include "decoder.h"
#include "GNU_UDP_client.hpp"

#define DATA_LEN            256
#define MAX_ERRORS          1
#define PREAMBLE_LEN        32
#define MAX_NOISE_NUM       1000
#define SUPER_X             2

#define TABLES_NUM      50
#define FFT_SIZE        (125000 / 100)

#ifdef QT_THREAD
class OpenUNBDemodulator : public QThread
#else
class OpenUNBDemodulator
#endif
{

#ifdef QT_THREAD
    Q_OBJECT
#endif

public:
    OpenUNBDemodulator(int inSempleRate, int xChan = 2, int id = 0);
    ~OpenUNBDemodulator();

    void addIQ(void* data, int sizeInSamples);
    void setCallback(void (*clb_f)(uint8_t* data, size_t size));
private:
    int numOfChannels;
    int decimationK;
    int id;

#ifdef QT_THREAD
    void run();
#endif

#ifdef WINDOW
    MainWindow* window = nullptr;
    std::thread* guiTh;

    void guiThread(int argc, char *argv[]);
#endif


    udp_client* udp;
    udp_client* udp2;

    unsigned int X_CHANNELS = 4;
    unsigned int CHANNEL_SAMPLE_RATE = 100;
    unsigned int SAMPLE_RATE = 125000;
    unsigned int SYM_RATE = 100;

    unsigned int SYM_LEN = (CHANNEL_SAMPLE_RATE/SYM_RATE);
    unsigned int CHANNELS = (SAMPLE_RATE / CHANNEL_SAMPLE_RATE * X_CHANNELS);

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
    std::complex<float>* iq1;
    std::complex<float>* iq2;

    float* i2, *q2;
    float* corr1, *corr2;

    uint32_t* corr;
    uint32_t* inv_corr;

    std::vector<PreamblePoint*> PreamblePointWithoutFullData;
    std::vector<PreamblePoint*> PreamblePointFullData;

    int wideSpectorDataSamples = 0;
    int decimatedDataSamples = 0;
    int channeslDataSamples = 0;
    bool invers = false;
    int shift = 0;
    int decimationPoint = 0;

#ifdef QT_THREAD
    bool isRunning = false;
#endif

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
    float** noise;
    float avg(float*, size_t);

    std::string uint32ToSring(uint32_t);

    const static int udpSemplSize = 100;
    fftwf_complex udpData[udpSemplSize];
    int it = 0;

    float** filterBank;
    int supX = 0;
};

#endif // DEMODULATOR_H
