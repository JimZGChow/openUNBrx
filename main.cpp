#include "SoapyEnum.h"
#include "SDRDevInfo.h"
#include "demodulator.h"
#include "udp_reciever.h"
#include "preamblepoint.h"
#include "GNU_UDP_client.hpp"

#include <algorithm>
#include <iostream>
#include <thread>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fftw3.h>
#include <time.h>
#include <csignal>

//#define FILE

int g_exitRecvThread = 0;

int count = 0;

void ctrlC(int) {
    std::cout << "Exit" << std::endl;
    g_exitRecvThread = 1;
    //exit (0);
}


void recvData(SoapySDR::Device* dev, SoapySDR::Stream* stream, std::vector<char>* recvBuf, int bytes_per_sample) {
    dev->activateStream(stream);
    int numElems = recvBuf->size() / bytes_per_sample;
    void* Buffs[1];
    int flags;
    long long timeNs;
    long maxTimeout = numElems / dev->getSampleRate(SOAPY_SDR_RX, 0) * 1000 * 1000;

    Buffs[0] = (void*)recvBuf->data();
    OpenUNBDemodulator dem(dev->getSampleRate(SOAPY_SDR_RX, 0), 4, 4);

    struct timespec begin, end;
    while(!g_exitRecvThread){
        int n_stream_read = dev->readStream(stream, Buffs, numElems, flags, timeNs, maxTimeout);

        if(n_stream_read < 0)
            std::cout << " Soapy read failed with code: " << n_stream_read << std::endl;
        else {
            //std::cout << " Read " << n_stream_read << std::endl;

            clock_gettime(CLOCK_REALTIME, &begin);
            dem.addIQ(Buffs[0], n_stream_read);
            clock_gettime(CLOCK_REALTIME, &end);

            double time = ((end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec)/1e9) * 1e6;

            std::cout << "Processing time: " << time << " ns (" << time / maxTimeout * 100 << " %) " << std::endl;
        }
    }

    dev->deactivateStream(stream);
}

int main(int argc, char *argv[]) {

    std::signal(SIGINT, ctrlC);


#ifndef FILE
    double freq = 868000000.0;

    if (argc == 4)
        freq = (double) atoi(argv[3]);

    SoapyEnum sdr_enum;

    std::string sdr_driver("hackrf");

    std::vector<SDRDevInfo*>* sdrDevices = sdr_enum.enumerateDevices();

    if(sdrDevices != 0) {

        std::vector<std::string> factories = sdr_enum.getFactories();

        std::cout << factories.size() << " devices" << std::endl;

        if(std::find(factories.begin(), factories.end(), sdr_driver) != factories.end())
        {
            std::cout << "Found " << sdr_driver << " factory " << std::endl;

            auto result = std::find_if(sdrDevices->begin(), sdrDevices->end(), [sdr_driver](SDRDevInfo* dev_i) { return dev_i->getDriver() == sdr_driver;});
            if(result == sdrDevices->end())
                std::cout << "Device " << sdr_driver << " not found " << std::endl;
            else
                std::cout << "Device " << sdr_driver << " found! " << std::endl;

            // Make RTL Device

            SDRDevInfo* devInfo = (*result);
            SoapySDR::Device* soapyDev = devInfo->getSoapyDevice();

            if(soapyDev)
            {
                SoapySDR::RangeList freqRange = soapyDev->getFrequencyRange(SOAPY_SDR_RX, 0);
                SoapySDR::RangeList sampleRateRange = soapyDev->getSampleRateRange(SOAPY_SDR_RX, 0);

                double sampleRate = soapyDev->getSampleRate(SOAPY_SDR_RX, 0);

                //std::cout << freqRange

                std::cout << "Got center frequency: " << freq / 1000000.0 << " MHz" << std::endl;
                std::cout << "Got sample rate: " << sampleRate / 1000.0 << " KHz" << std::endl;

                std::vector<std::string> formats = soapyDev->getStreamFormats(SOAPY_SDR_RX, 0);
                std::cout << "Supported formats: ";
                for(auto fmt : formats)
                    std::cout << fmt << ", ";
                std::cout << std::endl;

                int bytes_per_sample = 4;
                std::string fmt;
#ifdef CS16
                if(std::find(formats.begin(), formats.end(), "CS16") != formats.end()) {
                    fmt = "CS16";
                    bytes_per_sample = sizeof(short)*2;
                }
                else
#endif
                if(std::find(formats.begin(), formats.end(), "CF32") != formats.end()) {
                    fmt = "CF32";
                    bytes_per_sample = sizeof(float)*2;
                }
                else {
                    fmt = formats.back();
                    bytes_per_sample = sizeof(float)*2;
                }


                SoapySDR::Kwargs args = devInfo->getDeviceArgs();
                SoapySDR::Stream* pStream = soapyDev->setupStream(SOAPY_SDR_RX, fmt, std::vector<size_t>(), args);
                if(!pStream) {
                    std::cout << "Error setup stream !!! Exit.. " << std::endl;
                    exit(EXIT_FAILURE);
                }

                std::cout << "Format: " << fmt << std::endl;

                // Set new frequancy and samplerate
                std::cout << "Max frec.: " << freqRange.back().maximum() << std::endl;
                std::cout << "Min frec.: " << freqRange.front().minimum() << std::endl;

                std::cout << "Max sample frec.: " << sampleRateRange.back().maximum() << std::endl;
                std::cout << "Min sample frec.: " << sampleRateRange.front().minimum() << std::endl;

                //freq = 868000000.0;
                sampleRate = 1000000.0;

                soapyDev->setFrequency(SOAPY_SDR_RX, 0, freq);
                soapyDev->setSampleRate(SOAPY_SDR_RX, 0, sampleRate);

                soapyDev->setGain(SOAPY_SDR_RX, 0, 35.0);
                soapyDev->setGainMode(SOAPY_SDR_RX, 0, true);
                soapyDev->setDCOffsetMode(SOAPY_SDR_RX, 0, false);
                soapyDev->writeSetting("digital_agc", "true");
                soapyDev->writeSetting("direct_samp", "0");
                soapyDev->writeSetting("offset_tune", "true");

                std::cout << "Gain = " << soapyDev->getGain(SOAPY_SDR_RX, 0) << std::endl;

                std::cout << "Setting new center freq: " << soapyDev->getFrequency(SOAPY_SDR_RX, 0) / 1000000.0 << "MHz" << " new samplerate: " << soapyDev-> getSampleRate (SOAPY_SDR_RX, 0) / 1000 << "KHz" << std::endl;

                int streamMTU = soapyDev->getStreamMTU(pStream);

                std::cout << "streamMTU = " << soapyDev->getStreamMTU(pStream) << std::endl;

                std::vector<char> recvBuf(streamMTU * bytes_per_sample);

                std::cout << "Starting receiver thread.... " << std::endl;

                //std::thread t(recvThread, soapyDev, pStream, &recvBuf, bytes_per_sample);

                std::cout << "Press ENTER to exit receiving thread" << std::endl;
                recvData( soapyDev, pStream, &recvBuf, bytes_per_sample);

                soapyDev->closeStream(pStream);
            }
        }
    }
#else
    OpenUNBDemodulator dem(1000000, 2);

    std::fstream f("/home/def/workspace/OpenUNB/rxMy/build-OpenUNB_GUI_rx-Desktop_Qt_5_15_1_GCC_64bit-Debug/dump_s/1M.complex");

    if ( !f.is_open()) {
        return -1;
    }
    float k=0;

    f.seekg(0, std::ios_base::end);
    int e = f.tellg();
    f.seekg(0, std::ios_base::beg);

    std::cout << "File size: " << e << std::endl;
    std::cout << "File time: " << e/sizeof(std::complex<float>)/1000000 << " s" << std::endl;

    long startSeek = 0 * sizeof (fftwf_complex);
    long long seek = startSeek;//;1000000 * 2 * sizeof(float) + 4000000;
    long long step = 1024 * 128 * 2 * sizeof(float) * 1;
    //long long step = e;
    char* data = new char[step];
    memset(data, 1, step);

    const double mean = 0.0;
    double stddev = 0.5;
    std::default_random_engine generator;
    std::normal_distribution<float> dist(mean, 1.0);

    struct timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);
    clock_gettime(CLOCK_REALTIME, &end);

    int sr = 1000000;
    double fr = 25;
    double div = fr/sr;
    double delta = 0.0f;

    float c = -1.0;

    clock_gettime(CLOCK_REALTIME, &begin);
    while(!g_exitRecvThread){
        int n_stream_read;

        f.seekg(seek, std::ios_base::beg);
        //std::cout << " pos " << f.tellg() << std::endl;
        f.read(data, step);

        float* dataf = (float*)data;
        std::complex<float>* datac = (std::complex<float>*)data;
        for (int i=0; i<step / (sizeof(float)); i+=2) {
            std::complex<float> mul(cos(delta), sin(delta));
            //dataf[i] += dist(generator) * stddev;
            //dataf[i] = c;
            //dataf[i + 1] = c;
            //c += 0.01;

            //if (c > 1)
            //    c = -1.0;
            //datac[i/2] = mul;
            delta += div * 2 * M_PI;
        }

        dem.addIQ(data, step / (sizeof(float) * 2));
        //usleep(step / (sizeof(float) * 2));
        seek += step;
        if (seek > e - step) {
            clock_gettime(CLOCK_REALTIME, &end);
            double time = ((end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec)/1e9);
            seek = startSeek;

            //fr += 1;
            //div = fr/sr;
            stddev += 0.01;
            if (stddev > 1.2)
                stddev = 0.01;
            std::cout << "======= Time: " << time << " sec" << std::endl;
            clock_gettime(CLOCK_REALTIME, &begin);
        }
        //std::cout << "============" << std::endl;
    }

    //return a.exec();
    return 0;
#endif

    //return a.exec();
    return 0;
}
