#include "SoapyEnum.h"
#include "SDRDevInfo.h"
#include "demodulator.h"
#include "pcscl_noperm.h"
#include "udp_reciever.h"
#include "preamblepoint.h"

#include <algorithm>
#include <iostream>
#include <thread>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fftw3.h>
#include <time.h>
#include <csignal>


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
    Demodulator dem(dev->getSampleRate(SOAPY_SDR_RX, 0));

    struct timespec begin, end;
    UDP_Reciever udp(recvBuf->size());
    //udp.readData((uint8_t*)Buffs[0]);
    std::fstream f("/home/pi/1M_openUNB.complex");

    float k=0;

    int seek = 1000000 * 2 * sizeof(float) + 4000000;
    int step = 1000000 * 2 * sizeof(float) * 1;
    char* data = new char[step];

    while(!g_exitRecvThread){
        int n_stream_read = dev->readStream(stream, Buffs, numElems, flags, timeNs, maxTimeout);
        //int n_stream_read = udp.readData((uint8_t*)Buffs[0]);
//        int n_stream_read;

//        f.seekg(seek, std::ios_base::beg);
//        std::cout << " pos " << f.tellg() << std::endl;
//        f.read(data, step);
//        dem.addIQ(data, step / (sizeof(float) * 2));
//        usleep(1000000);
//        //seek += step;
//        if (seek > 38000000 - step)
//            seek = 0;
//        continue;
//        float* f = (float*)Buffs[0];
//        for (int i=0; i<n_stream_read; i++) {
//            f[i*2] = k;
//            f[i*2 + 1] = k;

//            k += 0.001;
//            if (k > 1)
//                k = 0;
//        }

        if(n_stream_read < 0)
            std::cout << " Soapy read failed with code: " << n_stream_read << std::endl;
        else {
            //std::cout << " Read " << n_stream_read << std::endl;

            clock_gettime(CLOCK_REALTIME, &begin);
            dem.addIQ(Buffs[0], n_stream_read);
            clock_gettime(CLOCK_REALTIME, &end);

            double time = ((end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec)/1e9) * 1e6;

            //std::cout << "Processing time: " << time << " ns (" << time / maxTimeout * 100 << " %) " << std::endl;
        }
    }

    dev->deactivateStream(stream);
}

int main(int argc, char *argv[]) {

    std::signal(SIGINT, ctrlC);

    double freq = 868000000.0;

    if (argc == 4)
        freq = (double) atoi(argv[3]);

    SoapyEnum sdr_enum;

    std::string sdr_driver("rtlsdr");

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

                soapyDev->setGain(SOAPY_SDR_RX, 0, 15.0);
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

    //return a.exec();
    return 0;
}
