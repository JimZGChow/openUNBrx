#include "demodulator.h"

Demodulator::Demodulator(int inSempleRate) {
    numOfChannels = SAMPLE_RATE / CHANNEL_SAMPLE_RATE * X2;
    decimationK = inSempleRate / SAMPLE_RATE;

    std::cout << "Demodulator init. Num of channels: " << numOfChannels << ". decimationK = " << decimationK << std::endl;

    channeslData = new std::vector<std::complex<float>>[numOfChannels];
    fftw_in = (fftwf_complex*) fftwf_malloc (numOfChannels * sizeof(fftwf_complex));
    fftw_out = (fftwf_complex*) fftwf_malloc (numOfChannels * sizeof(fftwf_complex));
    fftw_p = fftwf_plan_dft_1d(numOfChannels, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);

    i1 = new float[numOfChannels];
    q1 = new float[numOfChannels];

    guiTh = new std::thread(&Demodulator::guiThread, this, 0, nullptr);
}

Demodulator::~Demodulator() {
    fftwf_free(fftw_in);
    fftwf_free(fftw_out);
    fftwf_destroy_plan(fftw_p);
    delete[] i1;
    delete[] q1;
}

void Demodulator::addIQ(void* data, int size) {
    //printf("addIQ %d\n", size);
    float* dataPtr = (float*)data;

    wideSpectorData.insert(wideSpectorData.end(), dataPtr, dataPtr + size*2);
    wideSpectorDataSamples += size;


//    if (DbS.isConnected()) {
//        DbS.sendData(data, 1024*sizeof(float) * 2, dataType::data_1_mhz);
//    }

    decimation();
    //if (decimatidData.size() > 250*1024)
    channelize();

    char tmp[numOfChannels * sizeof(float) * 2];
    memcpy(tmp, fftw_out + numOfChannels/2, sizeof(tmp) / 2);
    memcpy(tmp + sizeof(tmp) / 2, fftw_out, sizeof(tmp) / 2);

    if (window) {
        window->push1MHzData(fftw_out, numOfChannels);
    }

    if (channeslData[0].size() > 100) {
        for (int i=0; i<numOfChannels; i++) {
            //std::cout << "ch = " << i << std::endl;
            findPreambles(i);
            channeslData[i].clear();
        }
    }
}

float dds = 0;

void Demodulator::decimation() {
    fftwf_complex* dataIn = (fftwf_complex*)wideSpectorData.data();

    int i;

    for (i=0; i < (int)(wideSpectorData.size()/2) - BL_1M_to_125K - decimationK; i += decimationK) {
        float iC = 0.0;
        float qC = 0.0;

        for (int j=0; j<BL_1M_to_125K; j++) {
            iC += dataIn[i + j][0] * B_1M_to_125K[j];
            qC += dataIn[i + j][1] * B_1M_to_125K[j];
        }

//        fftwf_complex tmp;
//        tmp[0] = iC;
//        tmp[1] = qC;
        dds += 0.000001;
        if (dds > 1)
            dds = 0;
        //DbS.sendData(&tmp, numOfChannels * 2 * sizeof(float), data_100_hz);
        decimatidData.push_back(iC);
        decimatidData.push_back(qC);
        decimatidDataSamples++;
    }

    if (i !=0) {
        i -= decimationK;
        wideSpectorData.erase(wideSpectorData.begin(), wideSpectorData.begin() + i*2);
    }
    //printf("decimatidData %d\n", wideSpectorData.size());
}

float ou = 0;

void Demodulator::channelize() {
    fftwf_complex* dataIn = (fftwf_complex*)decimatidData.data();

    int i;

    //std::cout << (int)(decimatidData.size()/2) - BL_125K_to_100 - numOfChannels*BL_125K_to_100/numOfChannels << std::endl;

    for (i=0; i < (int)(decimatidData.size()/2) - BL_125K_to_100 - numOfChannels*BL_125K_to_100/numOfChannels - numOfChannels/2; i += numOfChannels/2) {
        for (unsigned int ch=0; ch < numOfChannels; ch++) {
            float iC = 0.0;
            float qC = 0.0;

            for (int j=0; j < BL_125K_to_100/numOfChannels; j++) {
                iC += dataIn[i + ch + j*numOfChannels][0] * B_125K_to_100[ch + j*numOfChannels];
                qC += dataIn[i + ch + j*numOfChannels][1] * B_125K_to_100[ch + j*numOfChannels];
            }


            if (invers) {
                if (ch < numOfChannels/2) {
                    fftw_in[ch + numOfChannels/2][0] = iC;
                    fftw_in[ch + numOfChannels/2][1] = qC;
                }
                else {
                    fftw_in[ch - numOfChannels/2][0] = iC;
                    fftw_in[ch - numOfChannels/2][1] = qC;
                }
            }
            else {
                fftw_in[ch][0] = iC;
                fftw_in[ch][1] = qC;
            }


            //fftw_in[ch][0] = iC;
            //fftw_in[ch][1] = qC;
        }

        fftwf_execute(fftw_p);

        channeslDataSamples++;

        if (window) {
            window->push100HzData(fftw_out[window->getSelectedChannel()]);
        }

        for (int ch=0; ch < numOfChannels; ch++) {
            std::complex<float> c(fftw_out[ch][0], fftw_out[ch][1]);

            channeslData[ch].push_back(c);
        }

        invers = !invers;
    }

    //printf("channeslData %d %d\n", i / (numOfChannels/2), channeslData[0].size() * numOfChannels);

    if (i != 0) {
        //i -= numOfChannels/2;
        decimatidData.erase(decimatidData.begin(), decimatidData.begin() + i*2);
    }
}

int Demodulator::bitDif(uint64_t a, uint64_t b) {
    int ret = 0;
    uint64_t c = a ^ b;
    for (int i=0; i<64; i++) {
        if ((c >> i) & 1)
            ret++;
    }

    return ret;
}

void Demodulator::findPreambles(int ch) {
    //float i1=0, q1=0;
    float i2=0, q2=0;
    float corr1 = 0, corr2 = 0;

    uint32_t corr = 0;
    uint32_t inv_corr = 0;

    for (int x=0; x < channeslData[ch].size(); x += SYM_LEN) {
        i2 = channeslData[ch][x].imag();
        q2 = channeslData[ch][x].real();

        corr1 = i1[ch]*i2 + q1[ch]*q2;
        corr2 = -i1[ch]*q2 + q1[ch]*i2;

        corr = (corr << 1) + (corr1 > 0);
        inv_corr = (inv_corr << 1) + (corr2 > 0);

        i1[ch] = i2;
        q1[ch] = q2;

        int err =  bitDif(corr & 0xFFFFFFFF, 0xFFFF5555);
        if (err < 3) {
            std::cout << "!!!!!!! "<< (corr & 0xFFFF) << " " << err << " " << ch << std::endl;
        }

        if ((inv_corr) == 0x76732422) {
            std::cout << "!!!!!!!2 " << inv_corr << " " << ch << std::endl;
        }
    }
}

void Demodulator::guiThread(int argc, char *argv[]) {
    QApplication a(argc, argv);
    window = new MainWindow();
    window->show();
    a.exec();

    exit(0);
}

uint64_t Demodulator::toDifCoding(uint64_t a) {
    uint64_t ret = 0;
    for (int i=64; i <= 1; i++) {
        ret = (ret << 1) | ((a >> i) & (a >> (i - 1)));
    }


}

/*
void Demodulator::freqCounter() {
    const int m = 25;
    int lastWD = 0, lastDD = 0, lastCD = 0;
    int avgWD[m], avgDD[m], avgCD[m];
    int _avgWD, _avgDD, _avgCD;
    int index = 0;
    while (1) {
        avgWD[index] = wideSpectorDataSamples - lastWD;
        avgDD[index] = decimatidDataSamples - lastDD;
        avgCD[index] = channeslDataSamples - lastCD;
        index = (index + 1) % m;

        _avgWD = _avgDD = _avgCD = 0;
        for (int i=0; i<m; i++) {
            _avgWD += avgWD[i];
            _avgDD += avgDD[i];
            _avgCD += avgCD[i];
        }

        std::cout << " WD = " << _avgWD/m << std::endl;
        std::cout << " DD = " << _avgDD/m << std::endl;
        std::cout << " CD = " << _avgCD/m << std::endl;

        lastWD = wideSpectorDataSamples;
        lastDD = decimatidDataSamples;
        lastCD = channeslDataSamples;

        usleep(1000000);
    }
}*/
