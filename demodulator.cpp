#include "demodulator.h"

//#define DUMP_TO_FILE


OpenUNBDemodulator::OpenUNBDemodulator(int inSempleRate, int xChan, int _id) {
    X_CHANNELS = xChan * SUPER_X;
    CHANNELS = (SAMPLE_RATE / CHANNEL_SAMPLE_RATE * X_CHANNELS);
    numOfChannels = CHANNELS;
    decimationK = inSempleRate / SAMPLE_RATE;

    std::cout << "Demodulator init. Num of channels: " << numOfChannels << ". decimationK = " << decimationK << std::endl;

    channeslData = new std::vector<std::complex<float>>[numOfChannels];
    fftw_in = (fftwf_complex*) fftwf_malloc (numOfChannels / SUPER_X * sizeof(fftwf_complex));
    fftw_out = (fftwf_complex*) fftwf_malloc (numOfChannels / SUPER_X * sizeof(fftwf_complex));
    fftw_p = fftwf_plan_dft_1d(numOfChannels / SUPER_X, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);

    iq1 = new std::complex<float>[CHANNELS];
    iq2 = new std::complex<float>[CHANNELS];

    i2 = new float[CHANNELS];
    q2 = new float[CHANNELS];
    corr1 = new float[CHANNELS];
    corr2 = new float[CHANNELS];

    corr = new uint32_t[CHANNELS];
    inv_corr = new uint32_t[CHANNELS];

    noise = new double*[CHANNELS];
    for (int i=0; i<CHANNELS; i++) {
        noise[i] = new double[MAX_NOISE_NUM];
    }

    filterBank = new float*[CHANNELS / SUPER_X];
    for (int i=0; i<CHANNELS / SUPER_X; i++) {
        filterBank[i] = new float[BL_125K_to_100 / CHANNELS * SUPER_X];

        for (int j=0; j<BL_125K_to_100 / CHANNELS * SUPER_X; j++) {
            filterBank[i][j] = B_125K_to_100[i + j * CHANNELS / SUPER_X];
        }
    }


    //speedTh = new std::thread(&OpenUNBDemodulator::freqCounter, this);

    i1 = new float[numOfChannels];
    q1 = new float[numOfChannels];

    dec = new OpenUNBDecoder(SYM_LEN);

    id = _id;


#ifdef WINDOW
    guiTh = new std::thread(&OpenUNBDemodulator::guiThread, this, 0, nullptr);
#endif

    udp = new udp_client("192.168.159.1", 3333);
    udp2 = new udp_client("192.168.159.1", 3334);
    //mw = new MainWindow();
}

OpenUNBDemodulator::~OpenUNBDemodulator() {
    fftwf_free(fftw_in);
    fftwf_free(fftw_out);
    fftwf_destroy_plan(fftw_p);
    delete[] i1;
    delete[] q1;
    delete dec;

    delete[] iq1;
    delete[] iq2;

    delete[] i2;
    delete[] q2;
    delete[] corr1;
    delete[] corr2;

    delete[] corr;
    delete[] inv_corr;

    for (int i=0; i<CHANNELS; i++) {
        delete[] noise[i];
    }
    delete[] noise;

    for (int i=0; i<CHANNELS; i++) {
        delete[] filterBank[i];
    }
    delete[] filterBank;

}


void OpenUNBDemodulator::setCallback(void (*clb_f)(uint8_t* data, size_t size)) {
    dec->setCallback(clb_f);
}



void OpenUNBDemodulator::addIQ(void* data, int size) {
    struct timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    float* dataPtr = (float*)data;

    wideSpectorData.insert(wideSpectorData.end(), dataPtr, dataPtr + size*2);
    wideSpectorDataSamples += size;

#ifdef QT_THREAD
    while(isRunning);
    isRunning = true;
    start();
}

void OpenUNBDemodulator::run() {
#endif

    decimation();
    channelize();

    if (channeslData[0].size() > 100) {
        for (int i=0; i<PreamblePointWithoutFullData.size(); i++) {
            if (PreamblePointWithoutFullData[i]->data.size() >= (DATA_LEN + PREAMBLE_LEN + 1) * SYM_LEN) {
                //PreamblePointFullData.push_back(PreamblePointWithoutFullData[i]);
                dec->pushPreablePoint(PreamblePointWithoutFullData[i]);
                PreamblePointWithoutFullData.erase(PreamblePointWithoutFullData.begin() + i);
                i--;
            } else {
                PreamblePointWithoutFullData[i]->data.insert(PreamblePointWithoutFullData[i]->data.end(), channeslData[PreamblePointWithoutFullData[i]->channel].begin() + PREAMBLE_LEN, channeslData[PreamblePointWithoutFullData[i]->channel].end());
            }
        }

        int pr = 0;

        for (int i=0; i<numOfChannels; i++) {
            pr += findPreambles(i);
        }

#ifdef WINDOW
        if (window)
            window->setPreamble(dec->size());
#endif
        //if (pr)
        //    std::cout << id << ") Preambles: " << pr << "(" << dec->size() << ")" << std::endl;
    }
#ifdef QT_THREAD
    isRunning = false;

#endif
    clock_gettime(CLOCK_REALTIME, &end);

#ifdef WINDOW
    double time = ((end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec)/1e9) * 1e6;
    if (window)
        window->addProcessTime(time / timeout * 100);
    //std::cout << "Processing time: " << time << " ns (" << time / maxTimeout * 100 << " %) " << std::endl;
#endif
}

void OpenUNBDemodulator::decimation() {
    fftwf_complex* dataIn = (fftwf_complex*)wideSpectorData.data();
    lv_32fc_t dataToFir[BL_1M_to_125K];
    lv_32fc_t dataFromFir[BL_1M_to_125K];

    int i;

    for (i=0; i < (int)(wideSpectorData.size()/2) - BL_1M_to_125K - decimationK; i += decimationK) {
        float iC = 0.0;
        float qC = 0.0;

        //memcpy(dataToFir, dataIn[i], BL_1M_to_125K * sizeof(lv_32fc_t));
        volk_32fc_32f_multiply_32fc(dataFromFir, (lv_32fc_t*)dataIn[i], B_1M_to_125K, BL_1M_to_125K);
        for (int j=0; j<BL_1M_to_125K; j++) {
            iC += dataFromFir[j].real();
            qC += dataFromFir[j].imag();
        }

        decimatedData.push_back(iC);
        decimatedData.push_back(qC);

        decimatedDataSamples++;
    }

    if (i !=0) {
#ifdef DUMP_TO_FILE
        char p[256];
        sprintf(p, "dump/1M.complex", i);
        FILE* f = fopen(p, "a+");
        fwrite(wideSpectorData.data(), sizeof(float), i*2, f);
        fclose(f);
#endif
        wideSpectorData.erase(wideSpectorData.begin(), wideSpectorData.begin() + i*2);
    }

}



void OpenUNBDemodulator::channelize() {
    lv_32fc_t* dataIn = (lv_32fc_t*)decimatedData.data();
    lv_32fc_t dataToFir[BL_125K_to_100/numOfChannels];
    lv_32fc_t dataFromFir[BL_125K_to_100/numOfChannels];

    int i;

    int quarter = numOfChannels / X_CHANNELS;
    int iQuarter = 0;

    for (i=0; i < (int)(decimatedData.size()/2) - 2 * BL_125K_to_100 - numOfChannels/2; i += numOfChannels / X_CHANNELS / SUPER_X) {
        for (unsigned int ch=0; ch < numOfChannels/SUPER_X; ch++) {
            float iC = 0.0;
            float qC = 0.0;


            for (int j=0; j < BL_125K_to_100/numOfChannels*SUPER_X; j++) {
                iC += dataIn[i + ch + j*numOfChannels/SUPER_X].real() * filterBank[ch][j];
                qC += dataIn[i + ch + j*numOfChannels/SUPER_X].imag() * filterBank[ch][j];
            }

            /*
            for (int j=0; j < BL_125K_to_100/numOfChannels; j++) {
                dataToFir[j] = dataIn[i + ch + j*numOfChannels];
            }

            volk_32fc_32f_multiply_32fc(dataFromFir, dataToFir, filterBank[ch], BL_125K_to_100/numOfChannels);
            for (int j=0; j<BL_125K_to_100/numOfChannels; j++) {
                iC += dataFromFir[j].real();
                qC += dataFromFir[j].imag();
            }*/

            int iter = (numOfChannels / X_CHANNELS * shift + ch) % (numOfChannels / SUPER_X);
            fftw_in[iter][0] = iC;
            fftw_in[iter][1] = qC;

        }

        fftwf_execute(fftw_p);
        supX = (supX + 1) % SUPER_X;

#ifdef WINDOW
        int chOut = 0;

        if (window) {
            chOut = window->getSelectedChannel();
        }

        if (supX == chOut % SUPER_X) {
            udpData[it][0] = fftw_out[chOut/SUPER_X][0];
            udpData[it][1] = fftw_out[chOut/SUPER_X][1];
            it++;

            if (it == udpSemplSize) {
                it = 0;
                udp->send((char*)udpData, sizeof(udpData));
            }
        }
#endif


        channeslDataSamples++;

#ifdef DUMP_TO_FILE
        for(int i = 0; i < 20; i++) {
            if (i % SUPER_X == supX) {
                char p[256];
                sprintf(p, "dump/ch%d.complex", i);
                FILE* f = fopen(p, "a+");
                fwrite(fftw_out[i / SUPER_X], sizeof(fftwf_complex), 1, f);
                fclose(f);
            }
        }
#endif

        for (int ch=0; ch < numOfChannels / SUPER_X; ch++) {
            std::complex<float> c(fftw_out[ch][0], fftw_out[ch][1]);
            //noise += fftw_out[ch][0] * fftw_out[ch][0] + fftw_out[ch][1] * fftw_out[ch][1];
            //noise[ch * SUPER_X + supX][noiseNum] = fftw_out[ch][0] * fftw_out[ch][0] + fftw_out[ch][1] * fftw_out[ch][1];

            channeslData[ch * SUPER_X + supX].push_back(c);
        }

        noise[CHANNELS/2][noiseNum] = fftw_out[numOfChannels / SUPER_X/2][0] * fftw_out[numOfChannels / SUPER_X/2][0] + fftw_out[numOfChannels / SUPER_X/2][1] * fftw_out[numOfChannels / SUPER_X/2][1];
        invers = !invers;

        shift = (shift + 1) % X_CHANNELS;

        if (shift == 0) {
            noiseNum = (noiseNum + 1) % MAX_NOISE_NUM;
#ifdef WINDOW
            if (window)
                window->push1MHzData(fftw_out, numOfChannels / SUPER_X);
#endif
        }
    }

    if (i != 0) {
#ifdef DUMP_TO_FILE
        char p[256];
        sprintf(p, "dump/125K.complex", i);
        FILE* f = fopen(p, "a+");
        fwrite(decimatedData.data(), sizeof(float), i*2, f);
        fclose(f);
#endif

#ifdef WINDOW
        udp2->send((char*)decimatedData.data(), std::min(32*1024, i*2 * 4));
#endif
        decimatedData.erase(decimatedData.begin(), decimatedData.begin() + i*2);
    }

}

int OpenUNBDemodulator::bitDif(uint32_t a, uint32_t b) {
    uint32_t n = a ^ b;
    n -= (n >> 1) & 0x55555555;
    n = ((n >> 2) & 0x33333333) + (n & 0x33333333);
    n = ((((n >> 4) + n) & 0x0F0F0F0F) * 0x01010101) >> 24;
    return n;

}

//#define FP_LOG
//#define N_C

int OpenUNBDemodulator::findPreambles(int ch) {
    int ret = 0;
    int x;
    bool b = false;
    int bcount = 0;

    for (x=32; x < channeslData[ch].size(); x += SYM_LEN) {
        iq2[ch] = channeslData[ch][x];

        std::complex<float> diff(iq2[ch] * iq1[ch]);
        corr1[ch] = diff.imag();
        corr2[ch] = diff.real();

        corr[ch] = (corr[ch] >> 1) | ((corr1[ch] > 0 ? 1:0) << 31);
        inv_corr[ch] = (inv_corr[ch] >> 1) | ((corr2[ch] > 0 ? 1:0) << 31);

        iq1[ch] = std::complex<float>(iq2[ch].real(), -iq2[ch].imag());

        const uint32_t prea = 0x97157A6F;

        int err1 =  bitDif(corr[ch], prea);
        int err2 =  bitDif(inv_corr[ch], (~prea));

        int err3 =  bitDif(corr[ch], ~prea);
        int err4 =  bitDif(inv_corr[ch], (prea));

        double noiseD = 0.0;
        if (err1 <= MAX_ERRORS || err2 <= MAX_ERRORS || err3 <= MAX_ERRORS || err4 <= MAX_ERRORS) {
            b = true;
            noiseD = avg(noise[numOfChannels/2], MAX_NOISE_NUM);
        }

#ifdef FP_LOG
        if (ch == 20 && b) {
            std::cout << "Preamble!!!" << std::endl;
            std::cout << " Channel: " << ch << std::endl;
            std::cout << " Errors: " << err1 << std::endl;
            std::cout << " Data: " << uint32ToSring(corr[ch]) << std::endl;
            std::cout << " Expt: " << uint32ToSring(prea) << std::endl;
            float qAvg = 0;
            float iAvg = 0;
            for (int i=0; i<32; i++) {
                qAvg += abs(channeslData[ch][x - i].imag());
                iAvg += abs(channeslData[ch][x - i].real());
            }
            std::cout << " qAvg: " << qAvg/iAvg << std::endl;
            b = false;
        }
#endif
        if (err1 <= MAX_ERRORS) {
            PreamblePoint* pp1 = new PreamblePoint;
            pp1->channel = ch;
            pp1->inv = 0;
            pp1->preableErrors = err1;
            pp1->noise = noiseD;
            pp1->pos = x + totalSamples;
            pp1->batch = totalBatches;
            pp1->data.insert(pp1->data.begin(), channeslData[ch].begin() + x - 32,channeslData[ch].begin() + x - 32 + std::min(channeslData[ch].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(pp1);
            ret++;

#ifdef N_C
            PreamblePoint* ppn = new PreamblePoint(pp1);
            ppn->channel = ch - 1;
            PreamblePoint* ppp = new PreamblePoint(pp1);
            ppp->channel = ch + 1;
            ppn->data.insert(ppn->data.begin(), channeslData[ch-1].begin() + x - 32,channeslData[ch-1].begin() + x - 32 + std::min(channeslData[ch-1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppn);
            ppp->data.insert(ppp->data.begin(), channeslData[ch+1].begin() + x - 32,channeslData[ch+1].begin() + x - 32 + std::min(channeslData[ch+1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppp);
#endif
        }

        if (err2 <= MAX_ERRORS) {
            PreamblePoint* pp2 = new PreamblePoint;
            pp2->channel = ch;
            pp2->inv = PreamblePoint::invIQ;
            pp2->preableErrors = err2;
            pp2->noise = noiseD;
            pp2->pos = x + totalSamples;
            pp2->batch = totalBatches;
            pp2->data.insert(pp2->data.begin(), channeslData[ch].begin() + x - 32,channeslData[ch].begin() + x - 32 + std::min(channeslData[ch].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(pp2);
            ret++;

#ifdef N_C
            PreamblePoint* ppn = new PreamblePoint(pp2);
            ppn->channel = ch - 1;
            PreamblePoint* ppp = new PreamblePoint(pp2);
            ppp->channel = ch + 1;
            ppn->data.insert(ppn->data.begin(), channeslData[ch-1].begin() + x - 32,channeslData[ch-1].begin() + x - 32 + std::min(channeslData[ch-1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppn);
            ppp->data.insert(ppp->data.begin(), channeslData[ch+1].begin() + x - 32,channeslData[ch+1].begin() + x - 32 + std::min(channeslData[ch+1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppp);
#endif
        }

        if (err3 <= MAX_ERRORS) {
            PreamblePoint* pp3 = new PreamblePoint;
            pp3->channel = ch;
            pp3->inv = PreamblePoint::invBits;
            pp3->preableErrors = err3;
            pp3->noise = noiseD;
            pp3->pos = x + totalSamples;
            pp3->batch = totalBatches;
            pp3->data.insert(pp3->data.begin(), channeslData[ch].begin() + x - 32,channeslData[ch].begin() + x - 32 + std::min(channeslData[ch].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(pp3);

#ifdef N_C
            PreamblePoint* ppn = new PreamblePoint(pp3);
            ppn->channel = ch - 1;
            PreamblePoint* ppp = new PreamblePoint(pp3);
            ppp->channel = ch + 1;
            ppn->data.insert(ppn->data.begin(), channeslData[ch-1].begin() + x - 32,channeslData[ch-1].begin() + x - 32 + std::min(channeslData[ch-1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppn);
            ppp->data.insert(ppp->data.begin(), channeslData[ch+1].begin() + x - 32,channeslData[ch+1].begin() + x - 32 + std::min(channeslData[ch+1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppp);
#endif
            ret++;
        }

        if (err4 <= MAX_ERRORS) {
            PreamblePoint* pp4 = new PreamblePoint;
            pp4->channel = ch;
            pp4->inv = PreamblePoint::invIQ | PreamblePoint::invBits;
            pp4->preableErrors = err4;
            pp4->noise = noiseD;
            pp4->pos = x + totalSamples;
            pp4->batch = totalBatches;
            pp4->data.insert(pp4->data.begin(), channeslData[ch].begin() + x - 32,channeslData[ch].begin() + x - 32 + std::min(channeslData[ch].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(pp4);

#ifdef N_C
            PreamblePoint* ppn = new PreamblePoint(pp4);
            ppn->channel = ch - 1;
            PreamblePoint* ppp = new PreamblePoint(pp4);
            ppp->channel = ch + 1;
            ppn->data.insert(ppn->data.begin(), channeslData[ch-1].begin() + x - 32,channeslData[ch-1].begin() + x - 32 + std::min(channeslData[ch-1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppn);
            ppp->data.insert(ppp->data.begin(), channeslData[ch+1].begin() + x - 32,channeslData[ch+1].begin() + x - 32 + std::min(channeslData[ch+1].size() - x + 32, (size_t)(DATA_LEN + PREAMBLE_LEN + 1)));
            PreamblePointWithoutFullData.push_back(ppp);
#endif
            ret++;


        }

    }
    if (x != 32)
        channeslData[ch].erase(channeslData[ch].begin(), channeslData[ch].end() - 32);

    return ret;
}

void OpenUNBDemodulator::freqCounter() {
    const int m = 7;
    int lastWD = 0, lastDD = 0, lastCD = 0;
    int avgWD[m], avgDD[m], avgCD[m];
    int _avgWD, _avgDD, _avgCD;
    int index = 0;
    while (1) {
        avgWD[index] = wideSpectorDataSamples - lastWD;
        avgDD[index] = decimatedDataSamples - lastDD;
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
        lastDD = decimatedDataSamples;
        lastCD = channeslDataSamples;

        usleep(1000000);
    }
}

std::string OpenUNBDemodulator::uint32ToSring(uint32_t a) {
    std::string ret = "";

    for (int i=0; i<sizeof(uint32_t) * 8; i++) {
        ret.push_back((a >> (sizeof(uint32_t) * 8 - i - 1)) & 1 ? '1' : '0');
    }

    return ret;
}

double OpenUNBDemodulator::avg(double* data, size_t size) {

    double ret = 0;
    for (size_t i=0; i<size; i++) {
        ret += data[i];
    }
    //volk_32fc_
    //volk_32f_s32f_calc_spectral_noise_floor_32f

    return ret/size;
}

void OpenUNBDemodulator::dumpToFile(std::string fileName, void *data, size_t size) {
    std::fstream file(fileName);
    if (file.is_open()) {
        file.write((char*)data, size);
        file.close();
    }
}

#ifdef WINDOW
void OpenUNBDemodulator::guiThread(int argc, char *argv[]) {
    QApplication a(argc, argv);
    window = new MainWindow(CHANNELS/SUPER_X);
    window->show();
    a.exec();

    exit(0);
}
#endif

