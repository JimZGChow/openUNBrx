#include "demodulator.h"

OpenUNBDemodulator::OpenUNBDemodulator(int inSempleRate) {
    numOfChannels = SAMPLE_RATE / CHANNEL_SAMPLE_RATE * X2;
    decimationK = inSempleRate / SAMPLE_RATE;

    std::cout << "Demodulator init. Num of channels: " << numOfChannels << ". decimationK = " << decimationK << std::endl;

    channeslData = new std::vector<std::complex<float>>[numOfChannels];
    fftw_in = (fftwf_complex*) fftwf_malloc (numOfChannels * sizeof(fftwf_complex));
    fftw_out = (fftwf_complex*) fftwf_malloc (numOfChannels * sizeof(fftwf_complex));
    fftw_p = fftwf_plan_dft_1d(numOfChannels, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);

    //speedTh = new std::thread(&OpenUNBDemodulator::freqCounter, this);

    i1 = new float[numOfChannels];
    q1 = new float[numOfChannels];

    dec = new OpenUNBDecoder(SYM_LEN);
}

OpenUNBDemodulator::~OpenUNBDemodulator() {
    fftwf_free(fftw_in);
    fftwf_free(fftw_out);
    fftwf_destroy_plan(fftw_p);
    delete[] i1;
    delete[] q1;
    delete dec;
}


void OpenUNBDemodulator::setCallback(void (*clb_f)(uint8_t* data, size_t size)) {
    dec->setCallback(clb_f);
}

void OpenUNBDemodulator::addIQ(void* data, int size) {

    float* dataPtr = (float*)data;

    //std::cout << size*2 * sizeof(float) << std::endl;

    wideSpectorData.insert(wideSpectorData.end(), dataPtr, dataPtr + size*2);
    wideSpectorDataSamples += size;

    //udp1233->send((char*)dataPtr, sizeof(float)*size*2);

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
            //std::cout << "ch = " << i << std::endl;
            //if (i >= 145 && i <= 160)
                pr += findPreambles(i);
            if (i == numOfChannels - 1) {
                totalSamples += channeslData[i].size();
                totalBatches++;
            }
            //channeslData[i].erase(channeslData[i].begin(), channeslData[i].end() - 32);
        }

        for (int i=0; i<numOfChannels; i++) {
            channeslData[i].erase(channeslData[i].begin(), channeslData[i].end() - 32);
        }

        if (pr)
            std::cout << "Preambles: " << pr << "(" << dec->size() << ")" << std::endl;
    }
}

void OpenUNBDemodulator::decimation() {

    fftwf_complex* dataIn = (fftwf_complex*)wideSpectorData.data();
    int i;

    for (i=0; i < (int)(wideSpectorData.size()/2) - BL_1M_to_125K - decimationK; i += decimationK) {
        float iC = 0.0;
        float qC = 0.0;

        for (int j=0; j<BL_1M_to_125K; j++) {
            iC += dataIn[i + j][0] * B_1M_to_125K[j];
            qC += dataIn[i + j][1] * B_1M_to_125K[j];
        }

        decimatedData.push_back(iC);
        decimatedData.push_back(qC);

        decimatedDataSamples++;
    }

    if (i !=0) {
        wideSpectorData.erase(wideSpectorData.begin(), wideSpectorData.begin() + i*2);
    }

}

void OpenUNBDemodulator::channelize() {
    fftwf_complex* dataIn = (fftwf_complex*)decimatedData.data();

    int i;

    for (i=0; i < (int)(decimatedData.size()/2) - BL_125K_to_100 - numOfChannels*BL_125K_to_100/numOfChannels - numOfChannels/2; i += numOfChannels/2) {
        //udp1233->send((char*)dataIn[i], 1024);

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
        }

        fftwf_execute(fftw_p);
        //udp1233->send((char*)fftw_in, sizeof(float)*numOfChannels*2);

        channeslDataSamples++;

        for (int ch=0; ch < numOfChannels; ch++) {
            std::complex<float> c(fftw_out[ch][0], fftw_out[ch][1]);
            //noise += fftw_out[ch][0] * fftw_out[ch][0] + fftw_out[ch][1] * fftw_out[ch][1];
            noise[ch][noiseNum] = fftw_out[ch][0] * fftw_out[ch][0] + fftw_out[ch][1] * fftw_out[ch][1];
            channeslData[ch].push_back(c);

        }
        //noise /= numOfChannels;
        noiseNum = (noiseNum + 1) % MAX_NOISE_NUM;

        invers = !invers;
    }

    if (i != 0) {
        decimatedData.erase(decimatedData.begin(), decimatedData.begin() + i*2);
    }
}

int OpenUNBDemodulator::bitDif(uint32_t a, uint32_t b) {
    int ret = 0;
    uint64_t c = a ^ b;
    for (int i=0; i<64; i++) {
        if ((c >> i) & 1)
            ret++;
    }

    return ret;
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
        if (err1 <= MAX_ERRORS || err2 <= MAX_ERRORS || err3 <= MAX_ERRORS || err4 <= MAX_ERRORS) {
            b = true;
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
            pp1->noise = avg(noise[ch], MAX_NOISE_NUM);
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
            pp2->noise = avg(noise[ch], MAX_NOISE_NUM);
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
            pp3->noise = avg(noise[ch], MAX_NOISE_NUM);
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
            pp4->noise = avg(noise[ch], MAX_NOISE_NUM);
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

        /*
        if (b)
            bcount++;

        if (bcount > 32) {
            bcount = 0;
            b = false;

            std::cout << " Channel: " << ch << std::endl;
            std::cout << " Data: " << uint32ToSring(corr[ch]) << std::endl;
            std::cout << "~Data: " << uint32ToSring(inv_corr[ch]) << std::endl;
            std::cout << " Expt: " << uint32ToSring(prea) << std::endl;
            std::cout << " Err: " << err1 << " " << err2 << " " << err3 << " " << err4 << " " << std::endl;
        }
        */

    }

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

float OpenUNBDemodulator::avg(float* data, size_t size) {
    float ret = 0;
    for (size_t i=0; i<size; i++) {
        ret += data[i];
    }
    return ret/size;
}

void OpenUNBDemodulator::dumpToFile(std::string fileName, void *data, size_t size) {
    std::fstream file(fileName);
    if (file.is_open()) {
        file.write((char*)data, size);
        file.close();
    }
}
