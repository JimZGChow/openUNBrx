#include "decoder.h"

//#define DEC_LOG
//#define FILE_DUMP
#define TEST_COUNTER
#define OUTPUT_MSG

OpenUNBDecoder::OpenUNBDecoder(int symLen) {
    //initOpenUNBCodec();

    isRun = true;

    symlen = symLen;

    initOpenUNBCodec();

    th = new std::thread(&OpenUNBDecoder::run, this);

}

void OpenUNBDecoder::setCallback(void (*_clb_f)(uint8_t* data, size_t size)) {
    clb_f = _clb_f;
}

void OpenUNBDecoder::pushPreablePoint(PreamblePoint* pp) {
    mut.lock();
    preablePointVector.push(pp);
    mut.unlock();
}

int getErrors(std::vector<uint8_t> a, std::vector<uint8_t> b) {
    int err = 0;
    for (int i=0; i< std::min(a.size(), b.size()); i++) {

        if (a[i] != b[i])
            err++;
    }
    return  err;
}


void OpenUNBDecoder::run() {
    std::time_t result = std::time(nullptr);
    while (isRun) {
        if (preablePointVector.size() > 0) {
            mut.lock();
            PreamblePoint* pp = preablePointVector.front();
            preablePointVector.pop();
            mut.unlock();

#ifdef FILE_DUMP
            std::localtime(&result);
            std::string fname = "data/data" + std::to_string(result) + "_" + std::to_string(pp->channel) + ".complex";
            std::fstream file1(fname, std::ios_base::out | std::ios_base::binary);

            if (file1.is_open()) {
                file1.write((char*)pp->data.data(), pp->data.size() * sizeof(float) * 2);
                file1.close();
            }
            else {
                std::cout << "File error: " << std::strerror(errno) << std::endl;
            }
#endif

            std::vector<float> data;
            float max = 0.0;
            for (int i=33; i< pp->data.size(); i += symlen) {
                float corr;
                if (!(pp->inv & PreamblePoint::invIQ)) {
                    corr = (pp->data[i] * std::complex<float>(pp->data[i-1].real(), -pp->data[i-1].imag())).imag();
                } else {
                    corr = (pp->data[i] * std::complex<float>(pp->data[i-1].real(), -pp->data[i-1].imag())).real();
                }

                if (corr > max)
                    max = corr;

                if (pp->inv & PreamblePoint::invBits && !(pp->inv & PreamblePoint::invIQ)) {
                    data.push_back(-corr);
                } else {
                    data.push_back(corr);
                }

            }

            float qAvg = 0;
            float iAvg = 0;
            float rssi = 0;
            for (int i=32; i < 64; i++) {
                rssi += pp->data[i].real() * pp->data[i].real() + pp->data[i].imag() * pp->data[i].imag();
            }
            rssi /= pp->data.size();
            rssi = 10*log10(rssi);

#ifdef DEC_LOG
            std::cout << std::endl << "PP: " <<(int)pp->inv << "    Ch: " << pp->channel << std::endl;

            for (int i=0; i<data.size(); i++) {
                std::cout << (int)(data[i] > 0 ? 1 : 0);
            }
            std::cout << std::endl;

            std::cout << "x: " << pp->pos << std::endl;
            std::cout << "batch: " << pp->batch << std::endl;

#endif

            auto res = decode_64(data);

            if (res.size() > 0) {
                res.erase(res.end() - 10, res.end());
                std::cout << " OpenUNB data! SNR: " << rssi - 10*log10(pp->noise) << ", RSSI: " << rssi << ", NOISE: " << 10*log10(pp->noise) << ", ch: " << pp->channel << std::endl;
                std::cout << " Data: { " << getStringHexFromVector(res) << " }" << std::endl;
            }

            for (int i=0; i<64; i++) {
                //data.push_back(0);
            }
            res = decode_96(data);

            if (res.size() > 0) {
                res.erase(res.end() - 10, res.end());
                std::cout << " OpenUNB data! SNR: " << rssi - 10*log10(pp->noise) << ", RSSI: " << rssi << ", NOISE: " << 10*log10(pp->noise) << ", ch: " << pp->channel << ", iqAvgDiv: " << 20 * log10(iAvg/qAvg) << std::endl;
                std::cout << " Data: { " << getStringHexFromVector(res) << " }" << std::endl;

            }
            delete pp;

        } else {
            usleep(100);
        }
    }
}

std::vector<uint8_t> OpenUNBDecoder::bits_to_byties(const std::vector<uint8_t>& in) {
    std::vector<uint8_t> res(in.size()/8);

    for (int i=0; i < res.size(); i++) {
        res[i] = 0;
    }

    for (int i=0; i < in.size(); i++) {
        res[i / 8] |= in[i] << (i % 8);
    }

    return res;
}

int OpenUNBDecoder::size() {
    return preablePointVector.size();
}
