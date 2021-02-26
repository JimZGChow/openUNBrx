#include "decoder.h"

OpenUNBDecoder::OpenUNBDecoder(int symLen) {
    initOpenUNBCodec();

    isRun = true;

    symlen = symLen;

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

//#define DEC_LOG
//#define FILE_DUMP

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
            for (int i=1; i< pp->data.size(); i += symlen) {
                float corr;
                if (!(pp->inv & PreamblePoint::invIQ)) {
                    //corr = (pp->data[i - 1].real() * pp->data[i].real() + pp->data[i - 1].imag() * pp->data[i].imag()) > 0 ? 1:0;
                    //corr = pp->data[i - 1 * symlen].real() * pp->data[i].real() + pp->data[i - 1 * symlen].imag() * pp->data[i].imag();
                    corr = (pp->data[i] * std::complex<float>(pp->data[i-1].real(), -pp->data[i-1].imag())).imag();
                } else {
                    //corr = (-pp->data[i - 1].imag() * pp->data[i].real() + pp->data[i - 1].real() * pp->data[i].imag()) > 0 ? 1:0;
                    //corr = -pp->data[i - 1 * symlen].imag() * pp->data[i].real() + pp->data[i - 1 * symlen].real() * pp->data[i].imag();
                    corr = (pp->data[i] * std::complex<float>(pp->data[i-1].real(), -pp->data[i-1].imag())).real();
                }

                if (corr > max)
                    max = corr;

                if (pp->inv & PreamblePoint::invBits) {
                    data.push_back(corr);
                } else {
                    data.push_back(corr);
                }

            }
            std::vector<float> tmp(256);
            std::vector<uint8_t> tmpb(256);
            //std::vector<uint8_t> tmpexp_enc = getVectorFromStringHex("6F7A15976F7A15976F7A15976F7A15976F7A15976F7A15976F7A15976F7A1597");
            //std::vector<uint8_t> tmpexp_enc = getVectorFromStringHex("000000000000000007EA4AFD97B6DAA1234E0E39BD75ACC1A4F60755EED58E4A");
            //std::vector<uint8_t> tmpexp = getVectorFromStringHex("4BF2EC8D3838EF5AF49EC177");//{0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0};

            for (int i=0; i<256; i++) {
                tmp[i] = data[i + 32];
                tmpb[i] = data[i + 32] > 0 ? 0:1;
            }

            float qAvg = 0;
            float iAvg = 0;
            float rssi = 0;
            for (int i=32; i < 64; i++) {
                qAvg += fabs(pp->data[i].imag());
                iAvg += fabs(pp->data[i].real());
                rssi += pp->data[i].real() * pp->data[i].real() + pp->data[i].imag() * pp->data[i].imag();
            }
            rssi /= pp->data.size();
            rssi = 10*log10f(rssi);

            //std::cout << pp->channel << ") Errors: " << getErrors(tmpb, tmpexp_enc) << std::endl;
#ifdef DEC_LOG
            if (pp->channel == 20 || true) {
            std::cout << std::endl << "PP: " <<(int)pp->inv << "    Ch: " << pp->channel << std::endl;
            std::cout << "Recieved data: \n" << getStringBinFromVector(tmpb) << std::endl;
            //std::cout << "Expected data: \n" << getStringBinFromVector(tmpexp_enc) << std::endl;
            //std::cout << "Exp == Rec: " << (int) (tmpb == tmpexp_enc) << std::endl;
            //std::cout << "Errors: " << getErrors(tmpb, tmpexp_enc) << std::endl;
            //std::cout << "iqAvgDiv: " << std::max(qAvg, iAvg) / std::min(qAvg, iAvg) << std::endl;
            std::cout << "x: " << pp->pos << std::endl;
            std::cout << "batch: " << pp->batch << std::endl;
            }
#endif

            auto res = decode_64(tmp);

            if (res.size() > 0) {
                std::cout << " OpenUNB data! SNR: " << rssi - 10*log10f(pp->noise) << ", RSSI: " << rssi << ", NOISE: " << 10*log10f(pp->noise) << ", ch: " << pp->channel << ", iqAvgDiv: " << 20 * log10f(iAvg/qAvg) << std::endl;
                //std::cout << " Data: { " << getStringHexFromVector(res) << " }" << std::endl;

                std::vector<uint8_t> dataHex = bits_to_byties(res);
                uint32_t addr = 0;
                uint16_t payload = 0;
                uint32_t crcIn = 0;
                memcpy(&addr, dataHex.data(), 3);
                memcpy(&payload, dataHex.data() + 3, 2);
                memcpy(&crcIn, dataHex.data() + 5, 3);

                char txOut[128];
                sprintf(txOut, "Data: { addr: %.6X,\t payload: %.4X,\t crc: %.6X}", addr, payload, crcIn);
                std::cout << "     " << txOut << std::endl;

                uint32_t crcCalc = crc24(dataHex.data(), 5);

                if (crcCalc == crcIn) {
                    std::cout << "         CRC OK" <<std::endl;

                    std::localtime(&result);
                    std::string fname = "msg.txt";
                    std::fstream file1(fname, std::ios_base::out | std::ios_base::app);

                    if (file1.is_open()) {
                        file1 << std::asctime(std::localtime(&result)) << " \t" << txOut << std::endl;
                        file1.close();
                    }
                    else {
                        std::cout << "File error: " << std::strerror(errno) << std::endl;
                    }
                }
            }

//#endif
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
