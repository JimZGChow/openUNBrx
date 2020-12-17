#include "decoder.h"

uint8_t iwd_full_64[] = { 1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,1,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,1,0,1};
uint8_t iwd_full_96[] = {0,0,0,1,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,0,1,1,1,0,1,1,1,0,0,1,1,1,0,0,1,1,0,0,0,1,1,1,1,0,1,1,0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,1,1,1,0,1,1,1,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,1,0,1,0,1,0,1,0,1,1,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


uint8_t info_bit_pattern_64[] = {0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
const int num_of_nonzero_bits_64 = 74;
const int short_64 = 0;
uint8_t info_bit_pattern_96[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
const int num_of_nonzero_bits_96 = 170;
const int short_96 = 64;

int* short_mask;
uint8_t* frozen_indicator_64;
uint8_t* frozen_indicator_96;
uint16_t* crc_mask_64;
uint16_t* crc_mask_96;

Decoder::Decoder() {
    frozen_indicator_64 = new uint8_t[sizeof(info_bit_pattern_64)];
        frozen_indicator_96 = new uint8_t[sizeof(info_bit_pattern_96)];

        for (int i = 0; i < sizeof(info_bit_pattern_64); i++) {
            frozen_indicator_64[i] = !info_bit_pattern_64[i];
        }

        for (int i = 0; i < sizeof(info_bit_pattern_96); i++) {
            frozen_indicator_96[i] = !info_bit_pattern_96[i];
        }

        short_mask = new int[64];

        for (int i = 0, j = sizeof(info_bit_pattern_96) - 1; i < 64; j--) {
            if (info_bit_pattern_96[j] == 1) {
                short_mask[i++] = j + 1;
            }
        }

        for (int i = 0; i < 32; i++) {
            int tmp = short_mask[i];
            short_mask[i] = short_mask[63 - i];
            short_mask[63 - i] = tmp;
        }

        crc_mask_64 = new uint16_t[num_of_nonzero_bits_64 - short_64];
        crc_mask_96 = new uint16_t[num_of_nonzero_bits_96 - short_96];

        for (int i = 0, j = 0; j < (num_of_nonzero_bits_64 - short_64); i++) {
            if (info_bit_pattern_64[i]) {
                crc_mask_64[j++] = i;
            }
        }

        for (int i = 0, j = 0; j < (num_of_nonzero_bits_96 - short_96); i++) {
            if (info_bit_pattern_96[i]) {
                crc_mask_96[j++] = i;
            }
        }

    isRun = true;

    th = new std::thread(&Decoder::run, this);
}

void Decoder::pushPreablePoint(PreamblePoint* pp) {
    mut.lock();
    preablePointVector.push(pp);
    mut.unlock();
}

//#define DEC_LOG

void Decoder::run() {
    std::time_t result = std::time(nullptr);
    while (isRun) {
        if (preablePointVector.size() > 0) {
            mut.lock();
            PreamblePoint* pp = preablePointVector.front();
            preablePointVector.pop();
            mut.unlock();

            /*
            std::string fname = "data" + std::to_string(std::localtime(&result)->tm_gmtoff) + "_" + std::to_string(pp->channel) + ".complex";
            std::fstream file1(fname, std::ios_base::out | std::ios_base::binary);

            if (file1.is_open()) {
                file1.write((char*)pp->data.data(), pp->data.size() * sizeof(float) * 2);
                file1.close();
            }
            else {
                std::cout << "File error: " << std::strerror(errno) << std::endl;
            }
            */

            std::vector<float> data;
            for (int i=1; i< pp->data.size(); i++) {
                float corr;
                if (!(pp->inv & PreamblePoint::invIQ)) {
                    //corr = (pp->data[i - 1].real() * pp->data[i].real() + pp->data[i - 1].imag() * pp->data[i].imag()) > 0 ? 1:0;
                    corr = pp->data[i - 1].real() * pp->data[i].real() + pp->data[i - 1].imag() * pp->data[i].imag();
                } else {
                    //corr = (-pp->data[i - 1].imag() * pp->data[i].real() + pp->data[i - 1].real() * pp->data[i].imag()) > 0 ? 1:0;
                    corr = -pp->data[i - 1].imag() * pp->data[i].real() + pp->data[i - 1].real() * pp->data[i].imag();
                }

                if (pp->inv & PreamblePoint::invBits) {
                    data.push_back(corr);
                } else {
                    data.push_back(corr);
                }

            }
            std::vector<float> tmp(256);
            std::vector<uint8_t> tmpb(256);
            std::vector<uint8_t> tmpexp_enc = getVectorFromStringHex("000000000000000007EA4AFD97B6DAA1234E0E39BD75ACC1A4F60755EED58E4A");
            std::vector<uint8_t> tmpexp = getVectorFromStringHex("4BF2EC8D3838EF5AF49EC177");//{0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0};
            for (int i=0; i<256; i+=8) {
                for (int j=0; j<8; j++) {
                    tmp[i + j] = data[255 - (i+(7-j)) + 32];
                    tmpb[i + j] = data[255 - (i+(7-j)) + 32] > 0 ? 0:1;
                }
            }
            float qAvg = 0;
            float iAvg = 0;
            float rssi = 0;
            for (int i=0; i < pp->data.size(); i++) {
                qAvg += fabs(pp->data[i].imag());
                iAvg += fabs(pp->data[i].real());
                rssi += pp->data[i].real() * pp->data[i].real() + pp->data[i].imag() * pp->data[i].imag();
            }
            rssi /= pp->data.size();
            rssi = 10*log10f(rssi);

#ifdef DEC_LOG
            std::cout << "PP: " <<(int)pp->inv << "    Ch: " << pp->channel << std::endl;
            std::cout << "Recieved data: \n" << getStringBinFromVector(tmpb) << std::endl;
            std::cout << "Expected data: \n" << getStringBinFromVector(tmpexp_enc) << std::endl;
            std::cout << "Errors: " << getErrors(tmpb, tmpexp_enc) << std::endl;
            std::cout << "qAvg: " << qAvg/iAvg << std::endl;
            std::cout << "x: " << pp->pos << std::endl;
#endif
            for (int i=0; i<short_96; i++) {
                tmp[short_mask[i]] = 10000;
            }

            std::vector<std::vector<float>> prob;
            std::vector<std::vector<uint8_t>> dec = pcscl_prep(8, 16, tmp, &prob, info_bit_pattern_96);
            dec = polar_transform_noperm(dec);
            dec = extract_with_filter(dec, crc_mask_96, num_of_nonzero_bits_96 - short_96);

            std::vector<uint8_t> crc_err;
            crc_err = crc_ok_array(0x327, dec);

//#ifdef DEC_LOG
            int i;
            bool b = false;
            float ber;
            float maxProbe = -100000000;
            int maxProbeIndex = -1;
            int err;
            for (int i = 0; i < crc_err.size(); i++) {
                if (crc_err[i] == 1) {
                    if (maxProbe < prob[i][prob.size()-1]) {
                        maxProbe = prob[i][prob.size()-1];
                        maxProbeIndex = i;
                    }

                    /*
                    std::vector<uint8_t> dec_data = remove_crc(dec[i]);
                    std::cout << "Decoded data:  " << getStringHexFromVector(dec_data) << std::endl;
                    std::cout << "Expected data: " << getStringHexFromVector(tmpexp) << std::endl;
                    //std::cout << " SNR: " << rssi - 10*log10f(pp->noise) << ", Er: " << getErrors(tmpb, tmpexp_enc) << ", ber: " << ber << ", ch: " << pp->channel << std::endl;
                    if (dec_data == tmpexp)
                        b = true;
                    //}
                    */
                }
            }

            if (maxProbeIndex != -1) {
                err = getErrors(remove_crc(dec[maxProbeIndex]), tmpexp);
                ber = (float)err / 96.0;
                if (ber == 0.0)
                    b = true;

            } else {
                err = getErrors(remove_crc(dec[0]), tmpexp);
                ber = (float)err / 96.0;
            }

            std::fstream file("data5.data", std::ios_base::in | std::ios_base::app);
            if (file.is_open()) {
                file /* << std::asctime(std::localtime(&result))*/ << "\t\t" << rssi - 10*log10f(pp->noise) << "\t\t" << getStringHexFromVector(remove_crc(dec[i > crc_err.size() ? 0 : i])) << "\t\t" << ber << "\t\t" << ber << "\t\t" << b << std::endl;
                file.close();
            }
            if (b)
                std::cout << " SNR: " << rssi - 10*log10f(pp->noise) << ", Er: " << err << ", ber: " << ber << ", ch: " << pp->channel << ", b: " << b << std::endl;

//#endif
            delete pp;

        } else {
            usleep(100);
        }
    }
}

int Decoder::size() {
    return preablePointVector.size();
}
