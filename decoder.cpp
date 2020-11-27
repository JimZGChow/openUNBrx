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

void Decoder::run() {
    while (isRun) {
        if (preablePointVector.size() > 0) {
            mut.lock();
            PreamblePoint* pp = preablePointVector.front();
            preablePointVector.pop();
            mut.unlock();

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
            std::vector<float> tmp(128);
            std::vector<uint8_t> tmpb(128);
            std::vector<uint8_t> tmpexp = {0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0};
            for (int i=0; i<128; i+=8) {
                for (int j=0; j<8; j++) {
                    tmp[i + j] = data[127 - (i+(7-j)) + 32];
                    tmpb[i + j] = data[127 - (i+(7-j)) + 32] > 0 ? 0:1;
                }
            }
            std::cout << "PP: " <<(int)pp->inv << std::endl;
            std::cout << "Recieved data: " << getStringBinFromVector(tmpb) << std::endl;
            std::cout << "Expected data: " << getStringBinFromVector(tmpexp) << std::endl;
            std::cout << "Errors: " << getErrors(tmpb, tmpexp) << std::endl;

            std::vector<std::vector<float>> prob;
            std::vector<std::vector<uint8_t>> dec = pcscl_prep(7, 16, tmp, &prob, info_bit_pattern_64);
            dec = polar_transform_noperm(dec);
            dec = extract_with_filter(dec, crc_mask_64, num_of_nonzero_bits_64 - short_64);

            std::vector<uint8_t> crc_err;
            crc_err = crc_ok_array(0x327, dec);

            for (int i = 0; i < crc_err.size(); i++) {
                if (crc_err[i]) {
                    std::cout << "Decoded data: " << getStringHexFromVector(dec[i]) << std::endl;
                }
            }

            delete pp;

        } else {
            usleep(100);
        }
    }
}
