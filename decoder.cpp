#include "decoder.h"

Decoder::Decoder() {
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
                    data.push_back(-corr);
                } else {
                    data.push_back(corr);
                }

            }
            std::vector<float> tmp(128);
            for (int i=0; i<128; i+=8) {
                for (int j=0; j<8; j++) {
                    tmp[i + j] = data[127 - (i+(7-j)) + 32];
                }
            }
            for (int i=0; i<128; i+=1) {
                std::cout << (tmp[i] > 0 ? "0":"1");
            }
            std::cout << std::endl;

            delete pp;

        } else {
            usleep(100);
        }
    }
}
