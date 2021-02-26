#ifndef DECODER_H
#define DECODER_H

#include <thread>
#include <mutex>
#include <unistd.h>
#include <queue>
#include <iostream>
#include <fstream>
#include <cstring>

#include <OpenUNBCodecHL.h>

#include "preamblepoint.h"
#include "crc24.h"

class OpenUNBDecoder
{
public:
    OpenUNBDecoder(int symLen);

    void pushPreablePoint(PreamblePoint* pp);
    int size();
    void setCallback(void (*clb_f)(uint8_t* data, size_t size));
private:
    std::mutex mut;
    std::queue<PreamblePoint*> preablePointVector;
    std::thread* th;
    int symlen;
    std::vector<uint8_t> bits_to_byties(const std::vector<uint8_t>& in);

    void (*clb_f)(uint8_t* data, size_t size) = nullptr;

    bool isRun;

    void run();

};

#endif // DECODER_H
