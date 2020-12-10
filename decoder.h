#ifndef DECODER_H
#define DECODER_H

#include <thread>
#include <mutex>
#include <unistd.h>
#include <queue>
#include <iostream>
#include <fstream>
#include <cstring>

#include "preamblepoint.h"
#include "pcscl_noperm.h"
#include "VectorHelper.h"

class Decoder
{
public:
    Decoder();

    void pushPreablePoint(PreamblePoint* pp);
    int size();
private:
    std::mutex mut;
    std::queue<PreamblePoint*> preablePointVector;
    std::thread* th;

    bool isRun;

    void run();

};

#endif // DECODER_H
