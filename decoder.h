#ifndef DECODER_H
#define DECODER_H

#include <thread>
#include <mutex>
#include <unistd.h>
#include <queue>
#include <iostream>
#include "preamblepoint.h"

class Decoder
{
public:
    Decoder();

    void pushPreablePoint(PreamblePoint* pp);
private:
    std::mutex mut;
    std::queue<PreamblePoint*> preablePointVector;
    std::thread* th;

    bool isRun;

    void run();

};

#endif // DECODER_H
