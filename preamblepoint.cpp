#include "preamblepoint.h"

PreamblePoint::PreamblePoint()
{

}

PreamblePoint::PreamblePoint(PreamblePoint* pp) {
    preableErrors = pp->preableErrors;
    channel = pp->channel;
    inv = pp->inv;
    noise = pp->noise;
    pos = pp->pos;
    batch = pp->batch;
}
