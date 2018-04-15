#ifndef __SMG2S_H__
#define __SMG2S_H__

//#include "../parVector/parVector.cc"
//#include "../../parMatrix/parMatrixSparse.cc"
#include "specGen.h"
#include <math.h>
#include <complex.h>

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

template<typename T, typename S>
parMatrixSparse<T,S> *smg2s(S probSize, Nilpotency<S> nilp, S lbandwidth);

#endif