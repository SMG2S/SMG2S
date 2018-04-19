/*
   This file is part of SMG2S.
   Author(s): Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr or xinzhe.wu1990@gmail.com>
        Date: 2018-04-20
   Copyright (C) 2018-     Xinzhe WU
   
   SMG2S is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   SMG2S is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with SMG2S.  If not, see <http://www.gnu.org/licenses/>.
*/

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