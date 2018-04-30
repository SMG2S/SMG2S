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

#include "c_wrapper.h"
#include "../../utils/utils.h"
#include "../../parMatrix/parMatrixSparse.h"
#include "../../smg2s/smg2s.h"

#ifdef __cplusplus
extern "C" {
#endif

/*Nilpotency Matrix C Wrapper*/
struct NilpotencyInt{
  Nilpotency<int> nilp;
};

struct NilpotencyInt *newNilpotencyInt(void){
  return new struct NilpotencyInt;
}

void ReleaseNilpotencyInt(struct NilpotencyInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void NilpType1(struct NilpotencyInt *n, int num, int size){
  n->nilp.NilpType1(num,size);
}

void NilpType2(struct NilpotencyInt *n, int num, int size){
  n->nilp.NilpType2(num,size);
}

void showNilpotencyInt(struct NilpotencyInt *n){
  printf("Info ]> Show Nilptent Matrix Informations\n");
  printf("        diagPosition = %d\n",n->nilp.diagPosition);
  printf("        nbOne = %d\n",n->nilp.nbOne);
  printf("        matrix_size = %d\n",n->nilp.matrix_size);
  printf("        nilpotency = %d\n",n->nilp.nilpotency);
}

/*parVectorMap C wrapper*/
struct parVectorMapInt{
  parVectorMap<int> parMatrixMap;
};

/*parMatrixSparse C wrapper*/
struct parMatrixSparseDoubleInt{
  parMatrixSparse<double,int> parMatrix;
};

struct parMatrixSparseDoubleInt *newParMatrixSparseDoubleInt(void){
  return new struct parMatrixSparseDoubleInt;
}

void ReleaseParMatrixSparseDoubleInt(struct parMatrixSparseDoubleInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void LOC_MatView(struct parMatrixSparseDoubleInt *m){
  m->parMatrix.LOC_MatView();
}

void smg2s(struct parMatrixSparseDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum){
  m->parMatrix = *smg2s<double,int>(probSize, nilp->nilp, lbandwidth,spectrum);
}
#ifdef __cplusplus
};
#endif
