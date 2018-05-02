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

#if defined (__USE_COMPLEX__)
/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexDoubleInt{
  parMatrixSparse<std::complex<double>,int> parMatrix;
};
struct parMatrixSparseComplexDoubleInt *newParMatrixSparseComplexDoubleInt(void){
  return new struct parMatrixSparseComplexDoubleInt;
}
void ReleaseParMatrixSparseComplexDoubleInt(struct parMatrixSparseComplexDoubleInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexDoubleInt *m, int *size,int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<std::complex<double> >::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }

  *size2 = count2;
/*
  *rows = (int *)malloc(count*sizeof(int));
//  *rows = new int [count];
  for(int i = 0; i < count; i++){
//    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
    *(*rows+i) = 1;

  }

  *size = count;

  std::vector<std::complex<double> >::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }


  *cols = (int *)malloc(count2);
  *real = (double *)malloc(count2);
  *imag = (double *)malloc(count2);

  for(int i = 0; i < count2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
*/
}

void Loc_CSRGetRowsArray2(struct parMatrixSparseComplexDoubleInt *m, int size, int **rows){
  std::vector<int>::iterator it;
  for(int i = 0; i < size; i++){
//    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
    *(*rows+i) = 1;

  }
}

void Loc_CSRGetColsArray(struct parMatrixSparseComplexDoubleInt *m, int **cols, int *size){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }

  *cols = (int *)malloc(count);
  for(int i = 0; i < count; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
  }
  *size = count;
}

void smg2sComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<double>,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#else
/*parMatrixSparse double int C wrapper*/
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

void LOC_MatViewDoubleInt(struct parMatrixSparseDoubleInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeDoubleInt(struct parMatrixSparseDoubleInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRDoubleInt(struct parMatrixSparseDoubleInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}

void smg2sDoubleInt(struct parMatrixSparseDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<double,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}
#endif

#ifdef __cplusplus
};
#endif
