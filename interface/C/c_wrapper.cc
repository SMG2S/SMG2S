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

#if defined(__USE_COMPLEX__) && defined(__USE_DOUBLE__) && defined (__USE_64BIT__)
//complex double long int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexDoubleLongInt{
  parMatrixSparse<std::complex<double>,__int64_t> parMatrix;
};
struct parMatrixSparseComplexDoubleLongInt *newParMatrixSparseComplexDoubleLongInt(void){
  return new struct parMatrixSparseComplexDoubleLongInt;
}
void ReleaseParMatrixSparseComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexDoubleLongInt *m, int *size,int *size2){
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
}

void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexDoubleLongInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<double>,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#elif defined (__USE_COMPLEX__) && defined(__USE_DOUBLE__)
//complex double int

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
}

void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexDoubleInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<double>,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#elif defined (__USE_COMPLEX__) && defined(__USE_64BIT__)
//complex  single long int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexLongInt{
  parMatrixSparse<std::complex<int>,__int64_t> parMatrix;
};
struct parMatrixSparseComplexLongInt *newParMatrixSparseComplexLongInt(void){
  return new struct parMatrixSparseComplexLongInt;
}
void ReleaseParMatrixSparseComplexLongInt(struct parMatrixSparseComplexLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewComplexLongInt(struct parMatrixSparseComplexLongInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeComplexLongInt(struct parMatrixSparseComplexLongInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRComplexLongInt(struct parMatrixSparseComplexLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexLongInt *m, int *size,int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<std::complex<int> >::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }

  *size2 = count2;
}

void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexLongInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sComplexLongInt(struct parMatrixSparseComplexLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<int>,__int64_t>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#elif defined (__USE_DOUBLE__) && defined(__USE_64BIT__)
//real double long int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseDoubleLongInt{
  parMatrixSparse<double,__int64_t> parMatrix;
};
struct parMatrixSparseDoubleLongInt *newParMatrixSparseDoubleLongInt(void){
  return new struct parMatrixSparseDoubleLongInt;
}
void ReleaseParMatrixSparseDoubleLongInt(struct parMatrixSparseDoubleLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewDoubleLongInt(struct parMatrixSparseDoubleLongInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeDoubleLongInt(struct parMatrixSparseDoubleLongInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRDoubleLongInt(struct parMatrixSparseDoubleLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizes(struct parMatrixSparseDoubleLongInt *m, int *size,int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<double>::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }

  *size2 = count2;
}

void Loc_CSRGetRowsArrays(struct parMatrixSparseDoubleLongInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sDoubleLongInt(struct parMatrixSparseDoubleLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<double,__int64_t>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#elif defined (__USE_COMPLEX__)
//complex single int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexInt{
  parMatrixSparse<std::complex<int>,int> parMatrix;
};
struct parMatrixSparseComplexInt *newParMatrixSparseComplexInt(void){
  return new struct parMatrixSparseComplexInt;
}
void ReleaseParMatrixSparseComplexInt(struct parMatrixSparseComplexInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewComplexInt(struct parMatrixSparseComplexInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeComplexInt(struct parMatrixSparseComplexInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRComplexInt(struct parMatrixSparseComplexInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizes(struct parMatrixSparseComplexInt *m, int *size,int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<std::complex<int>>::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }

  *size2 = count2;
}

void Loc_CSRGetRowsArrays(struct parMatrixSparseComplexInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sComplexInt(struct parMatrixSparseComplexInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<int>,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#elif defined (__USE_DOUBLE__)
//real double int
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

void Loc_RealCSRGetRowsArraySizes(struct parMatrixSparseDoubleInt *m, int *size, int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<double>::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(*itm!= 0){
      count2++;
    }
  }

  *size2 = count2;
}
void Loc_RealCSRGetRowsArrays(struct parMatrixSparseDoubleInt *m, int size, int **rows, int size2, int **cols, double **vals){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*vals+i) = m->parMatrix.CSR_loc->vals[i];
  }
}


void smg2sDoubleInt(struct parMatrixSparseDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<double,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#elif defined (__USE_64BIT__)
//real single long int


/*parMatrixSparse long int C wrapper*/
struct parMatrixSparseLongInt{
  parMatrixSparse<long,int> parMatrix;
};

struct parMatrixSparseLongInt *newParMatrixSparseLongInt(void){
  return new struct parMatrixSparseDoubleInt;
}

void ReleaseParMatrixSparseLongInt(struct parMatrixSparseLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void LOC_MatViewLongInt(struct parMatrixSparseLongInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeLongInt(struct parMatrixSparseLongInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRLongInt(struct parMatrixSparseLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}

void Loc_LongCSRGetRowsArraySizes(struct parMatrixSparseLongInt *m, int *size, int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<long>::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(*itm!= 0){
      count2++;
    }
  }

  *size2 = count2;
}
void Loc_RealCSRGetRowsArrays(struct parMatrixSparseLongInt *m, int size, int **rows, int size2, int **cols, double **vals){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*vals+i) = m->parMatrix.CSR_loc->vals[i];
  }
}


void smg2sLongInt(struct parMatrixSparseLongInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<long,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#else

/*parMatrixSparse int C wrapper*/
struct parMatrixSparseInt{
  parMatrixSparse<int,int> parMatrix;
};

struct parMatrixSparseInt *newParMatrixSparseInt(void){
  return new struct parMatrixSparseInt;
}

void ReleaseParMatrixSparseInt(struct parMatrixSparseInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void LOC_MatViewInt(struct parMatrixSparseInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeInt(struct parMatrixSparseInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRInt(struct parMatrixSparseInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}

void Loc_RealCSRGetRowsArraySizes(struct parMatrixSparseInt *m, int *size, int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<int>::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(*itm!= 0){
      count2++;
    }
  }

  *size2 = count2;
}
void Loc_RealCSRGetRowsArrays(struct parMatrixSparseInt *m, int size, int **rows, int size2, int **cols, double **vals){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*vals+i) = m->parMatrix.CSR_loc->vals[i];
  }
}


void smg2sInt(struct parMatrixSparseInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<int,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}
#endif

#ifdef __cplusplus
};
#endif
