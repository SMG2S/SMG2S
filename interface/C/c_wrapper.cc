/*

MIT License

Copyright (c) 2019 Xinzhe WU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
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

struct NilpotencyLongInt{
  Nilpotency<__int64_t> nilp;
};

struct NilpotencyInt *newNilpotencyInt(void){
  return new struct NilpotencyInt;
}

struct NilpotencyLongInt *newNilpotencyLongInt(void){
  return new struct NilpotencyLongInt;
}

void ReleaseNilpotencyInt(struct NilpotencyInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void ReleaseNilpotencyLongInt(struct NilpotencyLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void NilpType1(struct NilpotencyInt *n, int num, int size){
  n->nilp.NilpType1(num,size);
}

void NilpType2(struct NilpotencyInt *n, int num, int size){
  n->nilp.NilpType2(num,size);
}

void NilpType1Long(struct NilpotencyLongInt *n, __int64_t num, __int64_t size){
  n->nilp.NilpType1(num,size);
}

void NilpType2Long(struct NilpotencyLongInt *n, __int64_t num, __int64_t size){
  n->nilp.NilpType2(num,size);
}

void showNilpotencyInt(struct NilpotencyInt *n){
  printf("Info ]> Show Nilptent Matrix Informations\n");
  printf("        diagPosition = %d\n",n->nilp.diagPosition);
  printf("        nbOne = %d\n",n->nilp.nbOne);
  printf("        matrix_size = %d\n",n->nilp.matrix_size);
  printf("        nilpotency = %d\n",n->nilp.nilpotency);
}

void showNilpotencyLongInt(struct NilpotencyLongInt *n){
  printf("Info ]> Show Nilptent Matrix Informations\n");
  printf("        diagPosition = %lli\n",n->nilp.diagPosition);
  printf("        nbOne = %lli\n",n->nilp.nbOne);
  printf("        matrix_size = %lli\n",n->nilp.matrix_size);
  printf("        nilpotency = %lli\n",n->nilp.nilpotency);
}

/*parVectorMap C wrapper*/
struct parVectorMapInt{
  parVectorMap<int> parMatrixMap;
};

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

void GetLocalSizeComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, __int64_t *rs, __int64_t *cs){
  __int64_t a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizesComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, __int64_t *size,__int64_t *size2){
  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<std::complex<double> >::iterator itm;
  __int64_t count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }

  *size2 = count2;
}

void Loc_CSRGetRowsArraysComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, __int64_t size, int **rows, __int64_t size2, int **cols, double **real, double **imag){
  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sComplexDoubleLongInt(struct parMatrixSparseComplexDoubleLongInt *m, __int64_t probSize, struct NilpotencyLongInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<double>,__int64_t>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

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


void Loc_CSRGetRowsArraySizesComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int *size,int *size2){
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

void Loc_CSRGetRowsArraysComplexDoubleInt(struct parMatrixSparseComplexDoubleInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag){
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

//complex  single long int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexSingleLongInt{
  parMatrixSparse<std::complex<float>,__int64_t> parMatrix;
};
struct parMatrixSparseComplexSingleLongInt *newParMatrixSparseComplexSingleLongInt(void){
  return new struct parMatrixSparseComplexSingleLongInt;
}
void ReleaseParMatrixSparseComplexSingleLongInt(struct parMatrixSparseComplexSingleLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewComplexSingleLongInt(struct parMatrixSparseComplexSingleLongInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeComplexSingleLongInt(struct parMatrixSparseComplexSingleLongInt *m, __int64_t *rs, __int64_t *cs){
  __int64_t a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRComplexSingleLongInt(struct parMatrixSparseComplexSingleLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizesComplexSingleLongInt(struct parMatrixSparseComplexSingleLongInt *m, __int64_t *size,__int64_t *size2){
  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<std::complex<float> >::iterator itm;
  __int64_t count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }

  *size2 = count2;
}

void Loc_CSRGetRowsArraysComplexSingleLongInt(struct parMatrixSparseComplexSingleLongInt *m, __int64_t size, int **rows, __int64_t size2, int **cols, double **real, double **imag){
  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sComplexSingleLongInt(struct parMatrixSparseComplexSingleLongInt *m, __int64_t probSize, struct NilpotencyLongInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<float>,__int64_t>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

//real double long int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseRealDoubleLongInt{
  parMatrixSparse<double,__int64_t> parMatrix;
};
struct parMatrixSparseRealDoubleLongInt *newParMatrixSparseRealDoubleLongInt(void){
  return new struct parMatrixSparseRealDoubleLongInt;
}
void ReleaseParMatrixSparseRealDoubleLongInt(struct parMatrixSparseRealDoubleLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewRealDoubleLongInt(struct parMatrixSparseRealDoubleLongInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeRealDoubleLongInt(struct parMatrixSparseRealDoubleLongInt *m, __int64_t *rs, __int64_t *cs){
  __int64_t a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRRealDoubleLongInt(struct parMatrixSparseRealDoubleLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizesRealDoubleLongInt(struct parMatrixSparseRealDoubleLongInt *m, __int64_t *size,__int64_t *size2){
  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<double>::iterator itm;
  __int64_t count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(*itm!= 0){
      count2++;
    }
  }

  *size2 = count2;
}

void Loc_CSRGetRowsArraysRealDoubleLongInt(struct parMatrixSparseRealDoubleLongInt *m, __int64_t size, int **rows, __int64_t size2, int **cols, double **vals){
  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*vals+i) = m->parMatrix.CSR_loc->vals[i];
  }
}

void smg2sRealDoubleLongInt(struct parMatrixSparseRealDoubleLongInt *m, __int64_t probSize, struct NilpotencyLongInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<double,__int64_t>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

//complex single int

/*parMatrixSparse complex<double> int C wrapper*/
struct parMatrixSparseComplexSingleInt{
  parMatrixSparse<std::complex<float>,int> parMatrix;
};
struct parMatrixSparseComplexSingleInt *newParMatrixSparseComplexSingleInt(void){
  return new struct parMatrixSparseComplexSingleInt;
}
void ReleaseParMatrixSparseComplexSingleInt(struct parMatrixSparseComplexSingleInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}
extern void LOC_MatViewComplexInt(struct parMatrixSparseComplexSingleInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeComplexInt(struct parMatrixSparseComplexSingleInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRComplexInt(struct parMatrixSparseComplexSingleInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}


void Loc_CSRGetRowsArraySizesComplexSingleInt(struct parMatrixSparseComplexSingleInt *m, int *size,int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<std::complex<float>>::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(itm->imag() !=0 || itm->real() != 0){
      count2++;
    }
  }

  *size2 = count2;
}

void Loc_CSRGetRowsArraysComplexSingleInt(struct parMatrixSparseComplexSingleInt *m, int size, int **rows, int size2, int **cols, double **real, double **imag){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*real+i) = m->parMatrix.CSR_loc->vals[i].real();
    *(*imag+i) = m->parMatrix.CSR_loc->vals[i].imag();
  }
}

void smg2sComplexSingleInt(struct parMatrixSparseComplexSingleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<std::complex<float>,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

//real double int
/*parMatrixSparse double int C wrapper*/
struct parMatrixSparseRealDoubleInt{
  parMatrixSparse<double,int> parMatrix;
};

struct parMatrixSparseRealDoubleInt *newParMatrixSparseRealDoubleInt(void){
  return new struct parMatrixSparseRealDoubleInt;
}

void ReleaseParMatrixSparseRealDoubleInt(struct parMatrixSparseRealDoubleInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void LOC_MatViewRealDoubleInt(struct parMatrixSparseRealDoubleInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeRealDoubleInt(struct parMatrixSparseRealDoubleInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRRealDoubleInt(struct parMatrixSparseRealDoubleInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}

void Loc_RealCSRGetRowsArraySizesRealDoubleInt(struct parMatrixSparseRealDoubleInt *m, int *size, int *size2){
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
void Loc_RealCSRGetRowsArraysRealDoubleInt(struct parMatrixSparseRealDoubleInt *m, int size, int **rows, int size2, int **cols, double **vals){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*vals+i) = m->parMatrix.CSR_loc->vals[i];
  }
}


void smg2sRealDoubleInt(struct parMatrixSparseRealDoubleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<double,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

//real single long int


/*parMatrixSparse long int C wrapper*/
struct parMatrixSparseRealSingleLongInt{
  parMatrixSparse<float,__int64_t> parMatrix;
};

struct parMatrixSparseRealSingleLongInt *newParMatrixSparseRealSingleLongInt(void){
  return new struct parMatrixSparseRealSingleLongInt;
}

void ReleaseParMatrixSparseRealSingleLongInt(struct parMatrixSparseRealSingleLongInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void LOC_MatViewRealSingleLongInt(struct parMatrixSparseRealSingleLongInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeRealSingleLongInt(struct parMatrixSparseRealSingleLongInt *m, __int64_t *rs, __int64_t *cs){
  __int64_t a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRRealSingleLongInt(struct parMatrixSparseRealSingleLongInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}

void Loc_LongCSRGetRowsArraySizesRealSingleLongInt(struct parMatrixSparseRealSingleLongInt *m, __int64_t *size, __int64_t *size2){
  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<float>::iterator itm;
  __int64_t count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(*itm!= 0){
      count2++;
    }
  }

  *size2 = count2;
}
void Loc_RealCSRGetRowsArraysRealSingleLongInt(struct parMatrixSparseRealSingleLongInt *m, __int64_t size, __int64_t **rows, __int64_t size2, __int64_t **cols, double **vals){
  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*vals+i) = m->parMatrix.CSR_loc->vals[i];
  }
}


void smg2sRealSingleLongInt(struct parMatrixSparseRealSingleLongInt *m, __int64_t probSize, struct NilpotencyLongInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<float,__int64_t>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

/*parMatrixSparse int C wrapper*/
struct parMatrixSparseRealSingleInt{
  parMatrixSparse<float,int> parMatrix;
};

struct parMatrixSparseRealSingleInt *newParMatrixSparseRealSingleInt(void){
  return new struct parMatrixSparseRealSingleInt;
}

void ReleaseParMatrixSparseRealSingleInt(struct parMatrixSparseRealSingleInt **ppInstance){
  delete *ppInstance;
  *ppInstance = 0;
}

void LOC_MatViewRealSingleInt(struct parMatrixSparseRealSingleInt *m){
  m->parMatrix.LOC_MatView();
}

void GetLocalSizeRealSingleInt(struct parMatrixSparseRealSingleInt *m, int *rs, int *cs){
  int a, b;
  m->parMatrix.GetLocalSize(a, b);
  *rs =a; *cs = b;
}

void Loc_ConvertToCSRRealSingleInt(struct parMatrixSparseRealSingleInt *m){
  m->parMatrix.Loc_ConvertToCSR();
}

void Loc_RealCSRGetRowsArraySizesRealSingleInt(struct parMatrixSparseRealSingleInt *m, int *size, int *size2){
  std::vector<int>::iterator it;
  int count = 0;
  for(it = m->parMatrix.CSR_loc->rows.begin(); it != m->parMatrix.CSR_loc->rows.end(); ++it){
      count++;
  }
  *size = count;

  std::vector<float>::iterator itm;
  int count2 = 0;
  for(itm =  m->parMatrix.CSR_loc->vals.begin(); itm !=  m->parMatrix.CSR_loc->vals.end(); ++itm){
    if(*itm!= 0){
      count2++;
    }
  }

  *size2 = count2;
}
void Loc_RealCSRGetRowsArraysRealSingleInt(struct parMatrixSparseRealSingleInt *m, int size, int **rows, int size2, int **cols, double **vals){
  for(int i = 0; i < size; i++){
    *(*rows+i) = m->parMatrix.CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = m->parMatrix.CSR_loc->cols[i];
    *(*vals+i) = m->parMatrix.CSR_loc->vals[i];
  }
}


void smg2sRealSingleInt(struct parMatrixSparseRealSingleInt *m, int probSize, struct NilpotencyInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  m->parMatrix = *smg2s<float,int>(probSize, nilp->nilp, lbandwidth,spectrum,comm);
}

#ifdef __cplusplus
};
#endif
