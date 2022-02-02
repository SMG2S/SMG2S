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
#include "../../parMatrix/parMatrixSparse.h"
#include "../../smg2s/smg2s.h"
#include "../../smg2s/smg2s_nonsymmetric.h"

struct NilpInt{
  void *obj;
};

NilpInt_t *newNilpInt(){
  
  NilpInt_t *nilp;
  Nilpotency<int> *obj;

  nilp = (NilpInt *)malloc(sizeof(*nilp));
  obj = new Nilpotency<int>();
  nilp->obj = obj;

  return nilp;

}

void nilpInt_destory(NilpInt_t *nilp){
  if(nilp == NULL)
    return;

  delete nilp;

}

void nilpIntType1(NilpInt_t *nilp, int num, int size){
  Nilpotency<int> *obj;
  if(nilp == NULL)
    return;

  obj = static_cast<Nilpotency<int> *>(nilp->obj);
  obj->NilpType1(num, size);
}

void nilpIntShow(struct NilpInt *nilp){

  Nilpotency<int> *obj;
  if(nilp == NULL)
    return;

  obj = static_cast<Nilpotency<int> *>(nilp->obj);

  printf("Info ]> Show Nilptent Matrix Informations\n");
  printf("        diagPosition = %d\n",obj->diagPosition);
  printf("        nbOne = %d\n",obj->nbOne);
  printf("        matrix_size = %d\n",obj->matrix_size);
  printf("        nilpotency = %d\n",obj->nilpotency);
}


////
struct NilpLong{
  void *obj;
};

NilpLong_t *newNilpLong(){
  
  NilpLong_t *nilp;
  Nilpotency<long> *obj;

  nilp = (NilpLong_t *)malloc(sizeof(*nilp));
  obj = new Nilpotency<long>();
  nilp->obj = obj;

  return nilp;

}

void nilpLong_destory(NilpLong_t *nilp){
  if(nilp == NULL)
    return;

  delete nilp;

}

void nilpLongType1(NilpLong_t *nilp, long num, long size){
  Nilpotency<long> *obj;
  if(nilp == NULL)
    return;

  obj = static_cast<Nilpotency<long> *>(nilp->obj);
  obj->NilpType1(num, size);
}

void nilpLongShow(struct NilpLong *nilp){

  Nilpotency<long> *obj;
  if(nilp == NULL)
    return;

  obj = static_cast<Nilpotency<long> *>(nilp->obj);

  printf("Info ]> Show Nilptent Matrix Informations\n");
  printf("        diagPosition = %ld\n",obj->diagPosition);
  printf("        nbOne = %ld\n",obj->nbOne);
  printf("        matrix_size = %ld\n",obj->matrix_size);
  printf("        nilpotency = %ld\n",obj->nilpotency);
}

///

/*parMatrixSparse C wrapper*/

//real double int
struct parMatrixSparseDoubleInt{
  void *obj;
};

parMatrixSparseDoubleInt_t *newparMatrixSparseDoubleInt(){

  parMatrixSparseDoubleInt_t *mat;
  parMatrixSparse<double,int> *obj;

  mat = (parMatrixSparseDoubleInt *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<double,int>();
  mat->obj = obj;

  return mat;

}

void parMatrixSparseDoubleInt_destory(parMatrixSparseDoubleInt_t *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseDoubleInt_LocMatView(parMatrixSparseDoubleInt_t *mat){
  parMatrixSparse<double,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,int> *>(mat->obj);
  obj->LOC_MatView();

}
void GetLocalSizeDoubleInt(parMatrixSparseDoubleInt_t *mat, int *rs, int *cs){
  parMatrixSparse<double,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,int> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseDoubleInt_LocToCSR(parMatrixSparseDoubleInt_t *mat){
  parMatrixSparse<double,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,int> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;

}
void parMatrixSparseDoubleInt_LocGetCSRSize(parMatrixSparseDoubleInt_t *mat, int *size, int *size2){
  
  parMatrixSparse<double,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,int> *>(mat->obj);

  std::vector<int>::iterator it;
  int count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<double>::iterator itm;
  int count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= 0){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
  

}
void parMatrixSparseDoubleInt_LocGetCSRArrays(parMatrixSparseDoubleInt_t *mat, int size, int size2, int **rows, int **cols, double **vals){
  
  parMatrixSparse<double,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,int> *>(mat->obj);

  for(int i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = obj->CSR_loc->vals[i];
  }

}

void parMatrixSparseDoubleInt_smg2s(parMatrixSparseDoubleInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<double,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s<double,int>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;

}

void parMatrixSparseDoubleInt_nonsym_smg2s(parMatrixSparseDoubleInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<double,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s_nonsymmetric<double,int>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;
}

void parMatrixSparseDoubleInt_smg2s_advanace(parMatrixSparseDoubleInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, double *spectrum, double *init, MPI_Comm comm){
  parMatrixSparse<double,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s<double,int>(probSize, *nilpobj, lbandwidth,spectrum, init, comm);

  mat->obj = matobj;

}

void parMatrixSparseDoubleInt_nonsym_smg2s_advance(parMatrixSparseDoubleInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, dcomplex_t *spectrum, double *init, MPI_Comm comm){
  parMatrixSparse<double,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  std::complex<double> *spectrum2 = new std::complex<double>[probSize];
  for(int i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);    
  }
  matobj = smg2s_nonsymmetric<double,int>(probSize, *nilpobj, lbandwidth,spectrum2,init,comm);

  mat->obj = matobj;
}


// double long
struct parMatrixSparseDoubleLong{
  void *obj;
};


parMatrixSparseDoubleLong_t *newparMatrixSparseDoubleLong(){
  parMatrixSparseDoubleLong_t *mat;
  parMatrixSparse<double,__int64_t> *obj;

  mat = (parMatrixSparseDoubleLong *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<double,__int64_t>();
  mat->obj = obj;

  return mat;
}

void parMatrixSparseDoubleLong_destory(parMatrixSparseDoubleLong_t *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseDoubleLong_LocMatView(parMatrixSparseDoubleLong_t *mat){
  parMatrixSparse<double,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);
  obj->LOC_MatView();
}

void GetLocalSizeDoubleLong(parMatrixSparseDoubleLong_t *mat, __int64_t *rs, __int64_t *cs){
  parMatrixSparse<double,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseDoubleLong_LocToCSR(parMatrixSparseDoubleLong_t *mat){
  parMatrixSparse<double,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;
}

void parMatrixSparseDoubleLong_LocGetCSRSize(parMatrixSparseDoubleLong_t *mat, __int64_t *size, __int64_t *size2){
  parMatrixSparse<double,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);

  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<double>::iterator itm;
  __int64_t count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= 0){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
}

void parMatrixSparseDoubleLong_LocGetCSRArrays(parMatrixSparseDoubleLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, double **vals){
  parMatrixSparse<double,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);

  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = obj->CSR_loc->vals[i];
  }

}
void parMatrixSparseDoubleLong_smg2s(parMatrixSparseDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<double,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s<double,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;  
}
void parMatrixSparseDoubleLong_nonsym_smg2s(parMatrixSparseDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<double,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s_nonsymmetric<double,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;  
}



void parMatrixSparseDoubleLong_smg2s_advance(parMatrixSparseDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, double *spectrum, double *init, MPI_Comm comm){
  parMatrixSparse<double,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s<double,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,init, comm);

  mat->obj = matobj;  
}
void parMatrixSparseDoubleLong_nonsym_smg2s_advance(parMatrixSparseDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, dcomplex_t *spectrum, double *init, MPI_Comm comm){
  parMatrixSparse<double,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<double,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);


  std::complex<double> *spectrum2 = new std::complex<double>[probSize];
  for(__int64_t i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);    
  }


  matobj = smg2s_nonsymmetric<double,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum2,init,comm);

  mat->obj = matobj;  
}



////

//real float int
struct parMatrixSparseFloatInt{
  void *obj;
};

parMatrixSparseFloatInt_t *newparMatrixSparseFloatInt(){

  parMatrixSparseFloatInt_t *mat;
  parMatrixSparse<float,int> *obj;

  mat = (parMatrixSparseFloatInt_t *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<float,int>();
  mat->obj = obj;

  return mat;

}

void parMatrixSparseFloatInt_destory(parMatrixSparseFloatInt_t *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseFloatInt_LocMatView(parMatrixSparseFloatInt_t *mat){
  parMatrixSparse<float,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,int> *>(mat->obj);
  obj->LOC_MatView();

}
void GetLocalSizeFloatInt(parMatrixSparseFloatInt_t *mat, int *rs, int *cs){
  parMatrixSparse<float,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,int> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseFloatInt_LocToCSR(parMatrixSparseFloatInt_t *mat){
  parMatrixSparse<float,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,int> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;

}
void parMatrixSparseFloatInt_LocGetCSRSize(parMatrixSparseFloatInt_t *mat, int *size, int *size2){
  
  parMatrixSparse<float,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,int> *>(mat->obj);

  std::vector<int>::iterator it;
  int count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<float>::iterator itm;
  int count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= 0){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
  

}
void parMatrixSparseFloatInt_LocGetCSRArrays(parMatrixSparseFloatInt_t *mat, int size, int size2, int **rows, int **cols, float **vals){
  
  parMatrixSparse<float,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,int> *>(mat->obj);

  for(int i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = obj->CSR_loc->vals[i];
  }

}

void parMatrixSparseFloatInt_smg2s(parMatrixSparseFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<float,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s<float,int>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;

}

void parMatrixSparseFloatInt_nonsym_smg2s(parMatrixSparseFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<float,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s_nonsymmetric<float,int>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;
}


void parMatrixSparseFloatInt_smg2s_advance(parMatrixSparseFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, float *spectrum, float *init, MPI_Comm comm){
  parMatrixSparse<float,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s<float,int>(probSize, *nilpobj, lbandwidth,spectrum,init, comm);

  mat->obj = matobj;

}

void parMatrixSparseFloatInt_nonsym_smg2s_advance(parMatrixSparseFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, fcomplex_t *spectrum, float *init, MPI_Comm comm){
  parMatrixSparse<float,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  std::complex<float> *spectrum2 = new std::complex<float>[probSize];
  for(int i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);    
  }

  matobj = smg2s_nonsymmetric<float,int>(probSize, *nilpobj, lbandwidth,spectrum2,init,comm);

  mat->obj = matobj;
}


// float long
struct parMatrixSparseFloatLong{
  void *obj;
};


parMatrixSparseFloatLong_t *newparMatrixSparseFloatLong(){
  parMatrixSparseFloatLong_t *mat;
  parMatrixSparse<float,__int64_t> *obj;

  mat = (parMatrixSparseFloatLong *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<float,__int64_t>();
  mat->obj = obj;

  return mat;
}

void parMatrixSparseFloatLong_destory(parMatrixSparseFloatLong_t *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseFloatLong_LocMatView(parMatrixSparseFloatLong_t *mat){
  parMatrixSparse<float,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);
  obj->LOC_MatView();
}

void GetLocalSizeFloatLong(parMatrixSparseFloatLong_t *mat, __int64_t *rs, __int64_t *cs){
  parMatrixSparse<float,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseFloatLong_LocToCSR(parMatrixSparseFloatLong_t *mat){
  parMatrixSparse<float,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;
}

void parMatrixSparseFloatLong_LocGetCSRSize(parMatrixSparseFloatLong_t *mat, __int64_t *size, __int64_t *size2){
  parMatrixSparse<float,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);

  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<float>::iterator itm;
  __int64_t count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= 0){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
}

void parMatrixSparseFloatLong_LocGetCSRArrays(parMatrixSparseFloatLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, float **vals){
  parMatrixSparse<float,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);

  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = obj->CSR_loc->vals[i];
  }

}
void parMatrixSparseFloatLong_smg2s(parMatrixSparseFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<float,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s<float,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;  
}
void parMatrixSparseFloatLong_nonsym_smg2s(parMatrixSparseFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<float,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s_nonsymmetric<float,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;  
}


void parMatrixSparseFloatLong_smg2s_advance(parMatrixSparseFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, float *spectrum, float *init, MPI_Comm comm){
  parMatrixSparse<float,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s<float,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,init,comm);

  mat->obj = matobj;  
}
void parMatrixSparseFloatLong_nonsym_smg2s_advance(parMatrixSparseFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, fcomplex_t *spectrum, float *init, MPI_Comm comm){
  parMatrixSparse<float,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<float,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  std::complex<float> *spectrum2 = new std::complex<float>[probSize];
  for(__int64_t i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);    
  }

  matobj = smg2s_nonsymmetric<float,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum2,init,comm);

  mat->obj = matobj;  
}

//////

//double complex int
struct parMatrixSparseCmplxDoubleInt{
  void *obj;
};

parMatrixSparseCmplxDoubleInt_t *newparMatrixSparseCmplxDoubleInt(){

  parMatrixSparseCmplxDoubleInt_t *mat;
  parMatrixSparse<std::complex<double>,int> *obj;

  mat = (parMatrixSparseCmplxDoubleInt_t *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<std::complex<double>,int>();
  mat->obj = obj;

  return mat;

}

void parMatrixSparseCmplxDoubleInt_destory(parMatrixSparseCmplxDoubleInt *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseCmplxDoubleInt_LocMatView(parMatrixSparseCmplxDoubleInt *mat){
  parMatrixSparse<std::complex<double>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,int> *>(mat->obj);
  obj->LOC_MatView();

}
void GetLocalSizeCmplxDoubleInt(parMatrixSparseCmplxDoubleInt *mat, int *rs, int *cs){
  parMatrixSparse<std::complex<double>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,int> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseCmplxDoubleInt_LocToCSR(parMatrixSparseCmplxDoubleInt *mat){
  parMatrixSparse<std::complex<double>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,int> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;

}
void parMatrixSparseCmplxDoubleInt_LocGetCSRSize(parMatrixSparseCmplxDoubleInt *mat, int *size, int *size2){
  
  parMatrixSparse<std::complex<double>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,int> *>(mat->obj);

  std::vector<int>::iterator it;
  int count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<std::complex<double>>::iterator itm;
  int count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= std::complex<double>(0)){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
  

}
void parMatrixSparseCmplxDoubleInt_LocGetCSRArrays(parMatrixSparseCmplxDoubleInt *mat, int size, int size2, int **rows, int **cols, dcomplex_t **vals){
  
  parMatrixSparse<std::complex<double>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,int> *>(mat->obj);

  for(int i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = {(obj->CSR_loc->vals[i]).real(), (obj->CSR_loc->vals[i]).imag()};
  }

}

void parMatrixSparseCmplxDoubleInt_smg2s(parMatrixSparseCmplxDoubleInt *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<std::complex<double>,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<double>,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s<std::complex<double>,int>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;

}


void parMatrixSparseCmplxDoubleInt_smg2s_advance(parMatrixSparseCmplxDoubleInt *mat, int probSize, struct NilpInt *nilp, int lbandwidth, dcomplex *spectrum, dcomplex *init, MPI_Comm comm){
  parMatrixSparse<std::complex<double>,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<double>,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  std::complex<double> *spectrum2 = new std::complex<double>[probSize];

  std::complex<double> *init2 = new std::complex<double>[probSize];

  for(int i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);

  }

  matobj = smg2s<std::complex<double>,int>(probSize, *nilpobj, lbandwidth,spectrum2, init2,comm);

  mat->obj = matobj;

}


// double complex long
struct parMatrixSparseCmplxDoubleLong{
  void *obj;
};


parMatrixSparseCmplxDoubleLong_t *newparMatrixSparseCmplxDoubleLong(){
  parMatrixSparseCmplxDoubleLong_t *mat;
  parMatrixSparse<std::complex<double>,__int64_t> *obj;

  mat = (parMatrixSparseCmplxDoubleLong_t *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<std::complex<double>,__int64_t>();
  mat->obj = obj;

  return mat;
}

void parMatrixSparseCmplxDoubleLong_destory(parMatrixSparseCmplxDoubleLong_t *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseCmplxDoubleLong_LocMatView(parMatrixSparseCmplxDoubleLong_t *mat){
  parMatrixSparse<std::complex<double>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,__int64_t> *>(mat->obj);
  obj->LOC_MatView();
}

void GetLocalSizeDoubleLong(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t *rs, __int64_t *cs){
  parMatrixSparse<std::complex<double>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,__int64_t> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseDoubleLong_LocToCSR(parMatrixSparseCmplxDoubleLong_t *mat){
  parMatrixSparse<std::complex<double>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,__int64_t> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;
}

void parMatrixSparseCmplxDoubleLong_LocGetCSRSize(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t *size, __int64_t *size2){
  parMatrixSparse<std::complex<double>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,__int64_t> *>(mat->obj);

  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<std::complex<double>>::iterator itm;
  __int64_t count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= std::complex<double>(0)){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
}

void parMatrixSparseCmplxDoubleLong_LocGetCSRArrays(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, std::complex<double> **vals){
  parMatrixSparse<std::complex<double>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<double>,__int64_t> *>(mat->obj);

  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = {(obj->CSR_loc->vals[i]).real(), (obj->CSR_loc->vals[i]).imag()};
  }

}
void parMatrixSparseCmplxDoubleLong_smg2s(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<std::complex<double>,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<double>,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s<std::complex<double>,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;  
}


void parMatrixSparseCmplxDoubleLong_smg2s_advance(parMatrixSparseCmplxDoubleLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, dcomplex_t *spectrum, double *init, MPI_Comm comm){
  parMatrixSparse<std::complex<double>,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<double>,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);


  std::complex<double> *spectrum2 = new std::complex<double>[probSize];

  std::complex<double> *init2 = new std::complex<double>[probSize];

  for(__int64_t i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);

  }

  matobj = smg2s<std::complex<double>,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum2, init2,comm);

  mat->obj = matobj;  
}


////

//complex float int
struct parMatrixSparseCmplxFloatInt{
  void *obj;
};

parMatrixSparseCmplxFloatInt_t *newparMatrixSparseCmplxFloatInt(){

  parMatrixSparseCmplxFloatInt_t *mat;
  parMatrixSparse<std::complex<float>,int> *obj;

  mat = (parMatrixSparseCmplxFloatInt_t *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<std::complex<float>,int>();
  mat->obj = obj;

  return mat;

}

void parMatrixSparseCmplxFloatInt_destory(parMatrixSparseCmplxFloatInt_t *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseCmplxFloatInt_LocMatView(parMatrixSparseCmplxFloatInt_t *mat){
  parMatrixSparse<std::complex<float>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,int> *>(mat->obj);
  obj->LOC_MatView();

}
void GetLocalSizeCmplxFloatInt(parMatrixSparseCmplxFloatInt_t *mat, int *rs, int *cs){
  parMatrixSparse<std::complex<float>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,int> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseCmplxFloatInt_LocToCSR(parMatrixSparseCmplxFloatInt_t *mat){
  parMatrixSparse<std::complex<float>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,int> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;

}
void parMatrixSparseCmplxFloatInt_LocGetCSRSize(parMatrixSparseCmplxFloatInt_t *mat, int *size, int *size2){
  
  parMatrixSparse<std::complex<float>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,int> *>(mat->obj);

  std::vector<int>::iterator it;
  int count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<std::complex<float>>::iterator itm;
  int count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= std::complex<float>(0)){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
  

}
void parMatrixSparseCmplxFloatInt_LocGetCSRArrays(parMatrixSparseCmplxFloatInt_t *mat, int size, int size2, int **rows, int **cols, std::complex<float> **vals){
  
  parMatrixSparse<std::complex<float>,int> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,int> *>(mat->obj);

  for(int i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(int i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = {(obj->CSR_loc->vals[i]).real(), (obj->CSR_loc->vals[i]).imag()};
  }

}

void parMatrixSparseCmplxFloatInt_smg2s(parMatrixSparseCmplxFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<std::complex<float>,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<float>,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  matobj = smg2s<std::complex<float>,int>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;

}


void parMatrixSparseCmplxFloatInt_smg2s_advance(parMatrixSparseCmplxFloatInt_t *mat, int probSize, struct NilpInt *nilp, int lbandwidth, fcomplex_t *spectrum, float *init, MPI_Comm comm){
  parMatrixSparse<std::complex<float>,int> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<float>,int> *>(mat->obj);

  Nilpotency<int> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<int> *>(nilp->obj);

  std::complex<float> *spectrum2 = new std::complex<float>[probSize];

  std::complex<float> *init2 = new std::complex<float>[probSize];

  for(int i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);

  }

  matobj = smg2s<std::complex<float>,int>(probSize, *nilpobj, lbandwidth,spectrum2, init2,comm);

  mat->obj = matobj;

}


// complex float long
struct parMatrixSparseCmplxFloatLong{
  void *obj;
};


parMatrixSparseCmplxFloatLong_t *newparMatrixSparseCmplxFloatLong(){
  parMatrixSparseCmplxFloatLong_t *mat;
  parMatrixSparse<std::complex<float>,__int64_t> *obj;

  mat = (parMatrixSparseCmplxFloatLong_t *)malloc(sizeof(*mat));
  obj = new parMatrixSparse<std::complex<float>,__int64_t>();
  mat->obj = obj;

  return mat;
}

void parMatrixSparseCmplxFloatLong_destory(parMatrixSparseCmplxFloatLong_t *mat){
  if(mat == NULL)
    return;

  delete mat;
}

void parMatrixSparseCmplxFloatLong_LocMatView(parMatrixSparseCmplxFloatLong_t *mat){
  parMatrixSparse<std::complex<float>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,__int64_t> *>(mat->obj);
  obj->LOC_MatView();
}

void GetLocalSizeCmplxFloatLong(parMatrixSparseCmplxFloatLong_t *mat, __int64_t *rs, __int64_t *cs){
  parMatrixSparse<std::complex<float>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,__int64_t> *>(mat->obj);
  obj->GetLocalSize(*rs, *cs);
}

void parMatrixSparseCmplxFloatLong_LocToCSR(parMatrixSparseCmplxFloatLong_t *mat){
  parMatrixSparse<std::complex<float>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,__int64_t> *>(mat->obj);
  obj->Loc_ConvertToCSR();
  mat->obj = obj;
}

void parMatrixSparseCmplxFloatLong_LocGetCSRSize(parMatrixSparseCmplxFloatLong_t *mat, __int64_t *size, __int64_t *size2){
  parMatrixSparse<std::complex<float>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,__int64_t> *>(mat->obj);

  std::vector<__int64_t>::iterator it;
  __int64_t count = 0;

  for(it = obj->CSR_loc->rows.begin(); it != obj->CSR_loc->rows.end(); ++it){
      count++;
  }  

  std::vector<std::complex<float>>::iterator itm;
  __int64_t count2 = 0;
  for(itm = obj->CSR_loc->vals.begin(); itm != obj->CSR_loc->vals.end(); ++itm){
      if(*itm!= std::complex<float>(0)){
        count2++;
      }
  }

  *size = count;
  *size2 = count2;
}

void parMatrixSparseCmplxFloatLong_LocGetCSRArrays(parMatrixSparseCmplxFloatLong_t *mat, __int64_t size, __int64_t size2, __int64_t **rows, __int64_t **cols, fcomplex_t **vals){
  parMatrixSparse<std::complex<float>,__int64_t> *obj;
  if(mat == NULL)
    return;

  obj = static_cast<parMatrixSparse<std::complex<float>,__int64_t> *>(mat->obj);

  for(__int64_t i = 0; i < size; i++){
    *(*rows+i) = obj->CSR_loc->rows[i];
  }

  for(__int64_t i = 0; i < size2; i++){
    *(*cols+i) = obj->CSR_loc->cols[i];
    *(*vals+i) = {(obj->CSR_loc->vals[i]).real(), (obj->CSR_loc->vals[i]).imag()};
  }

}
void parMatrixSparseCmplxFloatLong_smg2s(parMatrixSparseCmplxFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, char *spectrum, MPI_Comm comm){
  parMatrixSparse<std::complex<float>,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<float>,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  matobj = smg2s<std::complex<float>,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum,comm);

  mat->obj = matobj;  
}



void parMatrixSparseCmplxFloatLong_smg2s_advance(parMatrixSparseCmplxFloatLong_t *mat, __int64_t probSize, struct NilpInt *nilp, __int64_t lbandwidth, fcomplex_t *spectrum, float *init, MPI_Comm comm){
  parMatrixSparse<std::complex<float>,__int64_t> *matobj;
  if(mat == NULL)
    return;

  matobj = static_cast<parMatrixSparse<std::complex<float>,__int64_t> *>(mat->obj);

  Nilpotency<__int64_t> *nilpobj;
  if(nilp == NULL)
    return;

  nilpobj = static_cast<Nilpotency<__int64_t> *>(nilp->obj);

  std::complex<float> *spectrum2 = new std::complex<float>[probSize];

  std::complex<float> *init2 = new std::complex<float>[probSize];

  for(__int64_t i = 0; i < probSize; i++){
    spectrum2[i].real(spectrum[i].real);
    spectrum2[i].imag(spectrum[i].imag);

  }

  matobj = smg2s<std::complex<float>,__int64_t>(probSize, *nilpobj, lbandwidth,spectrum2, init2,comm);

  mat->obj = matobj;  
}






