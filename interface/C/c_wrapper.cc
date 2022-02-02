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



// double long
struct parMatrixSparseDoubleLong;
typedef struct parMatrixSparseDoubleLong parMatrixSparseDoubleLong_t;

parMatrixSparseDoubleLong_t *newparMatrixSparseDoubleLong();
void parMatrixSparseDoubleLong_destory(parMatrixSparseDoubleLong_t *mat);
void parMatrixSparseDoubleLong_LocMatView(parMatrixSparseDoubleLong_t *mat);
void GetLocalSizeDoubleLong(parMatrixSparseDoubleLong_t *mat, long *rs, long *cs);
void parMatrixSparseDoubleLong_LocToCSR(parMatrixSparseDoubleLong_t *mat);
void parMatrixSparseDoubleLong_LocGetCSRSize(parMatrixSparseDoubleLong_t *mat, long *size, long *size2);
void parMatrixSparseDoubleLong_LocGetCSRArrays(parMatrixSparseDoubleLong_t *mat, long size, long size2, long **rows, long **cols, double **vals);
void parMatrixSparseDoubleLong_smg2s(parMatrixSparseDoubleLong_t *mat, long probSize, struct NilpInt *nilp, long lbandwidth, char *spectrum, MPI_Comm comm);
void parMatrixSparseDoubleLong_nonsym_smg2s(parMatrixSparseDoubleLong_t *mat, long probSize, struct NilpInt *nilp, long lbandwidth, char *spectrum, MPI_Comm comm);








