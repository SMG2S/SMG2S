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

#ifndef __SPECGEN_H__
#define __SPECGEN_H__

#include "../parMatrix/parMatrixSparse.h"
#include "complex"
#include "../utils/utils.h"
#include <string>

template<>
void parVector<std::complex<double>,int>::specGen(std::string spectrum){

  int    size;
  size = GetGlobalSize();
  std::complex<double>    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(int i=0; i < size; i++){
        val.real(i*10+1);
        val.imag(i*10+1);
        SetValueGlobal(i, val);
      }
   }
   else{

      ReadExtVec(spectrum);
   }
}


template<>
void parVector<std::complex<double>,__int64_t>::specGen(std::string spectrum){

  __int64_t    size;
  size = GetGlobalSize();
  std::complex<double>    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(__int64_t i=0; i < size; i++){
        val.real(i*10+1);
        val.imag(i*10+1);
        SetValueGlobal(i, val);
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}



template<>
void parVector<double,int>::specGen(std::string spectrum){

	int    size;
	size = GetGlobalSize();
	double    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(int i=0; i < size; i++){
        val = i*10+1;
        SetValueGlobal(i, val);
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}


template<>
void parVector<double,__int64_t>::specGen(std::string spectrum){

  __int64_t    size;
  size = GetGlobalSize();
  double    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(__int64_t i=0; i < size; i++){
        val = i*10+1;
        SetValueGlobal(i, val);
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}


template<>
void parVector<float,int>::specGen(std::string spectrum){

  int    size;
  size = GetGlobalSize();
  float    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(int i=0; i < size; i++){
        val = i*10+1;
        SetValueGlobal(i, val);
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}


template<>
void parVector<float,__int64_t>::specGen(std::string spectrum){

  __int64_t    size;
  size = GetGlobalSize();
  float    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(__int64_t i=0; i < size; i++){
        val = i*10+1;
        SetValueGlobal(i, val);
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}


template<>
void parVector<std::complex<float>,int>::specGen(std::string spectrum){

  int    size;
  size = GetGlobalSize();
  std::complex<float>    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(int i=0; i < size; i++){
        val.real(i*10+1);
        val.imag(i*10+1);
        SetValueGlobal(i, val);
      }
   }
   else{

      ReadExtVec(spectrum);
   }
}


template<>
void parVector<std::complex<float>,__int64_t>::specGen(std::string spectrum){

  __int64_t    size;
  size = GetGlobalSize();
  std::complex<float>    val;

   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(__int64_t i=0; i < size; i++){
        val.real(i*10+1);
        val.imag(i*10+1);
        SetValueGlobal(i, val);
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}



template<typename T, typename S>
void matInit(parMatrixSparse<T,S> *Am, parMatrixSparse<T,S> *matAop, S probSize, S lbandwidth){

    T rnd;

    //T scale;

/*This part can be replaced by users provded func*/

/*
    if(std::is_same<T,std::complex<double> >::value || std::is_same<T,std::complex<float> >::value){
      scale.real(0.00001);
      scale.imag(0.0);      
    }else{
      scale = 0.00001;
    }
*/

    for(S i = 0; i < probSize; i++){
        for(S j = i - lbandwidth; j < i; j++){
            if(j >= 0){
              //rnd = scale * random<T,S>(0,10);
              rnd = 1;
              Am->Loc_SetValue(i,j,rnd);
              matAop->Loc_SetValue(i,j,rnd);
            }
        }
    }
}

template<typename T, typename S>
void matInit2(parMatrixSparse<T,S> *Am, parMatrixSparse<T,S> *matAop, S probSize, S lbandwidth){

    T rnd;

    //T scale;

/*This part can be replaced by users provded func*/

/*
    if(std::is_same<T,std::complex<double> >::value || std::is_same<T,std::complex<float> >::value){
      scale.real(0.00001);
      scale.imag(0.0);      
    }else{
      scale = 0.00001;
    }
*/

    for(S i = 0; i < probSize; i++){
        for(S j = i - lbandwidth; j < i - 1; j++){
            if(j >= 0){
              //rnd = scale * random<T,S>(0,10);
              rnd = 1;
              Am->Loc_SetValue(i,j,rnd);
              matAop->Loc_SetValue(i,j,rnd);
            }
        }
    }
}

#endif
