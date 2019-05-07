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

#ifndef __SPECGEN_NONSYMMETRIC_H__
#define __SPECGEN_NONSYMMETRIC_H__

#include "../parMatrix/parMatrixSparse.h"
#include "complex"
#include "../utils/utils.h"
#include <string>

/*Non symmetric case*/

template<>
void parVector<std::complex<double>,int>::specGen2(std::string spectrum){

  int    size;
  size = GetGlobalSize();
  std::complex<double>    val;
  std::complex<double>    val2;


   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(int i=0; i < size; i=i+2){
        if( i == 2 ){
            val.real(i*1+3);
            val.imag(0);
            SetValueGlobal(i, val);
            SetValueGlobal(i+1, -val);
        }
        else{
          val.real(i*1+3);
          val.imag(i*2+1);
          val2.real(i*1+3);
          val2.imag(-i*2-1);
          SetValueGlobal(i, val);
          SetValueGlobal(i+1, val2);
        }
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}

template<>
void parVector<std::complex<float>,int>::specGen2(std::string spectrum){

  int    size;
  size = GetGlobalSize();
  std::complex<float>    val;
  std::complex<float>    val2;


   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(int i=0; i < size; i=i+2){
        if( i == 2 ){
            val.real(i*1+3);
            val.imag(0);
            SetValueGlobal(i, val);
            SetValueGlobal(i+1, -val);
        }
        else{
          val.real(i*1+3);
          val.imag(i*2+1);
          val2.real(i*1+3);
          val2.imag(-i*2-1);
          SetValueGlobal(i, val);
          SetValueGlobal(i+1, val2);
        }
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}

template<>
void parVector<std::complex<double>,__int64_t>::specGen2(std::string spectrum){

  __int64_t    size;
  size = GetGlobalSize();
  std::complex<double>    val;
  std::complex<double>    val2;


   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(__int64_t i=0; i < size; i = i + 2){
        if( i == 2 ){
            val.real(i*1+3);
            val.imag(0);
            SetValueGlobal(i, val);
            SetValueGlobal(i+1, -val);
        }
        else{
          val.real(i*1+3);
          val.imag(i*2+1);
          val2.real(i*1+3);
          val2.imag(-i*2-1);
          SetValueGlobal(i, val);
          SetValueGlobal(i+1, val2);
        }
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}

template<>
void parVector<std::complex<float>,__int64_t>::specGen2(std::string spectrum){

  __int64_t    size;
  size = GetGlobalSize();
  std::complex<float>    val;
  std::complex<float>    val2;


   if (spectrum.compare(" ") == 0){
      if(GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      for(int i=0; i < size; i=i+2){
        if( i == 2 ){
            val.real(i*1+3);
            val.imag(0);
            SetValueGlobal(i, val);
            SetValueGlobal(i+1, -val);
        }
        else{
          val.real(i*1+3);
          val.imag(i*2+1);
          val2.real(i*1+3);
          val2.imag(-i*2-1);
          SetValueGlobal(i, val);
          SetValueGlobal(i+1, val2);
        }
      }
   }
   else{
      ReadExtVec(spectrum);
   }
}


template<typename T, typename S>
void matInit2(parMatrixSparse<T,S> *Am, parMatrixSparse<T,S> *matAop, S probSize, S lbandwidth, parVector<std::complex<double>,S> *spec){

    T rnd;
    std::complex<double> *array;

    array = spec->GetArray();
    S loc_indx;

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
              rnd = 0.01;
              Am->Loc_SetValue(i,j,rnd);
              matAop->Loc_SetValue(i,j,rnd);
            }
        }
    }    


    for(S i = 0; i < probSize; i++){
      
      loc_indx = spec->Glob2Loc(i);

      Am->Loc_SetValue(i,i,array[loc_indx].real());
      matAop->Loc_SetValue(i,i,array[loc_indx].real());

    }
    
    for(S i = 0; i < probSize - 1; i = i + 2){
      
      loc_indx = spec->Glob2Loc(i);
      
      if(array[loc_indx].imag() != 0){
        Am->Loc_SetValue(i,i+1,std::abs(array[loc_indx].imag()));
        matAop->Loc_SetValue(i,i+1,std::abs(array[loc_indx].imag()));
        Am->Loc_SetValue(i+1,i,-std::abs(array[loc_indx].imag()));
        matAop->Loc_SetValue(i+1,i,-std::abs(array[loc_indx].imag()));
      }
    }
    
}

#endif
