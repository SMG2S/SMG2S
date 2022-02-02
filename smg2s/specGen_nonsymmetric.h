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


template <typename T, typename S>
void specGen2(parVector<std::complex<T>, S> *spec, std::complex<Base<T>> *spectrum){
  S    size = spec->GetGlobalSize();
  S lb = spec->GetLowerBound();
  S ub = spec->GetUpperBound();
  
  for(S i=lb; i < ub; i++){
    spec->SetValueGlobal(i, spectrum[i]);
  }
}


template <typename T, typename S>
void specGen2(parVector<std::complex<T>, S> *spec, std::string spectrum){
    int size;
    size = spec->GetGlobalSize();
    int lb = spec->GetLowerBound();
    int ub = spec->GetUpperBound();

    std::complex<T>    val;
    std::complex<T> val2;

    if (spectrum.compare(" ") == 0){
        if(spec->GetVecMap()->GetRank() == 0){
            printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
        }
      
        if(size % 2 != 0){
            if(spec->GetVecMap()->GetRank() == 0){
                fprintf(stderr, "ERROR ]> For using default spectrum generator, the size of matrix should be even\n");
            }
            return;
        }

        for(S i = 0; i < size; i = i + 2){
          if(i >= lb && i < ub){
            val.real(i*1+3);
            val.imag(i*2+1);
            val2.real(i*1+3);
            val2.imag(-i*2-1);
            spec->SetValueGlobal(i, val);
            spec->SetValueGlobal(i+1, val2);
          }

        }


    }else{
        spec->ReadExtVec(spectrum);
    }
}


template<typename T, typename S>
void setSpectrum(parMatrixSparse<T,S> *Am, S probSize, parVector<std::complex<T>,S> *spec){

    int size1 = sizeof(T) / sizeof(Base<T>);

    if(size1 == 2){
        printf("Info ]> For generating non-Symmetric matrices, the matrices should be real\n");
        return;
    }

    std::complex<Base<T>> *array;

    array = spec->GetArray();
    S loc_indx;

    for(S i = 0; i < probSize; i++){
      
      loc_indx = spec->Glob2Loc(i);

      Am->Loc_SetValue(i,i,array[loc_indx].real());

    }
    
    for(S i = 0; i < probSize - 1; i = i + 2){
      
      loc_indx = spec->Glob2Loc(i);
      
      if(array[loc_indx].imag() != 0){
        Am->Loc_SetValue(i,i+1,std::abs(array[loc_indx].imag()));
        Am->Loc_SetValue(i+1,i,-std::abs(array[loc_indx].imag()));
      }
    }
}

template<typename T, typename S>
void matInit2(parMatrixSparse<T,S> *Am, parMatrixSparse<T,S> *matAop, S probSize, S lbandwidth, T *init){

    int size1 = sizeof(T) / sizeof(Base<T>);

    if(size1 == 2){
        printf("Info ]> For generating non-Symmetric matrices, the matrices should be real\n");
        return;
    }

    S cnt = 0;

    
    for(S i = 0; i < probSize; i++){
        S start = 0 > i - lbandwidth ? 0 : i - lbandwidth;
        for(S j = start; j < i - 1; j++){
            Am->Loc_SetValue(i,j,init[cnt]);
            matAop->Loc_SetValue(i,j,init[cnt]);
            cnt++;
        }
    }    
    
}


template<typename T, typename S>
void matInit2(parMatrixSparse<T,S> *Am, parMatrixSparse<T,S> *matAop, S probSize, S lbandwidth){

    int size1 = sizeof(T) / sizeof(Base<T>);

    if(size1 == 2){
        printf("Info ]> For generating non-Symmetric matrices, the matrices should be real\n");
        return;
    }

    T rnd;
    
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
    
}

#endif
