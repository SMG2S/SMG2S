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

template<typename T, typename S>
void specGen(parVector<T, S> *spec, std::string spectrum){

  S    size = spec->GetGlobalSize();
  S lb = spec->GetLowerBound();
  S ub = spec->GetUpperBound();

  //complex values?
  int size1 = sizeof(T) / sizeof(Base<T>);
 
  std::vector<std::vector<Base<T> > > gen;

  T* specs;
  specs = new T[size];

   if (spectrum.compare(" ") == 0){
      if(spec->GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

      gen.resize(size1);
      for(int i = 0; i < size1; i++){
	gen[i].resize(size);
	for(S j=0; j < size; j++){
          gen[i][j] = j * 10 + 1;
      	}
      }	

      for(int i = 0; i < size1; i++){
        for(S j=0; j < size; j++){
	  reinterpret_cast<Base<T>*>(specs)[ size1 * j + i] = gen[i][j];
	}
      }

      for(S i=lb; i < ub; i++){
        spec->SetValueGlobal(i, specs[i]);
      }
   }
   else{

      spec->ReadExtVec(spectrum);
   }
}

template<typename T, typename S>
void matInit(parMatrixSparse<T,S> *Am, parMatrixSparse<T,S> *matAop, S probSize, S lbandwidth){

    T rnd;

    //T scale;

/*This part can be replaced by users provded func*/
    for(S i = 0; i < probSize; i++){
        for(S j = i - lbandwidth; j < i; j++){
            if(j >= 0){
              rnd = 1;
              Am->Loc_SetValue(i,j,rnd);
              matAop->Loc_SetValue(i,j,rnd);
            }
        }
    }
}

#endif
