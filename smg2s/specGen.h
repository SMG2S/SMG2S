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

#ifndef __SPECGEN_H__
#define __SPECGEN_H__

//#include "../parVector/parVector.cc"
#include "../parMatrix/parMatrixSparse.h"
#include "complex"
#include "../utils/utils.h"
#include <string>

template<typename T, typename S>
void specGen(parVector<T,S> *vec, std::string spectrum){

	S    size;
	size = vec->GetGlobalSize();
	T    val;
	double    rx, ry, r;

   if (spectrum.compare(" ") == 0){
      if(vec->GetVecMap()->GetRank() == 0){
         printf("Info ]> Do not provide the outside given spectrum file, using the internel function to generate them.\n");
      }

/*This part can be replaced by users provded func*/
      for(S i=0; i < size; i++){
#ifdef __USE_COMPLEX__
         val.real(i*10+1);
         val.imag(i*10+1);
#else
	val = i*10+1;
#endif
         vec->SetValueGlobal(i, val);
      }
   }
   else{

      vec->ReadExtVec(spectrum);
   }
}

template<typename T, typename S>
void matInit(parMatrixSparse<T,S> *Am, parMatrixSparse<T,S> *matAop, S probSize, S lbandwidth){

    T rnd;
/*This part can be replaced by users provded func*/

    for(S i = 0; i < probSize; i++){
        for(S j = i - lbandwidth; j < i; j++){
            if(j >= 0){
               rnd = 0.00001*random<T,S>(0,10);
              Am->Loc_SetValue(i,j,rnd);
              matAop->Loc_SetValue(i,j,rnd);
            }
        }
    }
}

#endif
