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

#endif
