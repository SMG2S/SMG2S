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

   Part of basic data structures' implementation of this file refers to this technical report 
   (http://orbit.dtu.dk/files/51272329/tr12_10_Alexandersen_Lazarov_Dammann_1.pdf)
*/

#ifndef __PAR_VECTOR_H__
#define __PAR_VECTOR_H__

#include <mpi.h> // Input/output
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <string>

#include "parVectorMap.h"
#include "../utils/utils.h"

template<typename T, typename S>
class parVector{
	private:
		T	*array;
		S	array_size;
		S	local_size;

		parVectorMap<S> *index_map;

	public:
		parVector();
		parVector(MPI_Comm ncomm, S lbound, S ubound);
		~parVector();

		parVectorMap<S> *GetVecMap(){return index_map;};

		S GetLowerBound();
		S GetUpperBound();
		S GetGlobalSize();
		S GetLocalSize();
		S GetArraySize();
		T *GetArray(){return array;};

		S Loc2Glob(S local_index);
		S Glob2Loc(S global_index);
		
		void AddValueLocal(S row, T value);
		void AddValuesLocal(S nindex, S *rows, T *values);

		void SetValueLocal(S row, T value);
		void SetValuesLocal(S nindex, S *rows, T *values);
		void SetValueGlobal(S index, T value);
		void SetValuesGlobal(S nindex, S *rows, T *values);

		void SetTovalue(T value);
		void SetToZero();
		void SetRandomValues(T min, T max);
		
		void VecAdd(parVector *v);
		void VecScale(T scale);
		T    VecDot(parVector *v);
		void ReadExtVec();
                void VecView();

		void RestoreArray(){};
};

#endif
