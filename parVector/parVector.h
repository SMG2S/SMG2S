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
