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

#ifndef __PAR_VECTOR_H__
#define __PAR_VECTOR_H__

#include <mpi.h> // Input/output
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>

#include "parVectorMap.h"
#include "../utils/utils.h"
#include "../utils/MPI_DataType.h"

template<typename T, typename S>
class parVector{
    private:
    	T	*array;
	S	local_size;
	S       global_size;
	MPI_Comm  	comm;
	parVectorMap<S> index_map;

    public:
	parVector();
	parVector(MPI_Comm ncomm, S lbound, S ubound);
	parVector(parVectorMap<S> map);
	~parVector();

	parVectorMap<S> GetVecMap(){return index_map;};

	S GetLowerBound(){return index_map.GetLowerBound();};
	S GetUpperBound(){return index_map.GetUpperBound();};
	S GetGlobalSize(){return global_size;};
	S GetLocalSize(){return local_size;};
	S GetRank(){return index_map.GetRank();};
	T* GetArray(){return array;};
	MPI_Comm GetComm(){return comm;};

	S Loc2Glob(S local_index){return index_map.Loc2Glob(local_index);};
	S Glob2Loc(S global_index){return index_map.Glob2Loc(global_index);};

	void SetToValue(T value);
	void SetToZero();		
	void VecView();
	void SetValueLocal(S row, T value);
	void SetValuesLocal(S nindex, S *rows, T *values);
	void SetValueGlobal(S index, T value);
	void SetValuesGlobal(S nindex, S *rows, T *values);		
	void AddValueLocal(S row, T value);
	void AddValuesLocal(S nindex, S *rows, T *values);
	void VecAdd(parVector v);

	void VecScale(T scale);
	T    VecDot(parVector v);
	void ReadExtVec(std::string spectrum);
        

};

template<typename T,typename S>
parVector<T,S>::parVector(){}


template<typename T,typename S>
parVector<T,S>::parVector(MPI_Comm ncomm, S lbound, S ubound)
{
    MPI_Comm_dup(ncomm, &comm);
    index_map = parVectorMap<S>(ncomm, lbound, ubound);
    local_size = index_map.GetLocalSize();
    global_size = index_map.GetGlobalSize();
    array = new T[local_size];
}

template<typename T,typename S>
parVector<T,S>::parVector(parVectorMap<S> map)
{
    index_map = map;
    comm = index_map.GetCurrentComm();
    local_size = index_map.GetLocalSize();
    global_size = index_map.GetGlobalSize();
    array = new T[local_size];
}

template<typename T,typename S>
parVector<T,S>::~parVector(){}


template<typename T, typename S>
void parVector<T,S>::SetToValue(T value)
{
    for(S i= 0; i < local_size; i++) {
    	array[i] = value;
    }
}

template<typename T, typename S>
void parVector<T,S>::SetToZero()
{
    T val = 0;
    SetToValue(val);
}

template<typename T, typename S>
void parVector<T,S>::VecView()
{
    int r;
    r = index_map.GetRank();
    if (r == 0){
	std::cout << "Parallel Vector View: " << std::endl << std::endl;
    }
    T *array = GetArray();
    S global;
    for(S i = 0; i < local_size; i++){
    	global = Loc2Glob(i);
	std::cout << "[" << global << "]: " << array[i] << std::endl;
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValueLocal(S row, T value)
{
    if (row < local_size){
	array[row] = value;
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValuesLocal(S nindex, S *rows, T *values)
{
    for(S i = 0; i < nindex; i++){
    	SetValueLocal(rows[i],values[i]);
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValueGlobal(S index, T value)
{
    int lower_bound = GetLowerBound();
    int upper_bound = GetUpperBound();
    if((index >= lower_bound) && (index < upper_bound)){
    	SetValueLocal(index_map.Glob2Loc(index), value);
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValuesGlobal(S nindex, S *rows, T *values)
{
    for(S i = 0; i < nindex; i++){
	SetValueLocal(Glob2Loc(rows[i]),values[i]);
    }
}

template<typename T, typename S>
void parVector<T,S>::AddValueLocal(S row, T value)
{
    if (row < local_size){
	array[row] = array[row] + value;
    }
}

template<typename T, typename S>
void parVector<T,S>::AddValuesLocal(S nindex, S *rows, T *values)
{
    for(S i = 0; i < nindex; i++){
	AddValueLocal(rows[i],values[i]);
    }
}

template<typename T, typename S>
void parVector<T,S>::VecAdd(parVector<T,S> v)
{
    if(local_size != v.GetLocalSize()){
	std::cout << "vector size not coherant" << std::endl;
    }
    else{
	for(S i = 0; i < local_size; i++){
	    array[i] = array[i] + v.array[i];
	}
    }
}

template<typename T, typename S>
void parVector<T,S>::VecScale(T scale)
{
    for(S i = 0; i < local_size; i++){
	array[i] = scale*array[i];
    }
}

template<typename T, typename S>
T parVector<T,S>::VecDot(parVector v)
{
    T sum;
    for(S i = 0; i < local_size; i++){
	sum += array[i]*v.array[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, getMPI_Type<T>(), MPI_SUM, comm);
    return sum;
}


template<typename T, typename S>
void parVector<T, S>::ReadExtVec(std::string spectrum)
{
    std::ifstream file(spectrum);
    std::string line;

    int lower_bound = GetLowerBound();
    int upper_bound = GetUpperBound();

    int size1 = sizeof(T) / sizeof(Base<T>);

    S idx;

    Base<T> in_vals[size1];

    T val;

    while (std::getline(file,line)) {
    	idx = 0;
        for(int i = 0; i < size1; i++){
            in_vals[i] = 0.0;
	}

        std::stringstream linestream ( line ) ;
        linestream >> idx;
	for(int i = 0; i < size1; i++){
            linestream >> in_vals[i];
        }

	int cnt = 0;
        if (idx!= 0)
        {
            for(int i = 0; i < size1; i++){
            	if(in_vals[i] != 0){
            	    cnt ++;
		}
            }
            if(cnt == 2){
                break ;
            }
        }
    }

	
    while (std::getline(file,line)) {
	idx = 0;
        for(int i = 0; i < size1; i++){
            in_vals[i] = 0.0;
        }

        std::stringstream linestream ( line ) ;
        linestream >> idx;
        for(int i = 0; i < size1; i++){
            linestream >> in_vals[i];
        }
        idx = idx - 1;

        for(int i = 0; i < size1; i++){
            reinterpret_cast<Base<T>(&)[size1]>(val)[i] = in_vals[i];
	}

	if((idx >= lower_bound) && (idx < upper_bound)){
	    AddValueLocal(index_map->Glob2Loc(idx),val);
	}		
    }

}

#endif
