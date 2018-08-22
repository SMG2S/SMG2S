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
#include <complex>

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

		void VecAdd(parVector *v);
		void VecScale(T scale);
		T    VecDot(parVector *v);
		void ReadExtVec(std::string spectrum);
        void VecView();

		void RestoreArray(){};

		void specGen(std::string spectrum);
};

template<typename T,typename S>
parVector<T,S>::parVector(){
	array = NULL;
	array_size = 0;
	local_size = 0;
	index_map = NULL;
}

template<typename T,typename S>
parVector<T,S>::parVector(MPI_Comm ncomm, S lbound, S ubound)
{
	index_map = new parVectorMap<S>(ncomm, lbound, ubound);
	index_map->AddUser();

	local_size = index_map->GetLocalSize();
	array_size = index_map->GetLocTotSize();
	array = new T[array_size];
}

template<typename T,typename S>
parVector<T,S>::~parVector()
{
	if (index_map !=NULL){
		index_map->DeleteUser();
		if(index_map->GetUser() == 0){delete index_map;}
	}
	if (array != NULL){
		delete [] array;
	}
}

template<typename T,typename S>
S parVector<T,S>::GetLowerBound()
{
	if (index_map != NULL){
		return index_map->GetLowerBound();
	}
	else {return 0;}
}

template<typename T, typename S>
S parVector<T,S>::GetUpperBound()
{
	if (index_map != NULL){
		return index_map->GetUpperBound();
	}
	else {return 0;}
}

template<typename T, typename S>
S parVector<T,S>::GetLocalSize()
{
	return local_size;
}

template<typename T, typename S>
S parVector<T,S>::GetGlobalSize()
{
	if(index_map != NULL){
		return index_map->GetGlobalSize();
	}
	else {return 0;}
}

template<typename T, typename S>
S parVector<T,S>::GetArraySize()
{
	return array_size;
}

template<typename T, typename S>
S parVector<T,S>::Loc2Glob(S local_index)
{
	if ( index_map != NULL ) {
		return index_map -> Loc2Glob(local_index);
	}
	else return -1;
}

template<typename T, typename S>
S parVector<T,S>::Glob2Loc(S global_index)
{
	if ( index_map != NULL ) {
		return index_map ->Glob2Loc(global_index);
	} else return -1;
}

template<typename T, typename S>
void parVector<T,S>::AddValueLocal(S row, T value)
{
	if (row < array_size){
		array[row] = array[row] + value;
		//array[row] = value;
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
void parVector<T,S>::SetValueLocal(S row, T value)
{
	if (row < array_size){
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
	SetValueLocal(Glob2Loc(index), value);
}

template<typename T, typename S>
void parVector<T,S>::SetValuesGlobal(S nindex, S *rows, T *values)
{
	for(S i = 0; i < nindex; i++){
		SetValueLocal(Glob2Loc(rows[i]),values[i]);
	}
}

template<typename T, typename S>
void parVector<T,S>::SetTovalue(T value)
{
	for(S i= 0; i < array_size; i++) {
		array[i] = value;
	}
}

template<typename T, typename S>
void parVector<T,S>::SetToZero()
{
	T val = 0;
	SetTovalue(val);
}




template<typename T, typename S>
void parVector<T,S>::VecAdd(parVector<T,S> *v)
{
	if(array_size != v->array_size){std::cout << "vector size not coherant" << std::endl;}
	else{
		for(S i = 0; i < array_size; i++){
			array[i] = array[i] + v->array[i];
		}
	}
}

template<typename T, typename S>
void parVector<T,S>::VecScale(T scale)
{
	for(S i = 0; i < array_size; i++){
		array[i] = scale*array[i];
	}
}

template<typename T, typename S>
T parVector<T,S>::VecDot(parVector *v)
{
	T sum;
	for(S i = 0; i < array_size; i++){
		sum += array[i]*v->array[i];
	}

	return sum;
}
template<typename T, typename S>
void parVector<T,S>::VecView()
{
	int r;
	r = index_map->GetRank();
	if (r == 0){
		std::cout << "Parallel Vector View: " << std::endl << std::endl;
	}
	T *array = GetArray();
	S global;
	for(S i = 0; i < array_size; i++){
		global = Loc2Glob(i);
		std::cout << "[" << global << "]: " << array[i] << std::endl;
	}
}




/////
template<>
void parVector<std::complex<double>,int>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	int lower_bound = GetLowerBound();
	int upper_bound = GetUpperBound();

	int val1;

		
	double val2, val3;
	std::complex<double> val;

	while (std::getline(file,line)) {
		val1 = 0;
		val2 = 0.0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val2 >> val3;
		if (val1!= 0&&val2!= 0.0&&val3!= 0.0)
		{
			break ;
		}
	}

	while(std::getline(file, line))
	{
		val1 = 0; val2 = 0.0; val3 = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val2 >> val3;
		val1 = val1 - 1;
		val.real(val2);val.imag(val3);
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}
	}	
}

/////////

template<>
void parVector<std::complex<double>,__int64_t>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	__int64_t lower_bound = GetLowerBound();
	__int64_t upper_bound = GetUpperBound();

	__int64_t val1;

		
	double val2, val3;
	std::complex<double> val;

	while (std::getline(file,line)) {
		val1 = 0;
		val2 = 0.0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val2 >> val3;
		if (val1!= 0&&val2!= 0.0&&val3!= 0.0)
		{
			break ;
		}
	}

	while(std::getline(file, line))
	{
		val1 = 0; val2 = 0.0; val3 = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val2 >> val3;
		val1 = val1 - 1;
		val.real(val2);val.imag(val3);
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}
	}	
}

///////////////


template<>
void parVector<std::complex<float>,int>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	int lower_bound = GetLowerBound();
	int upper_bound = GetUpperBound();

	int val1;

		
	float val2, val3;
	std::complex<float> val;

	while (std::getline(file,line)) {
		val1 = 0;
		val2 = 0.0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val2 >> val3;
		if (val1!= 0&&val2!= 0.0&&val3!= 0.0)
		{
			break ;
		}
	}

	while(std::getline(file, line))
	{
		val1 = 0; val2 = 0.0; val3 = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val2 >> val3;
		val1 = val1 - 1;
		val.real(val2);val.imag(val3);
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}
	}	
}



template<>
void parVector<std::complex<float>,__int64_t>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	__int64_t lower_bound = GetLowerBound();
	__int64_t upper_bound = GetUpperBound();

	__int64_t val1;

		
	float val2, val3;
	std::complex<float> val;

	while (std::getline(file,line)) {
		val1 = 0;
		val2 = 0.0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val2 >> val3;
		if (val1!= 0&&val2!= 0.0&&val3!= 0.0)
		{
			break ;
		}
	}

	while(std::getline(file, line))
	{
		val1 = 0; val2 = 0.0; val3 = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val2 >> val3;
		val1 = val1 - 1;
		val.real(val2);val.imag(val3);
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}
	}	
}






/////
template<>
void parVector<double,int>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	int lower_bound = GetLowerBound();
	int upper_bound = GetUpperBound();

	int val1;

		
	double val;

	while (std::getline(file,line)) {
		val1 = 0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val;
		if (val1!= 0 &&val != 0.0)
		{
			break ;
		}
	}


	while(std::getline(file, line))
	{
		val1 = 0; val = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val;
		val1 = val1 - 1;
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}

	}	
}

/////////

template<>
void parVector<double,__int64_t>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	__int64_t lower_bound = GetLowerBound();
	__int64_t upper_bound = GetUpperBound();

	__int64_t val1;

		
	double val;

	while (std::getline(file,line)) {
		val1 = 0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val;
		if (val1!= 0 &&val != 0.0)
		{
			break ;
		}
	}


	while(std::getline(file, line))
	{
		val1 = 0; val = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val;
		val1 = val1 - 1;
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}

	}	
}





/////
template<>
void parVector<float,int>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	int lower_bound = GetLowerBound();
	int upper_bound = GetUpperBound();

	int val1;

		
	float val;

	while (std::getline(file,line)) {
		val1 = 0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val;
		if (val1!= 0 &&val != 0.0)
		{
			break ;
		}
	}


	while(std::getline(file, line))
	{
		val1 = 0; val = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val;
		val1 = val1 - 1;
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}

	}	
}

/////////

template<>
void parVector<float,__int64_t>::ReadExtVec(std::string spectrum)
{
	std::ifstream file(spectrum);
	std::string line;

	__int64_t lower_bound = GetLowerBound();
	__int64_t upper_bound = GetUpperBound();

	__int64_t val1;

		
	float val;

	while (std::getline(file,line)) {
		val1 = 0;
		std::stringstream linestream ( line ) ;
		linestream >> val1 >> val;
		if (val1!= 0 &&val != 0.0)
		{
			break ;
		}
	}


	while(std::getline(file, line))
	{
		val1 = 0; val = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> val;
		val1 = val1 - 1;
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val);
		}

	}	
}

#endif
