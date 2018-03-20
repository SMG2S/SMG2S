#include <mpi.h>
#include <iostream> 
#include <fstream> 
#include <sstream>
#include <string>

#include "parVector.h"
#include "parVectorMap.cc"
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
void parVector<T,S>::AddValueLocal(S row, T value)
{
	if (row < array_size){
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

template<typename T,typename S>
void parVector<T,S>::SetRandomValues(T min, T max)
{
        for(S i= 0; i < array_size; i++) {
                array[i] = random_unint(min, max);
        }
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
	}
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
template<typename T, typename S>
void parVector<T,S>::ReadExtVec()
{
	std::ifstream file("vector.vec");
	std::string line;

	S lower_bound = GetLowerBound();
	S upper_bound = GetUpperBound();

	S val1, dummy;
	T val2;

	int quit;

	// Read past first few lines of input file that are not numbers

	while (std::getline(file,line)) {
		val1 = 0; 
		val2 = 0.0;
		std::stringstream linestream ( line ) ; 
		linestream >> val1 >> dummy >> val2;
		if (dummy!= 0&&val1!= 0&&val2!= 0.0) 
		{
			break ;
		} 
	}

	// Start reading in coordinates and if local  > add to vector

	while(std::getline(file, line))
	{
		val1 = 0; val2 = 0.0;

		std::stringstream linestream(line);

		linestream >> val1 >> dummy >> val2;
		val1 = val1 - 1;
		
		if((val1 >= lower_bound) && (val1 < upper_bound)){
			AddValueLocal(index_map->Glob2Loc(val1),val2);
		}

	}	
}

