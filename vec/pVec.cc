#include <mpi.h> // Input/output
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <string>

#include <pVecMap.h>

template<typename T, typename S>
pVec::pVec(){
	 array = NULL;
	 array_size = 0;
	 local_size = 0;
	 index_map = NULL;
}
template<typename T, typename S>
pVec::pVec(MPI_Comm ncomm, int lbound, int ubound){
	 index_map = new pVecMap ( ncomm , lbound , ubound ) ; 
	 index_map >AddUser();
	
	local_size = index_map->GetLocalSize();

	array_size = index_map->GetLocTotSize();
	
	array = new S[array_size];

}

template<typename T, typename S>
pVec::~pVec(){
	if (index_map !=NULL){
		index_map->DeleteUser();
		if(index_map->GetUser() == 0){delete index_map;}
	}
	if (array != NULL){
		delete [] array;
	}
}

template<typename T, typename S>
S pVec::GetLowerBound(){
	if (index_map != NULL){
		return index_map->GetLowerBound();
	}
	else {return 0;}
}

template<typename T, typename S>
S pvec::GetUpperBound(){
	if (index_map != NULL){
		return index_map->GetUpperBound();
	}
	else {return 0;}
}

template<typename T, typename S>
S pVec::GetLocalSize(){
	return local_size;
}

template<typename T, typename S>
S pVec:GetGlocalSize(){
	if(index_map != NULL){
		return index_map->GetGlobalSize();
	}
	else {return 0;}

template<typename T, typename S>
S pVec::GetArraySize(){
	return array_size;
}

template<typename T, typename S>
void pVec::AddValueLocal(S row, T value){
	if (row < array_size){
		array[row] = array[row] + value;
	}
}

template<typename T, typename S>
void pVec::AddValuesLocal(S nindex, S *rows, T *values)
{
	for(S i = 0; i < nindex; i++){
		AddValueLocal(rows[i],values[i]);
	}
}

template<typename T, typename S>
void pVec::SetTovalue(T value){
	for(S i= 0; i < array_size; i++) { array[i] = value;}
 
}
template<typename T, typename S>
void SetToZero(){
	SetToValue<T,S>(0.0);
}

template<typename T, typename S>
S pVec::Loc2Glob(S local_index){
	if ( index_map != NULL ) {
		return index_map >Loc2Glob(local_index);
	} else {return  -1; }
}

template<typename T, typename S>
S pVec::Glob2Loc(S global_index){
	 if ( index_map != NULL ) {
		return index_map >Glob2Loc(global_index);
	} else {return -1;}
}

	void ReadExtVec();
	

	void RestoreArray(){};

};
