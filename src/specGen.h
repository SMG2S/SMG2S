#ifndef __SPECGEN_H__
#define __SPECGEN_H__

//#include "../parVector/parVector.cc"
#include "../parMatrix/parMatrixSparse.cc"

template<typename T, typename S>
void specGen(parVector<T,S> *vec){

	S    size;
	size = vec->GetGlobalSize();
	T    val;

	for(S i=0; i < size; i++){
		val = 2+i;
		vec->SetValueGlobal(i, val);
	}

}

#endif