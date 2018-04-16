#ifndef __SPECGEN_H__
#define __SPECGEN_H__

//#include "../parVector/parVector.cc"
#include "../../parMatrix/parMatrixSparse.cc"
#include <complex>

template<typename T, typename S>
void specGen(parVector<T,S> *vec){

	S    size;
	size = vec->GetGlobalSize();
	T    val;

	for(S i=0; i < size; i++){
		val.imag(1) ;
		val.real(2) ;
		vec->SetValueGlobal(i, val);
	}

}

#endif