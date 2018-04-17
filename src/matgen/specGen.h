#ifndef __SPECGEN_H__
#define __SPECGEN_H__

//#include "../parVector/parVector.cc"
#include "../../parMatrix/parMatrixSparse.cc"
#include "complex"
#include "../../utils/utils.h"

template<typename T, typename S>
void specGen(parVector<T,S> *vec){

	S    size;
	size = vec->GetGlobalSize();
	T    val;
	double    rx, ry, r;
        
/*
        val.real(0.22);
	val.imag(0.22);

                vec->SetValueGlobal(0, val);
*/
/*

	for(S i=0; i < size; i++){
		val.real((i+1)*10);
		val.imag(10+i*10);
		vec->SetValueGlobal(i, val);
	}
*/

	vec->ReadExtVec();
	vec->VecView();

}

#endif
