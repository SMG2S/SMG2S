#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

#include <smg2s-interface.hpp>
#include <parVector/parVectorMap.hpp>
#include <parVector/parVector.hpp>

template<typename T, typename S>
S checkNonSymmSpec(parVector<std::complex<Base<T>>, S> spec){
    auto nu = spec.GetUpperNeighbor(1);
    auto nl = spec.GetLowerNeighbor(1);
    auto nu_2 = spec.GetUpperNeighbor(2);
    auto nl_2 = spec.GetLowerNeighbor(2);

    auto lower_bound = spec.GetLowerBound();
    auto upper_bound = spec.GetUpperBound();

    auto *array = spec.GetArray();
    auto local_size = spec.GetLocalSize();

    int MyPID, nProcs;
    MPI_Comm_rank(spec.GetComm(), &MyPID);
    MPI_Comm_size(spec.GetComm(), &nProcs);

    S overlap = 0;

    S idx = 0;
    S step;
    S nzeros;

    if(MyPID == 0){
    	nzeros = 0;
    	while(idx < (local_size - 1) ){
    	    if(array[idx].imag() == 0){
    	    	step = 1;
    	    	nzeros++;
    	    }else{
    	    	
    	    	try{

    	    	    if(array[idx].imag() + array[idx+1].imag() != T(0) || array[idx].real() != array[idx+1].real()){
    	    	        throw 505;
    	    	    }

    	    	}catch(...){
    	    	    
    	    	    std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
    	    	    
    	    	}

    	    	step = 2;
    	    }
    	    idx += step;
    	}

    	if( (local_size - nzeros) % 2 != 0 ){ //should compare the element from next proc
    	    if(array[local_size - 1].imag() != 0){

    	    	try{

    	    	    if(array[local_size - 1].imag() + nl.imag() != Base<T>(0) || array[local_size - 1].real() != nl.real()){
    	    	        throw 505;
    	    	    }

    	    	}catch(...){
    	    	    std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
    	    	}
    	    }
    	}
    }

    if(MyPID == nProcs - 1){
    	nzeros = 0;
    	idx = local_size - 1;
    	while(idx > 0){
    	    if(array[idx].imag() == 0){
    	    	step = 1;
    	    	nzeros++;
	    }else{
	    	try{
	    	    if(array[idx].imag() + array[idx-1].imag() != Base<T>(0) || array[idx].real() != array[idx-1].real()){
    	    	        throw 505;
    	    	    }
	    	}catch(...){
    	    	    std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
	    	}

    	    	step = 2;
	    }
	    idx -= step;
    	}

    	if( (local_size - nzeros) % 2 != 0 ){ //should compare the element from next proc
    	    if(array[0].imag() != 0){

    	    	try{
    	    	    if(array[0].imag() + nu.imag() != Base<T>(0) || array[0].real() != nu.real()){
    	    	        throw 505;
    	    	        overlap = 1;
    	    	    }

    	    	}catch(...){
    	    	    std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
    	    	}
    	    }
    	}

    }
    //
    if(MyPID != 0 && MyPID != nProcs - 1){
    	
    	S size = local_size;

    	if(array[0].imag() == 0){
    	    idx = 1;
    	    size = local_size - 1;
    	}else{
    	    if(array[0].imag() + array[1].imag() != Base<T>(0) || array[0].real() != array[1].real()){
    	    	if(array[0].imag() + nu.imag() != Base<T>(0) || array[0].real() != nu.real()){
    	    	    try{
    	    	    	throw 505;
    	    	    }catch(...){
    	    	    	std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
    	    	    }
    	    	}else{
    	    	    if(nu_2.imag() + nu.imag() != Base<T>(0) || nu_2.real() != nu.real()){
    	    	    	idx = 1;
    	    	    	size = local_size - 1;
    	    	    	overlap = 1;
    	    	    }else{
    	    	    	try{
    	    	    	    throw 505;
    	    	    	}catch(...){
	    	    	    std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
    	    	    	}
    	    	    }
    	    	}    	    	
    	    }else{
    	    	idx = 2;
    	    	size = local_size - 2;
    	    }
    	}

    	nzeros = 0;
    	while(idx < (local_size - 1) ){
    	    if(array[idx].imag() == 0){
    	    	step = 1;
    	    	nzeros++;
    	    }else{
    	    	if(array[idx].imag() + array[idx+1].imag() != Base<T>(0) || array[idx].real() != array[idx+1].real()){
    	    	    try{
    	    	        throw 505;
    	    	    }catch(...){
	    	    	std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
    	    	    }
    	    	}
    	    	step = 2;
    	    }
    	    idx += step;
    	}

    	if( (size - nzeros) % 2 != 0 ){ //should compare the element from next proc
    	    if(array[local_size - 1].imag() != 0){
    	    	if(array[local_size - 1].imag() + nl.imag() != Base<T>(0) || array[local_size - 1].real() != nl.real()){
    	    	    try{
    	    	    	throw 505;
    	    	    }catch(...){
	    	    	std::cout << "SMG2S]> Caught Exception: input spectrum is invalid for non-symmetric matrix" << std::endl;
    	    	    }
    	    	}
    	    }
    	}    	
    }

    return overlap;
}


template<typename T, typename S>
parVector<T, S> specNonHerm(parVectorMap<S> index_map, std::string spectrum){
    parVector<T, S> spec =  parVector<T, S>(index_map);
    spec.ReadExtVec(spectrum);

    return spec;
}


template<typename T, typename S>
parVector<T, S> specNonSymm(parVectorMap<S> index_map, std::string spectrum){
    parVector<T, S> spec =  parVector<T, S>(index_map);
    spec.ReadExtVec(spectrum);

    return spec;
}


template<typename T, typename S>
parVector<std::complex<Base<T>>, S> specNonSymmCplex(parVectorMap<S> index_map, std::string spectrum){
    parVector<std::complex<Base<T>>, S> spec =  parVector<std::complex<Base<T>>, S>(index_map);
    spec.ReadExtVec(spectrum);

    return spec;
}

#endif