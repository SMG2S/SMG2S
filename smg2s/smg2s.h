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

#ifndef __SMG2S_H__
#define __SMG2S_H__

#include "../parVector/parVector.h"
#include "../parMatrix/parMatrixSparse.h"
#include "../parMatrix/initMat.h"

#include <math.h>
#include <complex.h>
#include <string>


#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif


template<typename T, typename S>
parMatrixSparse<T,S> nonherm(S probSize, Nilpotent<S> nilp, initMat<S> init, std::string spectrum, MPI_Comm comm){
	int world_size;
    int world_rank;

	double start, end;

    MPI_Comm_size(comm, &world_size);
    MPI_Comm_rank(comm, &world_rank);

    S span, lower_b, upper_b;

    if(sizeof(T) / sizeof(Base<T>) != 2){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-hermitian matrix, please use complex scalar type" << std::endl;
        }  
    }
    span = S(ceil(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    auto parVecMap = parVectorMap<S>(comm, lower_b, upper_b);

	auto spec = specNonHerm<T, S>(parVecMap, spectrum);

    return nonherm(probSize, nilp, init, spec);
}


template<typename T, typename S>
parMatrixSparse<T,S> nonherm(S probSize, Nilpotent<S> nilp, initMat<S> init, parVector<T, S> spec){

    if(sizeof(T) / sizeof(Base<T>) != 2){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-hermitian matrix, please use complex scalar type" << std::endl;
        }  
    }

    auto Am = parMatrixSparse<T, S>(spec.GetVecMap());

    Am.initMat(init.diag_l, init.diag_u, Base<T>(init.scale), T(0.0), Base<T>(init.sparsity) );

    nonherm(probSize, nilp, &Am, spec);

    return Am;
}


template<typename T, typename S>
void nonherm(S probSize, Nilpotent<S> nilp, parMatrixSparse<T,S> *Am, parVector<T,S> spec){
	int world_size;
    int world_rank;

    MPI_Comm comm = Am->GetComm();

	double start, end;

    MPI_Comm_size(comm, &world_size);
    MPI_Comm_rank(comm, &world_rank);


    if(sizeof(T) / sizeof(Base<T>) != 2){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-hermitian matrix, please use complex scalar type" << std::endl;
        }   
    }

    auto matAop = parMatrixSparse<T, S>(Am->GetMap());

    start = MPI_Wtime();

    Am->setSpecNonHerm(spec);
    matAop.setSpecNonHerm(spec);

    end = MPI_Wtime();

    double t2 = end - start;

    __int64_t my_factorielle_bornes = 1;

    my_factorielle_bornes = factorial(1,2*(nilp.getDegree()-1) );

    Am->MatScale((T)my_factorielle_bornes);

    for (S k=1; k<=2*(nilp.getDegree()-1); k++){
    	auto MA = matAop.MA(nilp);
    	auto AM = matAop.AM(nilp);
    	matAop.copy(AM);
    	//matAop.MatAYPX(AM, 0);
    	matAop.MatAXPY(MA, -1);
    	my_factorielle_bornes = factorial(k+1,2*(nilp.getDegree()-1));
    	Am->MatAXPY(matAop, (double)my_factorielle_bornes);
    }

    my_factorielle_bornes = factorial(1,2*(nilp.getDegree()-1));

    double fac = (double)my_factorielle_bornes;
    double inv = 1/fac;

    Am->MatScale((T)inv);    
}


///		
		

template<typename T, typename S>
parMatrixSparse<T,S> nonsymm(S probSize, Nilpotent<S> nilp, initMat<S> init, std::string spectrum, MPI_Comm comm){
	int world_size;
    int world_rank;

	double start, end;

    MPI_Comm_size(comm, &world_size);
    MPI_Comm_rank(comm, &world_rank);

    S span, lower_b, upper_b;

    if(sizeof(T) / sizeof(Base<T>) != 1){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-symmetric matrix, please use real scalar type" << std::endl;
        }  
    }
    span = S(ceil(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    auto parVecMap = parVectorMap<S>(comm, lower_b, upper_b);

    std::string word = "complex";

    std::ifstream file(spectrum);
    std::string line;

	std::getline(file,line);
	auto  pos = line.find(word);

	if ( pos != std::string::npos){
        auto spec = specNonSymmCplex<T, S>(parVecMap, spectrum);
        return nonsymmconj(probSize, nilp, init, spec);	
	}else{
	    auto spec = specNonSymm<T, S>(parVecMap, spectrum);
        return nonsymm(probSize, nilp, init, spec);	
	}
    
}


template<typename T, typename S>
parMatrixSparse<T,S> nonsymm(S probSize, Nilpotent<S> nilp, initMat<S> init, parVector<T, S> spec){

    if(sizeof(T) / sizeof(Base<T>) != 1 ){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-symmetric matrix, please use real scalar type" << std::endl;
        }      
    }
    auto Am = parMatrixSparse<T, S>(spec.GetVecMap());

    Am.initMat(init.diag_l, init.diag_u, Base<T>(init.scale), T(0.0), Base<T>(init.sparsity) );

    nonsymm(probSize, nilp, &Am, spec);

    return Am;
}


template<typename T, typename S>
void nonsymm(S probSize, Nilpotent<S> nilp, parMatrixSparse<T,S> *Am, parVector<T,S> spec){
	int world_size;
    int world_rank;

    MPI_Comm comm = Am->GetComm();

	double start, end;

    MPI_Comm_size(comm, &world_size);
    MPI_Comm_rank(comm, &world_rank);

    if(sizeof(T) / sizeof(Base<T>) != 1){

        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-symmetric matrix, please use real scalar type" << std::endl;
        }  

    }

    auto matAop = parMatrixSparse<T, S>(Am->GetMap());

    start = MPI_Wtime();

    Am->setSpecNonSymm(spec);
    matAop.setSpecNonSymm(spec);

    end = MPI_Wtime();

    double t2 = end - start;

    __int64_t my_factorielle_bornes = 1;

    my_factorielle_bornes = factorial(1,2*(nilp.getDegree()-1) );

    Am->MatScale((T)my_factorielle_bornes);

    for (S k=1; k<=2*(nilp.getDegree()-1); k++){
    	auto MA = matAop.MA(nilp);
    	auto AM = matAop.AM(nilp);
    	matAop.copy(AM);
    	//matAop.MatAYPX(AM, 0);
    	matAop.MatAXPY(MA, -1);
    	my_factorielle_bornes = factorial(k+1,2*(nilp.getDegree()-1));
    	Am->MatAXPY(matAop, (double)my_factorielle_bornes);
    }

    my_factorielle_bornes = factorial(1,2*(nilp.getDegree()-1));

    double fac = (double)my_factorielle_bornes;
    double inv = 1/fac;

    Am->MatScale((T)inv);
    
}


template<typename T, typename S>
parMatrixSparse<T,S> nonsymmconj(S probSize, Nilpotent<S> nilp, initMat<S> init, parVector<std::complex<T>, S> spec){

    if(sizeof(T) / sizeof(Base<T>) != 1 ){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-symmetric matrix, please use real scalar type" << std::endl;
        }  
    }
    auto Am = parMatrixSparse<Base<T>, S>(spec.GetVecMap());

    Am.initMat(init.diag_l, init.diag_u, Base<T>(init.scale), Base<T>(0.0), Base<T>(init.sparsity) );

    nonsymmconj(probSize, nilp, &Am, spec);

    return Am;
}


template<typename T, typename S>
void nonsymmconj(S probSize, Nilpotent<S> nilp, parMatrixSparse<T,S> *Am, parVector<std::complex<T>, S> spec){
	int world_size;
    int world_rank;

    MPI_Comm comm = Am->GetComm();

	double start, end;

    MPI_Comm_size(comm, &world_size);
    MPI_Comm_rank(comm, &world_rank);

    if(sizeof(T) / sizeof(Base<T>) != 1){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: this is for non-symmetric matrix, please use real scalar type" << std::endl;
        }  
    }

    auto matAop = parMatrixSparse<T, S>(Am->GetMap());

    start = MPI_Wtime();

    Am->setSpecNonSymmCmplx(spec);
    matAop.setSpecNonSymmCmplx(spec);

    end = MPI_Wtime();

    double t2 = end - start;

    __int64_t my_factorielle_bornes = 1;

    my_factorielle_bornes = factorial(1,2*(nilp.getDegree()-1) );

    Am->MatScale((T)my_factorielle_bornes);

    for (S k=1; k<=2*(nilp.getDegree()-1); k++){
    	auto MA = matAop.MA(nilp);
    	auto AM = matAop.AM(nilp);
    	matAop.copy(AM);
    	//matAop.MatAYPX(AM, 0);
    	matAop.MatAXPY(MA, -1);
    	my_factorielle_bornes = factorial(k+1,2*(nilp.getDegree()-1));
    	Am->MatAXPY(matAop, (double)my_factorielle_bornes);
    }

    my_factorielle_bornes = factorial(1,2*(nilp.getDegree()-1));

    double fac = (double)my_factorielle_bornes;
    double inv = 1/fac;

    Am->MatScale((T)inv);

}




#endif