/*
MIT License
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France
Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
                                     Materials,  Forschungszentrum Juelich GmbH.
                                     
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

#include <parVector/parVectorMap.hpp>
#include <parVector/parVector.hpp>
#include <parMatrix/parMatrixSparse.hpp>
#include <smg2s/initMat.hpp>
#include <smg2s/nilpotent.hpp>
#include <smg2s/spectrum.hpp>

#include <math.h>

/** @defgroup group2 SMG2S
 *  This module relates to multiple implementation of SMG2S to generate a non-Symmetric/Hermitian sparse matrix with it. 
 *  @{
 */
//! Generating a non-Hermitian sparse matrix using the spectrum stored in local file
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] init a initMat object which determines the way of initialization of matrix, which will be further operated by SMG2S
  * @param[in] spectrum path and file name of a local file which provides the spectrum
  * @param[in] comm the working MPI communicator
  
  - The distribution of sparse matrix and vector is determined internally by this function. If you can cosider
  to use your own distribution scheme, please consider other two implementations of SMG2S for non-Hermitian matrix.
*/  
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

//! Generating a non-Hermitian sparse matrix using the spectrum stored in a parVector object
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] init a initMat object which determines the way of initialization of matrix, which will be further operated by SMG2S
  * @param[in] spec a parVector object which stores the given spectrum by users
  
  - The distribution of generated matrix across MPI procs is the same as the one of `spec` which provides the spectrum.
*/ 
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


//! Generating a non-Hermitian sparse matrix using the spectrum stored in a parVector object and user-provided initial matrix stored in a parMatrixSparse object
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] Am a parMatrixSparse provided by the users as a initial matrix of SMG2S, it can only be any sparse lower-triangular matrix
  * @param[in] spec a parVector object which stores the given spectrum by users
  * @param[out] Am for the output, it is overwritten by the generated non-Hermtian matrix
  
  - The distribution of generated matrix across MPI procs is the same as the one of `spec` which provides the spectrum.
  - The distribution scheme of `Am` and `spec` should be the same.
*/ 
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
        
//! Generating a non-Symmetric sparse matrix using the spectrum stored in local file
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] init a initMat object which determines the way of initialization of matrix, which will be further operated by SMG2S
  * @param[in] spectrum path and file name of a local file which provides the spectrum
  * @param[in] comm the working MPI communicator
  
  - The distribution of sparse matrix and vector is determined internally by this function. If you can cosider
  to use your own distribution scheme, please consider other two implementations of SMG2S for non-Hermitian matrix.
  - It can manage both the case with eigenvalues in real and the case with conjugate eigenvalues. The function will internally check whhich type
  of spectrum is provided in the local file
*/  
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

//! Generating a non-Symmetric sparse matrix using the spectrum stored in a parVector object
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] init a initMat object which determines the way of initialization of matrix, which will be further operated by SMG2S
  * @param[in] spec a parVector object which stores the given spectrum by users
  
  - The distribution of generated matrix across MPI procs is the same as the one of `spec` which provides the spectrum.
  - It works only with eigenvalues in real.
*/ 
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

//! Generating a non-Symmetric sparse matrix using the spectrum stored in a parVector object and user-provided initial matrix stored in a parMatrixSparse object
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] Am a parMatrixSparse provided by the users as a initial matrix of SMG2S, it can only be any sparse lower-triangular matrix
  * @param[in] spec a parVector object which stores the given spectrum by users
  * @param[out] Am for the output, it is overwritten by the generated matrix
  
  - The distribution of generated matrix across MPI procs is the same as the one of `spec` which provides the spectrum.
  - The distribution scheme of `Am` and `spec` should be the same.
  - It works only with eigenvalues in real.  
*/ 
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

//! Generating a non-Symmetric sparse matrix using the spectrum stored in a parVector object
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] init a initMat object which determines the way of initialization of matrix, which will be further operated by SMG2S
  * @param[in] spec a parVector object which stores the given spectrum by users
  
  - The distribution of generated matrix across MPI procs is the same as the one of `spec` which provides the spectrum.
  - It targets the case with conjugate eigenvalues.
*/ 
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

//! Generating a non-Symmetric sparse matrix using the spectrum stored in a parVector object and user-provided initial matrix stored in a parMatrixSparse object
/*!
  * @param[in] probSize row and column number of the sparse matrix to be generated
  * @param[in] nilp a Nilpotent object which determines the nilpotent matrix used by SMG2S
  * @param[in] Am a parMatrixSparse provided by the users as a initial matrix of SMG2S, it can only be any sparse lower-triangular matrix
  * @param[in] spec a parVector object which stores the given spectrum by users
  * @param[out] Am for the output, it is overwritten by the generated matrix
  
  - The distribution of generated matrix across MPI procs is the same as the one of `spec` which provides the spectrum.
  - The distribution scheme of `Am` and `spec` should be the same.
  - It targets the case with conjugate eigenvalues.
*/ 
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

/** @} */ // end of group2
#endif