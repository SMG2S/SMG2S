#ifndef __MPI_DTATYPE_H__
#define __MPI_DTATYPE_H__

//#include <random>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <mpi.h>
#include <complex>

template<class T>
MPI_Datatype MPI_Scalar(){
	
	int rank, size;
	MPI_Datatype MPI_Scalar;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(std::is_same<T,std::complex<double> >::value || std::is_same<T,std::complex<float> >::value){
		if(rank == 0){
			std::cout << "INFO ]>: Using Complex values" << std::endl;
		}
	} else if(std::is_same<T,double>::value || std::is_same<T,float>::value){
		if(rank == 0){
			std::cout << "INFO ]>: Using Real values" << std::endl; 
		}
	}

	if(std::is_same<T,std::complex<double> >::value){
		MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_Scalar);
	} else if(std::is_same<T,std::complex<float> >::value){
		MPI_Type_contiguous(2, MPI_FLOAT, &MPI_Scalar);
	} else if (std::is_same<T,double>::value){
		MPI_Type_contiguous(1, MPI_DOUBLE, &MPI_Scalar);
	} else if (std::is_same<T,float>::value){
		MPI_Type_contiguous(1, MPI_FLOAT, &MPI_Scalar);
	}

	MPI_Type_commit(&MPI_Scalar);

	return MPI_Scalar;

}; 


#endif
