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
 */

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
