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

#ifndef __MPI_DTATYPE_H__
#define __MPI_DTATYPE_H__

//#include <random>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <mpi.h>
#include <complex>
#include <stdint.h>

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

template <typename T>
MPI_Datatype getMPI_Type();

template <>
MPI_Datatype getMPI_Type<float>() {
  return MPI_FLOAT;
}

template <>
MPI_Datatype getMPI_Type<double>() {
  return MPI_DOUBLE;
}

template <>
MPI_Datatype getMPI_Type<std::complex<float> >() {
  return MPI_COMPLEX;
}

template <>
MPI_Datatype getMPI_Type<std::complex<double> >() {
  return MPI_DOUBLE_COMPLEX;
}

template <>
MPI_Datatype getMPI_Type<int>() {
  return MPI_INT;
}

template <>
MPI_Datatype getMPI_Type<__int64_t>() {
  return MPI_INT64_T;
}

template <>
MPI_Datatype getMPI_Type<long>() {
  return MPI_LONG;
}



#endif
