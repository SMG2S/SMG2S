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

#ifndef __PARVECTORMAP_H__
#define __PARVECTORMAP_H__

#include <mpi.h>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>

template<typename S>
class parVectorMap
{

    private:
	MPI_Comm  	comm;
	int	  	nproc;
	int	  	rank;
	S	  	lower_bound;
	S	  	upper_bound;
	S	  	local_size;
	S     	   	global_size;

    public:
	//constructor
	parVectorMap();
	parVectorMap(MPI_Comm ncomm, S lbound, S ubound);
	//destructor
	~parVectorMap();

	MPI_Comm GetCurrentComm(){return comm;};

	S Loc2Glob(S local_index);
	S Glob2Loc(S global_index);

	//get
	int GetRank(){return rank;};
	S GetLowerBound(){return lower_bound;};
	S GetUpperBound(){return upper_bound;};
	S GetLocalSize(){return local_size;};
	S GetGlobalSize(){return global_size;};
};

template<typename S>
parVectorMap<S>::parVectorMap(){}

template<typename S>
parVectorMap<S>::parVectorMap(MPI_Comm ncomm, S lbound, S ubound)
{
    MPI_Comm_dup(ncomm, &comm);
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);

    //setting lower and upper bound
    lower_bound = lbound;
    upper_bound = ubound;
    //get local vector size
    local_size = upper_bound - lower_bound;
    //initial global size
    global_size = 0;

    MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, comm);	
}

template<typename S>
parVectorMap<S>::~parVectorMap(){
}

template<typename S>
S parVectorMap<S>::Loc2Glob(S local_index){

    S global_index = lower_bound + local_index;

    if (global_index > global_size)
    {
    	throw "index out of range";
    }else{
    	return global_index;	
    }
    
}

template<typename S>
S parVectorMap<S>::Glob2Loc(S global_index){

    S local_index = global_index - lower_bound;
    
    if (local_index < 0 || local_index >= local_size)
    {
    	throw "index out of range";
    }else{
    	return local_index;	
    }
}


#endif
