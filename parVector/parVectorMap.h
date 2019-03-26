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
//#include "../config/config.h"


template<typename S>
class parVectorMap
{

	private:
		MPI_Comm  comm;
		int	  nproc;

		int	  rank;

		S	  lower_bound;
		S	  upper_bound;

		S	  local_size;
		S	  loctot_size;
		S     global_size;

		std::map<int,S> vectormap;
		int	  *lprocbound_map, *uprocbound_map;

		std::map<S,S> loc2glob;
		std::map<S,S> glob2loc;

		int users;

	public:
		//constructor
		parVectorMap(MPI_Comm ncomm, S lbound, S ubound);
		//destroyer
		~parVectorMap();

		MPI_Comm GetCurrentComm(){return comm;};

		S Loc2Glob(S local_index);
		S Glob2Loc(S global_index);

		//get
		int GetOwner(S index);
		int GetRank(){return rank;};

		S GetLowerBound(){return lower_bound;};
		S GetUpperBound(){return upper_bound;};
		S GetLocalSize(){return local_size;};
		S GetGlobalSize(){return global_size;};
		S GetLocTotSize(){return loctot_size;};


		//Active usrs

		int AddUser(){users++; return users;};
		int DeleteUser(){users--;return users;};
		int GetUser(){return users;};

};


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

	lprocbound_map = NULL;
	uprocbound_map = NULL;

	{
		lprocbound_map = new int[nproc] ;
		uprocbound_map = new int [nproc] ;
		// Gather all lower bounds and all upper bounds
		MPI_Allgather(&lower_bound , 1 , MPI_INT , lprocbound_map , 1 , MPI_INT, comm ) ;
		MPI_Allgather(&upper_bound , 1 , MPI_INT , uprocbound_map , 1 , MPI_INT , comm ) ;

		for (int i=0; i<nproc; i++) {
			vectormap [i] = uprocbound_map [i] ;
			if ( uprocbound_map [ i ] > global_size ) {
				global_size = uprocbound_map[i];
			}
	        }

	}

	loctot_size = local_size;
}

template<typename S>
parVectorMap<S>::~parVectorMap(){
	MPI_Comm_free(&comm);
	if(lprocbound_map!=NULL){
		delete [] lprocbound_map;
	}

        if(uprocbound_map!=NULL){
                delete [] uprocbound_map;
        }
}

template<typename S>
S parVectorMap<S>::Loc2Glob(S local_index){
	if((local_index < global_size) && (local_index >= 0))
		{return lower_bound + local_index;}
	else
		{return -1;}
}

template<typename S>
S parVectorMap<S>::Glob2Loc(S global_index){
	if((global_index >= lower_bound) && (global_index < upper_bound))
		{return global_index - lower_bound;}
	else
		{
			return -1;
		}
}

//get
template<typename S>
int parVectorMap<S>::GetOwner(S index)
{
	int m;
	if((index < global_size) && (index >= 0)){
		for (int i = 0; i <nproc; i++){
			if(index < uprocbound_map[i]){
				m = i;
				break;
			}
			else{
				m = -1;
			}
		}
	}
	else{
		m = -1;
	}
	return m;
}


#endif
