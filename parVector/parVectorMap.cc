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

   Part of basic data structures' implementation of this file refers to this technical report 
   (http://orbit.dtu.dk/files/51272329/tr12_10_Alexandersen_Lazarov_Dammann_1.pdf)
*/

#include "parVectorMap.h"

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



