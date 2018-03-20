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
		// Allocate process boundary maps to number of processes
		lprocbound_map = new int[nproc] ; 
		uprocbound_map = new int [nproc] ;
		// Gather all lower bounds and all upper bounds
		MPI_Allgather(&lower_bound , 1 , MPI_INT , lprocbound_map , 1 , MPI_INT, comm ) ;
		MPI_Allgather(&upper_bound , 1 , MPI_INT , uprocbound_map , 1 , MPI_INT , comm ) ;
		// Find global size and set vectormap to describe upper - bounds

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

//convert local index to gobal or vice versa
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
	else return -1;
}

//get
template<typename S>
int parVectorMap<S>::GetOwner(S index)
{
	if((index < global_size) && (index >= 0)){
		for (int i = 0; i <nproc; i++){
			if(index < uprocbound_map[i]){return i;}
		}
	}
	else {return -1;}

}



