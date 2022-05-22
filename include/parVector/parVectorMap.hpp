#ifndef __PARVECTORMAP_H__
#define __PARVECTORMAP_H__

#include "utils/MPI_DataType.hpp"

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
	std::vector<S>	lprocbound_map, uprocbound_map;

    public:
	//constructor
	parVectorMap();
	parVectorMap(MPI_Comm ncomm, S lbound, S ubound);
	//destructor
	~parVectorMap();

	bool operator == (const parVectorMap &map1){
	    bool ifsamecom;
	    int flag;
	    MPI_Comm_compare(comm, map1.comm, &flag);
	    if(flag == MPI_IDENT){
	    	ifsamecom = true;
	    }else{
	    	ifsamecom = false;
	    }

	    return (ifsamecom && nproc == map1.nproc && rank == map1.rank && lower_bound == map1.lower_bound && upper_bound == map1.upper_bound && local_size == map1.local_size && global_size == map1.global_size);
	};

	bool operator != (const parVectorMap &map1){
	    bool ifdiffcom;
	    int flag;
	    MPI_Comm_compare(comm, map1.comm, &flag);
	    if(flag == MPI_IDENT){
	    	ifdiffcom = false;
	    }else{
	    	ifdiffcom = true;
	    }

	    return (ifdiffcom || nproc != map1.nproc || rank != map1.rank || lower_bound != map1.lower_bound || upper_bound != map1.upper_bound || local_size != map1.local_size || global_size != map1.global_size);
	};

	MPI_Comm GetCurrentComm(){return comm;};

	S Loc2Glob(S local_index);
	S Glob2Loc(S global_index);

	//get
	int GetRank(){return rank;};
	S GetLowerBound(){return lower_bound;};
	S GetUpperBound(){return upper_bound;};
	S GetLocalSize(){return local_size;};
	S GetGlobalSize(){return global_size;};
	std::vector<S> GetLBoundMap(){return lprocbound_map;};
	std::vector<S> GetUBoundMap(){return uprocbound_map;};

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

    lprocbound_map.resize(nproc);
    uprocbound_map.resize(nproc);

    MPI_Allgather(&lower_bound, 1 , getMPI_Type<S>() , lprocbound_map.data() , 1 , getMPI_Type<S>(), comm ) ;
    MPI_Allgather(&upper_bound, 1 , getMPI_Type<S>() , uprocbound_map.data(), 1 ,  getMPI_Type<S>() , comm ) ;

}

template<typename S>
parVectorMap<S>::~parVectorMap(){
}

template<typename S>
S parVectorMap<S>::Loc2Glob(S local_index){

    S global_index = lower_bound + local_index;

    try{
    	if (global_index > global_size){
    	    global_index = -1;	
    	    throw (local_index);
    	}
    }
    catch(S local_index){
    	std::cout << "The given local index <" << local_index <<"> is out of bound." << std::endl;
    }
    
    return global_index;
}

template<typename S>
S parVectorMap<S>::Glob2Loc(S global_index){

    S local_index = global_index - lower_bound;

    try{
    	if (local_index < 0 || local_index >= local_size){
    	    local_index = -1;
    	    throw(global_index);	
    	}
    }catch(S global_index){
    	std::cout << "The given global index <" << global_index <<"> is out of bound." << std::endl;	
    }

    return local_index;
}


#endif