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

#ifndef __PARVECTORMAP_H__
#define __PARVECTORMAP_H__

#include "utils/MPI_DataType.hpp"

//!  @brief A class which determines the way to distribute a vector across MPI procs.
/*!
   * @ingroup group6 
  - This class is to create a mapping from a fixed-size vector to multiple MPI procs in 1D grid.
  - This class can also be used to create more distributed vectors and sparse matrices following
  the same way.
  - For each MPI proc, a piece of vector with indexing `[lower_bound, upper_bound)` is allocated.

  @tparam S type of integer to describes the dimension of vector to be generated. 
*/
template<typename S>
class parVectorMap
{

    private:
    //! The working MPI Communicator
	MPI_Comm  	comm;
	//! number of MPI procs within the working MPI communicator parVectorMap#comm
	int	  	nproc;
	//! rank of each MPI procs within the working MPI communicator parVectorMap#comm
	int	  	rank;
	//! the smallest index of a distributed vector on each MPI proc
	S	  	lower_bound;
	//! `upper_bound-1 = `  the largest index of a distributed vector on each MPI proc 
	S	  	upper_bound;
	//! The number of elements of vector stored on each MPI proc
	S	  	local_size;
	//! Global size of this distributed vector
	S     	   	global_size;
	//! A `std::vector` which stores the parVectorMap#lower_bound of all MPI procs together
	std::vector<S>	lprocbound_map;
	//! A `std::vector` which stores the parVectorMap#upper_bound of all MPI procs together	
	std::vector<S> uprocbound_map;

    public:
	//constructor
	parVectorMap();
    //! A constructor of `parVectorMap`. 
    /*!
      * @param[in] ncomm the working MPI Communicator
      * @param[in] lbound the smallest index of a distributed vector on each MPI proc
      * @param[in] ubound `ubound-1 = `  the largest index of a distributed vector on each MPI proc 
    */		
	parVectorMap(MPI_Comm ncomm, S lbound, S ubound);
	//! A destructor of `parVectorMap`
	~parVectorMap();

	//! Compare if this map is identical to another one `map1`
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

	//! Compare if this map is different with another one `map1`
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

	//! Convert a index of local vector on each MPI proc into its index in the global distributed vector
    /*!
      * @param[in] local_index the index local vector to be converted

      - Attention, this function is in distributed manner, each MPI proc can only convert the local index of vector stored on itself.
    */		
	S Loc2Glob(S local_index);
	//! Convert a index of global vector into its index in the local vector on each MPI proc.
    /*!
      * @param[in] global_index the index global vector to be converted

      - Attention, this function is in distributed manner, each MPI proc can only convert the global index of vector in the range `[lower_bound, upper_bound)`.
    */		
	S Glob2Loc(S global_index);

	//get
	//! Return parVectorMap<S>#comm
	MPI_Comm GetCurrentComm(){return comm;};
	//! Return parVectorMap<S>#rank
	int GetRank(){return rank;};
	//! Return parVectorMap<S>#lower_bound
	S GetLowerBound(){return lower_bound;};
	//! Return parVectorMap<S>#upper_bound	
	S GetUpperBound(){return upper_bound;};
	//! Return parVectorMap<S>#local_size		
	S GetLocalSize(){return local_size;};
	//! Return parVectorMap<S>#global_size		
	S GetGlobalSize(){return global_size;};
	//! Return parVectorMap<S>#lprocbound_map	
	std::vector<S> GetLBoundMap(){return lprocbound_map;};
	//! Return parVectorMap<S>#uprocbound_map		
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

    }

    return local_index;
}


#endif