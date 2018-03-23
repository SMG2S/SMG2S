#ifndef __PARVECTORMAP_H__
#define __PARVECTORMAP_H__

#include <mpi.h>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

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
		S         global_size;

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

		S Loc2Glob(S local_index);
		S Glob2Loc(S global_index);

		//get
		int GetOwner(S index);
		int GetRank(){return rank;};

		S GetLowerBound(){return lower_bound;};
		S GetUpperBound(){return upper_bound;};
		S GetLocalSize(){printf("------> Proc %d: local size = %d \n\n", rank, local_size);return local_size;};
		S GetGlobalSize(){return global_size;};
		S GetLocTotSize(){return loctot_size;};

		//Active usrs

		int AddUser(){users++; return users;};
		int DeleteUser(){users--;return users;};
		int GetUser(){return users;};
};

#endif
