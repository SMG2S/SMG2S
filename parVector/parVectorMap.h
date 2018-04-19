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

#ifndef __PARVECTORMAP_H__
#define __PARVECTORMAP_H__

#include <mpi.h>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include "../config/config.h"


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

#endif
