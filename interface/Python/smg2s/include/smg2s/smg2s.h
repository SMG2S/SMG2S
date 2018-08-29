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

#ifndef __SMG2S_H__
#define __SMG2S_H__

#include "../parVector/parVector.h"
#include "../parMatrix/parMatrixSparse.h"
#include "specGen.h"
#include <math.h>
#include <complex.h>
#include <string>

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

template<typename T, typename S>
parMatrixSparse<T,S> *smg2s(S probSize, Nilpotency<S> nilp, S lbandwidth, std::string spectrum, MPI_Comm comm){

	int world_size;
	int world_rank;

	double start, end;

	MPI_Comm_size(comm, &world_size);

	 // Get the rank of the process
    MPI_Comm_rank(comm, &world_rank);

    S span, lower_b, upper_b;

    span = S(ceil(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
		lower_b = world_rank * span;
		upper_b = probSize - 1 + 1;
    }else{
		lower_b = world_rank * span;
		upper_b = (world_rank + 1) * span - 1 + 1;
    }


	  parVector<T,S> *vec = new parVector<T,S>(comm, lower_b, upper_b);
    parVector<T,S> *prod = new parVector<T,S>(comm, lower_b, upper_b);

    MPI_Barrier(comm);

    //generate vec containing the given spectra
    vec->specGen(spectrum);

    //Matrix Initialization

    parMatrixSparse<T,S> *Am = new parMatrixSparse<T,S>(vec,prod);
    parMatrixSparse<T,S> *MA = new parMatrixSparse<T,S>(vec,prod);
    parMatrixSparse<T,S> *AM = new parMatrixSparse<T,S>(vec,prod);

    parMatrixSparse<T,S> *matAop = new parMatrixSparse<T,S>(vec,prod);

    MPI_Barrier(comm);

    //setup the lower part of initial matrix

    start = MPI_Wtime();

    matInit(Am, matAop, probSize, lbandwidth);

    Am->Loc_SetDiagonal(vec);
    matAop->Loc_SetDiagonal(vec);

    end = MPI_Wtime();

    double t2 = end - start;

    if(world_rank == 0) {printf("Initial matrix generation time = %1.6f\n", t2);}

    MPI_Barrier(comm);

    __int64_t my_factorielle_bornes = 1;

    my_factorielle_bornes = factorial(1,2*nilp.nbOne);

    Am->Loc_MatScale((T)my_factorielle_bornes);

    for (S k=1; k<=2*nilp.nbOne; k++){

    	matAop->MA(nilp, MA);
  	  matAop->AM(nilp, AM);
    	matAop->Loc_MatAYPX(AM, 0);
    	matAop->Loc_MatAXPY(MA, -1);

  	  my_factorielle_bornes = factorial(k+1,2*nilp.nbOne);
    	Am->Loc_MatAXPY(matAop, (double)my_factorielle_bornes);
    	MA->Loc_ZeroEntries();
    	AM->Loc_ZeroEntries();

    }


	my_factorielle_bornes = factorial(1,2*(nilp.nbOne));

  double fac = (double)my_factorielle_bornes;
  double inv = 1/fac;

	Am->Loc_MatScale((T)inv);

  //Am->LOC_MatView();

    return Am;
}

#endif
