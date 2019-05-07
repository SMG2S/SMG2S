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

#ifndef __SMG2S_NONSYMMETRIC_H__
#define __SMG2S_NONSYMMETRIC_H__

#include "../parVector/parVector.h"
#include "../parMatrix/parMatrixSparse.h"
#include "specGen_nonsymmetric.h"
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
parMatrixSparse<T,S> *smg2s_nonsymmetric(S probSize, Nilpotency<S> nilp, S lbandwidth, std::string spectrum, MPI_Comm comm){


  int world_size;
	int world_rank;

	double start, end;

	MPI_Comm_size(comm, &world_size);

	 // Get the rank of the process
    MPI_Comm_rank(comm, &world_rank);

    if(std::is_same<T,std::complex<double> >::value || std::is_same<T,std::complex<float> >::value){
        if(world_rank == 0){printf("ERROR ]: For the nonsymmetric case, the scalar type of matrix cannot be complex\n");}
        return 0;
    }

    MPI_Barrier(comm);

    S span, lower_b, upper_b;

    span = S(ceil(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
		lower_b = world_rank * span;
		upper_b = probSize - 1 + 1;
    }else{
		lower_b = world_rank * span;
		upper_b = (world_rank + 1) * span - 1 + 1;
    }


    parVector<std::complex<T>,S> *spec = new parVector<std::complex<T>,S>(comm, lower_b, upper_b);
    parVector<T,S> *vec = new parVector<T,S>(comm, lower_b, upper_b);


    MPI_Barrier(comm);

    //generate vec containing the given spectra
    spec->specGen2(spectrum);
    
    //spec->VecView();

    //Matrix Initialization

    parMatrixSparse<T,S> *Am = new parMatrixSparse<T,S>(vec,vec);
    parMatrixSparse<T,S> *MA = new parMatrixSparse<T,S>(vec,vec);
    parMatrixSparse<T,S> *AM = new parMatrixSparse<T,S>(vec,vec);
    parMatrixSparse<T,S> *matAop = new parMatrixSparse<T,S>(vec,vec);

    MPI_Barrier(comm);

    //setup the lower part of initial matrix

    start = MPI_Wtime();


    matInit2(Am, matAop, probSize, lbandwidth, spec);

    //Am->LOC_MatView();

    end = MPI_Wtime();

    double t2 = end - start;

    if(world_rank == 0) {printf("Initial matrix generation time = %1.6f\n", t2);}

    MPI_Barrier(comm);

    
    __int64_t my_factorielle_bornes = 1;

    my_factorielle_bornes = factorial(1,2*nilp.nbOne);

    Am->Loc_MatScale((double)my_factorielle_bornes);


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
