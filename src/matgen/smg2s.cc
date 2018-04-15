
#include "smg2s.h"

template<typename T, typename S>
parMatrixSparse<T,S> *smg2s(S probSize, Nilpotency<S> nilp, S lbandwidth){

	int world_size;
	int world_rank;

	double start, end;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	 // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    S span, lower_b, upper_b;

    span = S(floor(double(probSize)/double(world_size)));

    printf("span = %d\n", span);


    if(world_rank == world_size - 1){
		lower_b = world_rank * span;
		upper_b = probSize - 1 + 1;
    }else{
		lower_b = world_rank * span;
		upper_b = (world_rank + 1) * span - 1 + 1;
    }

    printf("Proc. %d   Lower bound = %d   Upper bound = %d \n",world_rank, lower_b , upper_b );


	parVector<T,S> *vec = new parVector<T,S>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<T,S> *prod = new parVector<T,S>(MPI_COMM_WORLD, lower_b, upper_b);

    MPI_Barrier(MPI_COMM_WORLD);

    //generate vec containing the given spectra
    specGen<T,S>(vec);

    //vec->VecView();

    //Matrix Initialization

    parMatrixSparse<T,S> *Am = new parMatrixSparse<T,S>(vec,prod);
    parMatrixSparse<T,S> *MA = new parMatrixSparse<T,S>(vec,prod);
    parMatrixSparse<T,S> *AM = new parMatrixSparse<T,S>(vec,prod);

    parMatrixSparse<T,S> *matAop = new parMatrixSparse<T,S>(vec,prod);

    if(world_rank == 0){printf("Matrix Initialized\n");}

    MPI_Barrier(MPI_COMM_WORLD);

    //setup the lower part of initial matrix

    start = MPI_Wtime();

    T rnd;

    for(S i = 0; i < probSize; i++){
        for(S j = i - lbandwidth; j < i; j++){
            if(j >= 0){
            	rnd = 0.5*random<T,S>(0,10);
            	//if(world_rank == 0) printf("rnd = %f\n", rnd);
                Am->Loc_SetValue(i,j,rnd);   
                matAop->Loc_SetValue(i,j,rnd);   

            }
        }
    }

    //insert the diagonal of initial matrix with given spectra.

    Am->Loc_SetDiagonal(vec);
    matAop->Loc_SetDiagonal(vec);

    end = MPI_Wtime();

    double t2 = end - start;

    if(world_rank == 0) {printf("Initial matrix generation time = %1.6f\n", t2);}

    MPI_Barrier(MPI_COMM_WORLD);

    S my_factorielle_bornes = 1;

    my_factorielle_bornes = factorial(1,(2*nilp.nbOne-2));

    printf("my_factorielle_bornes = %d\n", my_factorielle_bornes);
    
    Am->Loc_MatScale((T)my_factorielle_bornes);

    for (S k=1; k<=(2*nilp.nbOne-2); k++){
    	matAop->MA(nilp, MA);
    	matAop->AM(nilp, AM);
    	matAop->Loc_MatAYPX(AM, 0);
    	matAop->Loc_MatAXPY(MA, -1);

    	my_factorielle_bornes =  factorial(k+1,2*nilp.nbOne-2);
    	Am->Loc_MatAXPY(matAop, (T)my_factorielle_bornes);

    	MA->Loc_ZeroEntries();
    	AM->Loc_ZeroEntries();
    }
    
    Am->LOC_MatView();

	my_factorielle_bornes = factorial(1,(2*nilp.nbOne-2));

	//T inv = 1/(T)my_factorielle_bornes;
#ifdef __USE_DOUBLE__
    double fac = (double)my_factorielle_bornes;
    double inv = 1/fac;
#else
    float fac = (float)my_factorielle_bornes;
    float inv = 1/fac;
#endif

	Am->Loc_MatScale(inv);

    Am->LOC_MatView();

    return Am;
}
