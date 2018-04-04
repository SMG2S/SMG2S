#include "../parMatrix/parMatrixSparse.cc"
#include <math.h>
#include <complex.h>
#include <string.h>

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

int main(int argc, char** argv){

	int printEnds = 1, printMids = 1;

	int i,j,k;

	int span, lower_b, upper_b;

	const int tag1 = 0, tag2 = 1;

	int rank, size, intBuf1, intBuf2, number;

	double doubBuf1, start, finish, time;

	double a = 5;

	MPI_Status status;


	int commMethod = 0;
	int probSize = 320000;

	int maxCount = 500;

//	        printf("HELLO\n");

/*
	if (argc>=2){
		for (int i =0; i < argc; i++){
			if (strcasecmp(argv[i],"-comm")){
				commMethod = atoi(argv[i+1]) ;
			}
                        if (strcasecmp(argv[i],"-size")==0){
                                probSize = atoi(argv[i+1]) ;
                        }
                        if (strcasecmp(argv[i],"-count")==0){
                                maxCount = atoi(argv[i+1]) ;
                        }

		}
	}
*/
//	printf("HELLO\n");
#ifndef _OPENMP

	MPI_Init(&argc,&argv) ;
	MPI_Comm_rank ( MPI_COMM_WORLD , &rank);
	MPI_Comm_size ( MPI_COMM_WORLD , &size);

	if(rank == 0)(std::cout << ">>>> pure MPI run \n" << std::endl);

#else
	int dummy; dummy = 0;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&dummy);

	MPI_Comm_rank ( MPI_COMM_WORLD , &rank);
        MPI_Comm_size ( MPI_COMM_WORLD , &size);

        if(rank == 0)(std::cout << ">>>> Hybrid MPI + OpenMP  run \n" << std::endl);

#endif

	span = int(floor(double(probSize)/double(size)));
	
	if (rank == size -1){

		lower_b = rank*span;
		upper_b = probSize;
	}else {

                lower_b = rank*span;
                upper_b = (rank+1)*span;
	}


	if(rank == 0){
		printf ( "------------------------------------\n" );
		printf ( "-------------  SPMV TEST -----------\n" );
		printf ( "-------%d Procs - %d x %d matrix----\n", size, probSize, probSize );
		printf ( "------------------------------------\n" );
	}	           

	MPI_Barrier(MPI_COMM_WORLD);

	parVector<double,int> *vec = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);
	vec->SetTovalue(a);

	if(rank == 0)(std::cout << ">>>> xVector done !!! \n" << std::endl);

	parVector<double,int> *prod = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);
	prod->SetTovalue(0.0);

	if(rank == 0)(std::cout << ">>>> yVector done !!! \n" << std::endl);

	parMatrixSparse<double,int> *Am = new parMatrixSparse<double,int>(vec,prod);	

	Am->ReadExtMat();

	//Am->MatView();

        if(rank == 0)(std::cout << ">>>> matrix done !!! \n" << std::endl);

	Am->ConvertToCSR();
	
        if(rank == 0)(std::cout << ">>>> matrix converted to CSR !!! \n" << std::endl);

	Am->FindColsToRecv();
	
        if(rank == 0)(std::cout << ">>>> matrix communication mapped !!! \n" << std::endl);

	Am->SetupDataTypes();

        if(rank == 0)(std::cout << ">>>> matrix datatype done !!! \n" << std::endl);

        MPI_Barrier(MPI_COMM_WORLD);

	start = MPI_Wtime();

	for (i=0;i<maxCount;i++){
		Am->CSR_MatVecProd(vec, prod);
	}

	finish  = MPI_Wtime();

	time = finish - start ;

	if(rank == 0){
		printf ( "------------------------------------\n" );
                printf ( "---- SPMV Time is %f seconds --------\n", time );
                printf ( "------------------------------------\n" );
	}

	MPI_Barrier(MPI_COMM_WORLD);

/*
	delete Am;
	delete vec;
	delete prod;
*/

	MPI_Finalize();

	return 0;
}

