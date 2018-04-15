//#include "parVectorMap.cc"
//#include "../parMatrix/parMatrixSparse.cc"
#include "src/matgen/smg2s.cc"
#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include "utils/logo.h"
#include "config/config.h"
#include "utils/utils.h"

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

int main(int argc, char** argv) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int size;

    double start, end, time;

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) {logo(1.0);}

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    //MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0) printf("INFO ]> Starting ... \n");
    
    // Print off a hello world message
    if(rank == 0) printf("INFO ]> The MPI Comm World Size is %d\n", size);

    if (argc < 5) {
        // Tell the user how to run the program
        if(rank == 0){
            show_usage(argv[0]);
        }
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        MPI_Finalize();
        return 0;
    }

    char *dim, *l, *c;

    for (int i =0; i < argc; i++){

        if (strcasecmp(argv[i],"-SIZE")==0){
                dim = argv[i+1] ;
        }
        if (strcasecmp(argv[i],"-L")==0){
                l = argv[i+1] ;
        }

        if (strcasecmp(argv[i],"-C")==0){
                c = argv[i+1] ;
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);



#if defined(__USE_COMPLEX__) && defined(__USE_DOUBLE__) && defined (__USE_64BIT__)
//complex double long int

    __int64_t probSize, lbandwidth, length;

    probSize = atoll(dim);
    lbandwidth = atoll(l);
    length = atoll(c);

    Nilpotency<__int64_t> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<std::complex<double>,__int64_t> *Mt;

    Mt = smg2s<std::complex<double>,__int64_t>(probSize, nilp, lbandwidth);


#elif defined (__USE_COMPLEX__) && defined(__USE_DOUBLE__)
//complex double int

    int probSize, lbandwidth, length;

    probSize = atoi(dim);
    lbandwidth = atoi(l);
    length = atoi(c);

    Nilpotency<int> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<std::complex<double>,int> *Mt;

//    smg2s<std::complex<double>,int> (probSize, nilp, lbandwidth);
    start = MPI_Wtime();
    
    Mt =  smg2s<std::complex<double>,int> (probSize, nilp,lbandwidth);

    end = MPI_Wtime();

    time = end - start ;

    if(rank == 0){
        printf ( "------------------------------------\n" );
                printf ( "---- SMG2S Time is %f seconds --------\n", time );
                printf ( "------------------------------------\n" );
    }

#elif defined (__USE_COMPLEX__) && defined(__USE_64BIT__)
//complex  single long int

    __int64_t probSize, lbandwidth, length;

    probSize = atoll(dim);
    lbandwidth = atoll(l);
    length = atoll(c);

    Nilpotency<__int64_t> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<std::complex<float>,__int64_t> *Mt;

    Mt = smg2s<std::complex<float>,__int64_t>(probSize, nilp, lbandwidth);

#elif defined (__USE_COMPLEX__)
//complex single int

    int probSize, lbandwidth, length;

    probSize = atoi(dim);
    lbandwidth = atoi(l);
    length = atoi(c);

    Nilpotency<int> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<std::complex<float>,int> *Mt;

    Mt = smg2s<std::complex<float>,int>(probSize, nilp, lbandwidth);

#elif defined (__USE_DOUBLE__) && defined(__USE_64BIT__)
//real double long int

    __int64_t probSize, lbandwidth, length;

    probSize = atoll(dim);
    lbandwidth = atoll(l);
    length = atoll(c);

    Nilpotency<__int64_t> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<double,__int64_t> *Mt;

    Mt = smg2s<double,__int64_t>(probSize, nilp, lbandwidth);

#elif defined (__USE_DOUBLE__)
//real double int

    int probSize, lbandwidth, length;

    probSize = atoi(dim);
    lbandwidth = atoi(l);
    length = atoi(c);

    Nilpotency<int> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<double,int> *Mt;

    start = MPI_Wtime();
    Mt = smg2s<double,int>(probSize, nilp, lbandwidth);
    end = MPI_Wtime();


    time = end - start ;

    if(rank == 0){
        printf ( "------------------------------------\n" );
                printf ( "---- SMG2S Time is %f seconds --------\n", time );
                printf ( "------------------------------------\n" );
    }

   // Mt->LOC_MatView();

#elif defined (__USE_64BIT__)
//real single long int

    __int64_t probSize, lbandwidth, length;

    probSize = atoll(dim);
    lbandwidth = atoll(l);
    length = atoll(c);

    Nilpotency<__int64_t> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<float,__int64_t> *Mt;

    Mt = smg2s<float,__int64_t>(probSize, nilp, lbandwidth);
else
//real single int
    int probSize, lbandwidth, length;

    probSize = atoi(dim);
    lbandwidth = atoi(l);
    length = atoi(c);


    Nilpotency<int> nilp;
    
    nilp.NilpType1(length,probSize);

    parMatrixSparse<float,int> *Mt;


    Mt = smg2s<float,int>(probSize, nilp, lbandwidth);

#endif

    MPI_Finalize();

    return 0;
}

