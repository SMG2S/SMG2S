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

//#include "parVectorMap.cc"
//#include "../parMatrix/parMatrixSparse.cc"
#include "smg2s/smg2s.h"
#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include "utils/logo.h"
#include "config/config.h"
#include "utils/utils.h"
#include <string>

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
        MPI_Finalize();
        return 0;
    }

    char *dim, *l, *c;

    std::string spectrum = " ";

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

        if (strcasecmp(argv[i],"-SPTR")==0){
                spectrum.assign(argv[i+1]);
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

    start = MPI_Wtime();

    Mt = smg2s<std::complex<double>,__int64_t>(probSize, nilp, lbandwidth,spectrum, MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start ;

    if(rank == 0){
            printf ( "----------------------------------------------------\n" );
            printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

            printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
            printf ( "----------------------------------------------------\n" );
    }

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

    Mt =  smg2s<std::complex<double>,int> (probSize, nilp,lbandwidth, spectrum,MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start ;

    if(rank == 0){
            printf ( "----------------------------------------------------\n" );
            printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

            printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
            printf ( "----------------------------------------------------\n" );
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

    start = MPI_Wtime();

    Mt = smg2s<std::complex<float>,__int64_t>(probSize, nilp, lbandwidth,spectrum,MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start ;

    if(rank == 0){
            printf ( "----------------------------------------------------\n" );
            printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

            printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
            printf ( "----------------------------------------------------\n" );
    }

#elif defined (__USE_COMPLEX__)
//complex single int

    int probSize, lbandwidth, length;

    probSize = atoi(dim);
    lbandwidth = atoi(l);
    length = atoi(c);

    Nilpotency<int> nilp;

    nilp.NilpType1(length,probSize);

    parMatrixSparse<std::complex<float>,int> *Mt;

    start = MPI_Wtime();

    Mt = smg2s<std::complex<float>,int>(probSize, nilp, lbandwidth,spectrum,MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start ;

    if(rank == 0){
            printf ( "----------------------------------------------------\n" );
            printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

            printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
            printf ( "----------------------------------------------------\n" );
    }


#elif defined (__USE_DOUBLE__) && defined(__USE_64BIT__)
//real double long int

    __int64_t probSize, lbandwidth, length;

    probSize = atoll(dim);
    lbandwidth = atoll(l);
    length = atoll(c);

    Nilpotency<__int64_t> nilp;

    nilp.NilpType1(length,probSize);

    parMatrixSparse<double,__int64_t> *Mt;

    start = MPI_Wtime();

    Mt = smg2s<double,__int64_t>(probSize, nilp, lbandwidth,spectrum,MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start ;

    if(rank == 0){
            printf ( "----------------------------------------------------\n" );
            printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

            printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
            printf ( "----------------------------------------------------\n" );
    }


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
    Mt = smg2s<double,int>(probSize, nilp, lbandwidth,spectrum,MPI_COMM_WORLD);
    end = MPI_Wtime();


    time = end - start ;

    if(rank == 0){

        	printf ( "----------------------------------------------------\n" );
	        printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

                printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
                printf ( "----------------------------------------------------\n" );
    }

//    Mt->LOC_MatView();

#elif defined (__USE_64BIT__)
//real single long int

    __int64_t probSize, lbandwidth, length;

    probSize = atoll(dim);
    lbandwidth = atoll(l);
    length = atoll(c);

    Nilpotency<__int64_t> nilp;

    nilp.NilpType1(length,probSize);

    parMatrixSparse<float,__int64_t> *Mt;

    start = MPI_Wtime();

    Mt = smg2s<float,__int64_t>(probSize, nilp, lbandwidth,spectrum,MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start ;

    if(rank == 0){
            printf ( "----------------------------------------------------\n" );
            printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

            printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
            printf ( "----------------------------------------------------\n" );
    }

#else
//real single int
    int probSize, lbandwidth, length;

    probSize = atoi(dim);
    lbandwidth = atoi(l);
    length = atoi(c);


    Nilpotency<int> nilp;

    nilp.NilpType1(length,probSize);

    parMatrixSparse<float,int> *Mt;

    start = MPI_Wtime();

    Mt = smg2s<float,int>(probSize, nilp, lbandwidth,spectrum,MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start ;
    
    if(rank == 0){
            printf ( "----------------------------------------------------\n" );
            printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
                            printf ( "----------------------------------------------------\n" );

            printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
            printf ( "----------------------------------------------------\n" );
    }


#endif

    //delete Mt;

    MPI_Finalize();

    return 0;
}

