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
    //std::string mattype = "non-herm";
    std::string mattype = " ";
    //non-sym

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

        if (strcasecmp(argv[i],"-mattype")==0){
                mattype.assign(argv[i+1]);
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

    __int64_t probSize, lbandwidth, length;

    probSize = atoi(dim);
    lbandwidth = atoi(l);
    length = atoi(c);

    Nilpotency<__int64_t> nilp;

    nilp.NilpType2(length,probSize);

/*Non symmetric case*/

    parMatrixSparse<float,__int64_t> *Mt2;

    start = MPI_Wtime();

    Mt2 =  smg2s<float,__int64_t>(probSize, nilp,lbandwidth, spectrum, mattype, MPI_COMM_WORLD);

    end = MPI_Wtime();

    time = end - start;

    if(rank == 0){
        printf ( "----------------------------------------------------\n" );
        printf ( "----- SMG2S Finish the Matrix Generation------------\n" );
        printf ( "--------------- Non symmetric Matrix ---------------\n" );
        printf ( "----------------------------------------------------\n" );
        printf ( "---------- SMG2S Time is %f seconds ----------\n", time );
        printf ( "----------------------------------------------------\n" );
    }



    //delete Mt;

    MPI_Finalize();

    return 0;
}
