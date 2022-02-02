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
#include "smg2s/smg2s_nonsymmetric.h"
#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include "utils/logo.h"
#include "utils/utils.h"
#include "utils/cmdline/cmdline.h"
#include <string>
#include <typeinfo>  

using S = int;

int main(int argc, char** argv) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int size;

    double start, end, time;

    bool non_sym = false;

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) {logo(1.1);}

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Command Line parser
    cmdline::parser parser;

    parser.add<S>("dim", 'D', "Dimension of matrix to be generated", false, 1000);
    parser.add<S>("lbandwidth", 'L', "low bandwidth of initial matrix", false, 5);
    parser.add<S>("continous", 'C', "Continuous length in Nilpotent matrix", false, 2);
    parser.add<std::string>("spectrum", 'S', "local file with given spectrum", false, " ");
    parser.add<std::string>("mattype", 'M', "Matrix type to be generated: non-symmetric or non-Hermitian", false, "non-herm", cmdline::oneof<std::string>("non-herm", "non-sym"));

    parser.parse_check(argc, argv);

    //parser value from cmdline
    S probSize = parser.get<S>("dim");
    S lbandwidth = parser.get<S>("lbandwidth");
    S length = parser.get<S>("continous");
    std::string spectrum = parser.get<std::string>("spectrum");
    std::string mattype = parser.get<std::string>("mattype");

    if(rank == 0) printf("INFO ]> Starting ... \n");

    if(rank == 0) printf("INFO ]> The MPI Comm World Size is %d\n", size);

    Nilpotency<S> nilp;
    nilp.NilpType1(length, probSize);

    if(mattype == "non-herm"){
        using T = std::complex<double>;

        parMatrixSparse<T, S> *Mt;

        start = MPI_Wtime();

        Mt =  smg2s<T, S>(probSize, nilp,lbandwidth, spectrum, MPI_COMM_WORLD);

        end = MPI_Wtime();

        time = end - start;

        if(rank == 0){
            border_print2();
            center_print ( "SMG2S Finish the Generation of Non Hermitian Matrix",100 );
            std::cout << "\n                              Size = "<< demical<int>(probSize) <<"e^" << pw<int>(probSize) << ", L = " << lbandwidth << ", C = " << length << ", Proc = " << size << "\n" << std::endl;
            std::cout <<  "                               Data Types for the test: " << typeid(T).name() <<", "<< typeid(S).name() <<  "\n" << std::endl;
            printf ( "                                  SMG2S Time is %f seconds \n", time );
            border_print2();
        }
    }
   else
   {

        using P = double;

        parMatrixSparse<P, S> *Mt2;

        start = MPI_Wtime();

        Mt2 =  smg2s_nonsymmetric<P, S>(probSize, nilp,lbandwidth, spectrum, MPI_COMM_WORLD);

        end = MPI_Wtime();

        time = end - start;

        if(rank == 0){
            border_print2();
            center_print ( "SMG2S Finish the Generation of Non Symmetric Matrix",100 );

            std::cout << "\n                              Size = "<< demical<int>(probSize) <<"e^" << pw<int>(probSize) << ", L = " << lbandwidth << ", C = " << length << ", Proc = " << size << "\n" << std::endl;
            std::cout <<  "                               Data Types for the test: " << typeid(P).name() <<", "<< typeid(S).name() <<  "\n" << std::endl;
            printf ( "                                  SMG2S Time is %f seconds \n", time );
            border_print2();
        }
   }

   MPI_Finalize();

    return 0;
}
