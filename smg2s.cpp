/*
MIT License
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France
Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
                                     Materials,  Forschungszentrum Juelich GmbH.
                                     
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

#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include <string>
#include <typeinfo>  
#include <random>
#include <cmath>
#include <smg2s-interface.hpp>
#include <utils/logo.hpp>
#include <utils/cmdline/cmdline.hpp>

const double pi = std::acos(-1);

using S = int;
using T1 = std::complex<double>;
using T2 = double;

void specGenNonHerm(parVector<T1, S> *spec){
    std::mt19937_64 rd(12234);
    std::uniform_real_distribution<> d(0, 1);

    auto lb = spec->GetLowerBound();
    auto ub = spec->GetUpperBound();
    auto n = spec->GetGlobalSize();

    for(auto i = lb; i < ub; i++){
        T1 v(std::cos(i * 2 * pi / n) * 100 + 100 + 50 * d(rd), 100 * std::sin(i * 2 * pi / n) );
        spec->SetValueGlobal(i, v);
    }
}

void specGenNonSymmConj(parVector<std::complex<T2>, S> *spec){
    std::mt19937_64 rd(12234);
    std::uniform_real_distribution<> d(0, 1);

    auto lb = spec->GetLowerBound();
    auto ub = spec->GetUpperBound();
    auto n = spec->GetGlobalSize();
    
    for(auto i = 0; i < n; i = i + 2){
        if(i >= lb && i < ub){
            std::complex<T2> v(std::cos(i * pi / n) * 1000 + 1000 , 1000 * std::sin(i * pi / n)+0.001 );
            spec->SetValueGlobal(i, v);
            std::complex<T2> v2(std::cos(i * pi / n) * 1000 + 1000 , -1000 * std::sin(i * pi / n)-0.001 );
            spec->SetValueGlobal(i+1, v2);
        }
    }

}


int main(int argc, char** argv) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    int nProcs;
    int MyPID;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);

    if(MyPID == 0) {logo(1.1);}

    // Command Line parser
    cmdline::parser parser;

    parser.add<S>("dim", 'D', "Dimension of matrix to be generated", false, 1000);
    parser.add<S>("diagL", 'L', "offset of lower diagonal of initial matrix", false, -10);
    parser.add<S>("diagU", 'U', "offset of upper diagonal of initial matrix", false, -5);
    parser.add<S>("nilpOffset", 'O', "offset of diagonal of a nilpotent", false, 5);
    parser.add<S>("continous", 'C', "Continuous length in Nilpotent matrix", false, 2);

    parser.add<std::string>("mattype", 'M', "Matrix type to be generated: non-symmetric or non-Hermitian", false, "non-herm", cmdline::oneof<std::string>("non-herm", "non-symm"));

    parser.parse_check(argc, argv);

    //parser value from cmdline
    S probSize = parser.get<S>("dim");
    S diag_l = parser.get<S>("diagL");
    S diag_u = parser.get<S>("diagU");
    S length = parser.get<S>("continous");
    S offset = parser.get<S>("nilpOffset");
    std::string mattype = parser.get<std::string>("mattype");

    S span, lower_b, upper_b;

    span = S(floor(double(probSize)/double(nProcs)));

    if(MyPID == nProcs - 1){
        lower_b = MyPID * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = MyPID * span;
        upper_b = (MyPID + 1) * span - 1 + 1;
    }

    auto parVecMap = parVectorMap<S>(MPI_COMM_WORLD, lower_b, upper_b);
    
    if(MyPID == 0) printf("INFO ]> Starting ... \n");

    if(MyPID == 0) printf("INFO ]> The MPI Comm World Size is %d\n", nProcs);

    Nilpotent<S> nilp = Nilpotent<S>(length, offset, probSize);

    double start, end, time;

    if(mattype == "non-herm"){
        auto spec1 = parVector<T1, S>(parVecMap);

        specGenNonHerm(&spec1);

        start = MPI_Wtime();
        auto mat = nonherm<T1, S>(probSize, nilp, initMat<S>(diag_l, diag_u, 0.1, 0.95), spec1);
        end = MPI_Wtime();

        time = end - start;

        if(MyPID == 0){
            border_print2();
            center_print ( "SMG2S Finish the Generation of Non Hermitian Matrix",100 );
            std::cout << "\n                              Size = "<< demical<int>(probSize) <<"e^" << pw<int>(probSize) << ", L = " << diag_l << ", U = " << diag_u << ", O = " << offset << ", C = " << length << ", Proc = " << nProcs << "\n" << std::endl;
            std::cout <<  "                               Data Types for the test: " << typeid(T1).name() <<", "<< typeid(S).name() <<  "\n" << std::endl;
            printf ( "                                  SMG2S Time is %f seconds \n", time );
            border_print2();
        }
    }
   else
   {

        auto spec2 = parVector<std::complex<T2>, S>(parVecMap);

        specGenNonSymmConj(&spec2);

        start = MPI_Wtime();
        auto mat2 = nonsymmconj<T2, S>(probSize, nilp, initMat<S>(diag_l, diag_u, 0.1, 0.95), spec2);
        end = MPI_Wtime();

        time = end - start;
        
        if(MyPID == 0){
            border_print2();
            center_print ( "SMG2S Finish the Generation of Non Symmetric Matrix",100 );
            std::cout << "\n                              Size = "<< demical<int>(probSize) <<"e^" << pw<int>(probSize) << ", L = " << diag_l << ", U = " << diag_u << ", O = " << offset << ", C = " << length << ", Proc = " << nProcs << "\n" << std::endl;
            std::cout <<  "                               Data Types for the test: " << typeid(T2).name() <<", "<< typeid(S).name() <<  "\n" << std::endl;
            printf ( "                                  SMG2S Time is %f seconds \n", time );
            border_print2();
        }
   }

   MPI_Finalize();

}
