#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include <string>
#include <typeinfo>  
#include <mpi.h>
#include <random>
#include <cmath>

#include <smg2s-interface.hpp>

const double pi = std::acos(-1);

void specGenNonHerm1(parVector<std::complex<double>, int> *spec){
    std::mt19937_64 rd(12234);
    std::uniform_real_distribution<> d(0, 1);

    auto lb = spec->GetLowerBound();
    auto ub = spec->GetUpperBound();
    auto n = spec->GetGlobalSize();

    for(auto i = lb; i < ub; i++){
        std::complex<double> v(std::cos(i * 2 * pi / n) * 1 + 1 + d(rd), 1 * std::sin(i * 2 * pi / n) );
        spec->SetValueGlobal(i, v);
    }
}

//four clustered eigenvalues
void specGenNonHerm2(parVector<std::complex<double>, int> *spec){
    std::mt19937_64 rd(12234);
    std::uniform_real_distribution<> d(0, 1);

    auto lb = spec->GetLowerBound();
    auto ub = spec->GetUpperBound();
    auto n = spec->GetGlobalSize();

    for(auto i = lb; i < ub; i++){
        if(i < 0.25 * n){
            std::complex<double> v(std::cos(i * 8 * pi / n) * 1 + 1, 1 * std::sin(i * 8 * pi / n)  + 1);
            spec->SetValueGlobal(i, v);
        }else if(i >= 0.25 * n && i < 0.5 * n){
            std::complex<double> v(std::cos(i * 8 * pi / n) * 100 + 100, 100 * std::sin(i * 8 * pi / n) - 150);
            spec->SetValueGlobal(i, v);
        }else if(i >= 0.5 * n && i < 0.75 * n){
            std::complex<double> v(std::cos(i * 8 * pi / n) * 1000 + 1000, 1000 * std::sin(i * 8 * pi / n) + 50000);
            spec->SetValueGlobal(i, v);
        }else{
            std::complex<double> v(std::cos(i * 8 * pi / n) * 10000 + 10000, 10000 * std::sin(i * 8 * pi / n) - 20000);
            spec->SetValueGlobal(i, v);
        }

    }
}


void specGenNonSymm(parVector<double, int> *spec, double scale, double shift){
    std::mt19937_64 rd(12234);
    std::uniform_real_distribution<> d(0, 1);

    auto lb = spec->GetLowerBound();
    auto ub = spec->GetUpperBound();
    auto n = spec->GetGlobalSize();

    for(auto i = lb; i < ub; i++){
        spec->SetValueGlobal(i, scale * d(rd) + shift);
    }
}


void specGenNonSymm2(parVector<std::complex<double>, int> *spec, double scale, double shift){
    std::mt19937_64 rd(12234);
    std::uniform_real_distribution<> d(0, 1);

    auto lb = spec->GetLowerBound();
    auto ub = spec->GetUpperBound();
    auto n = spec->GetGlobalSize();

    for(auto i = lb; i < ub; i++){
        std::complex<double> v( scale * d(rd) + shift, 0);        
        spec->SetValueGlobal(i, v);
    }
}


void specGenNonSymmConj(parVector<std::complex<double>, int> *spec){
    std::mt19937_64 rd(12234);
    std::uniform_real_distribution<> d(0, 1);

    auto lb = spec->GetLowerBound();
    auto ub = spec->GetUpperBound();
    auto n = spec->GetGlobalSize();
    
    try{
        if(n % 2 != 0){
            throw 505;
        }
    }catch(...){
        std::cout << "for this generator of spectrum, the size of matrix needs to be a even number" << std::endl;
    }
    
    
    for(auto i = 0; i < n; i = i + 2){
        if(i >= lb && i < ub){
            std::complex<double> v(std::cos(i * 1 * pi / n) * 1000 + 1000 , 1000 * std::sin(i * 1 * pi / n)+0.001 );
            spec->SetValueGlobal(i, v);
            std::complex<double> v2(std::cos(i * 1 * pi / n) * 1000 + 1000 , -1000 * std::sin(i * 1 * pi / n)-0.001 );
            spec->SetValueGlobal(i+1, v2);
        }
    }

}


int main(int argc, char** argv) 
{
    MPI_Init(&argc, &argv);
    
    int world_size;
    int world_rank;
    int probSize = 100;
    int l_diag = -20;
    int u_diag = -5;
    int nbOne = 10;
    int offset = 20;
    double sparsity = 0.9;

    Nilpotent<int> nilp = Nilpotent<int>(nbOne, offset, probSize);

    nilp.show();

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int span, lower_b, upper_b;

    span = int(floor(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    auto parVecMap = parVectorMap<int>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<std::complex<double>, int> spec1 = parVector<std::complex<double>, int>(parVecMap);
    specGenNonHerm2(&spec1);
    spec1.writeToTxtCmplx("specNonHerm1.txt");

    auto mat = nonherm<std::complex<double>, int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 0.1, sparsity), spec1);
    mat.writeToMatrixMarketCmplx("matGenNonHerm1.mtx");

    parVector<double, int> spec2 = parVector<double, int>(parVecMap);
    specGenNonSymm(&spec2, 100000.0, 0.0);
    spec2.writeToTxt("specNonSymm1.txt");
    auto mat2 = nonsymm<double , int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 0.1, sparsity), spec2);
    mat2.writeToMatrixMarket("matGenNonSymm1.mtx");

    parVector<std::complex<double>, int> spec3 = parVector<std::complex<double>, int>(parVecMap);
    specGenNonSymmConj(&spec3);
    spec3.writeToTxtCmplx("specNonSymmConj1.txt");    
//    spec3.VecView();
    auto mat3 = nonsymmconj<double , int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 0.1, sparsity), spec3);
    mat3.writeToMatrixMarket("matGenNonSymmConj1.mtx");

    parVector<std::complex<double>, int> spec4 = parVector<std::complex<double>, int>(parVecMap);
    specGenNonSymm2(&spec4, 100000.0, 0.0);
    spec4.writeToTxtCmplx("specNonSymm2.txt");    
//    spec3.VecView();
    auto mat4 = nonsymmconj<double , int>(probSize, nilp, initMat<int>(-l_diag, -u_diag, 0.1, sparsity), spec4);
    mat4.writeToMatrixMarket("matGenNonSymm2.mtx");


	MPI_Finalize();
    

	return 0;
}
