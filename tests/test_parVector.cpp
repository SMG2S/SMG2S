#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include <string>
#include <typeinfo>  
#include <mpi.h>

#include "../parVector/parVectorMap.h"
#include "../parVector/parVector.h"

int main(int argc, char** argv) 
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    int world_size;
    int world_rank;
    int probSize = 7;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int span, lower_b, upper_b;

    span = int(ceil(double(probSize)/double(world_size)));

    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    auto parVecMap = parVectorMap<int>(MPI_COMM_WORLD, lower_b, upper_b);

    std::cout << parVecMap.GetRank() << " GetLowerBound: " << parVecMap.GetLowerBound() << std::endl;
    std::cout << parVecMap.GetRank() << " GetUpperBound: " << parVecMap.GetUpperBound() << std::endl;
    std::cout << parVecMap.GetRank() << " GetGlobalSize: " << parVecMap.GetGlobalSize() << std::endl;
 	std::cout << parVecMap.GetRank() << " GetLocalSize: " << parVecMap.GetLocalSize() << std::endl;

    parVector<int, int> vec = parVector<int, int>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<int, int> vec2 = parVector<int, int>(parVecMap);

    std::cout << world_rank << " GetLowerBound: " << vec.GetLowerBound() << std::endl;
    std::cout << world_rank << " GetUpperBound: " << vec.GetUpperBound() << std::endl;
    std::cout << world_rank << " GetGlobalSize: " << vec.GetGlobalSize() << std::endl;
 	std::cout << world_rank << " GetLocalSize: " << vec.GetLocalSize() << std::endl;

    for(auto i = lower_b; i < upper_b; i++){
    	std::cout <<  vec.GetRank() << " Glob2Loc " << i << ": " << vec.Glob2Loc(i) << std::endl;
    }

    for(auto i = 0; i < parVecMap.GetLocalSize(); i++){
    	std::cout <<  vec.GetRank() << " Loc2Glob " << i << ": " << vec.Loc2Glob(i) << std::endl;
    }
    vec.SetToValue(1);
 	vec2.SetToValue(4);
 	int *vals = new int[5];
 	int *rows = new int[5];

 	for(auto i = 0; i < 5; i++){
 		//vec.SetValueLocal(i, i+1);
 		vals[i] = i + 2;
 		rows[i] = i;
 	}

 	vec.SetValuesLocal(5, rows, vals);

	for(auto i = lower_b; i < upper_b; i++){
		vec.SetValueGlobal(i, i+2);
	}

	for(auto i = 0; i < parVecMap.GetLocalSize(); i++){
		vec.AddValueLocal(i, 2);
	}

	vec.VecAdd(vec2);
	vec.VecScale(2);
	int dot = vec.VecDot(vec2);
 	//vec.VecView();

 	std::cout << "dot = " << dot << std::endl;

 	//auto spec1 = specNonHerm<std::complex<double>, int>(parVecMap, "v2.txt");
 	auto spec1 = specNonHerm<std::complex<double>, int>(parVecMap, "v2.txt");
	//spec1.VecView();
    
    auto spec2 = specNonSymm<double, int>(parVecMap, "v1.txt");
    //spec2.VecView();

    auto spec3 = specNonSymmCplex<double, int>(parVecMap, "v3.txt");
    spec3.VecView();

    //for(auto i = 0 ; i < spec3.size(); i++){
    //	std::cout << spec3[i] << std::endl;
    //}

	MPI_Finalize();

	return 0;
}