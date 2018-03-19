//#include "parVectorMap.cc"
#include "../parMatrix/parMatrixSparse.cc"
#include <math.h>
#include <complex.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

 //   parVectorMap<int> pm(MPI_COMM_WORLD, 1,2);

    //MPI_Barrier(MPI_COMM_WORLD);
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);


    int probSize = 10;
    int span, lower_b, upper_b;

    span = int(floor(double(probSize)/double(world_size)));

    printf("span = %d\n", span);
    if(world_rank == world_size - 1){
	lower_b = world_rank * span;
	upper_b = probSize - 1 + 1;
    }
    else{
	lower_b = world_rank * span;
	upper_b = (world_rank + 1) * span - 1 + 1;
    }

    printf("Proc. %d   Lower bound = %d   Upper bound = %d \n",world_rank, lower_b , upper_b ) ;

    double a = 1.0, c = 2.0; float b =2.0;

    parVector<double,int> *vec = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);
    parVector<double,int> *prod = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);

    vec->SetTovalue(a);

    prod->SetTovalue(c);

    MPI_Barrier(MPI_COMM_WORLD);

    vec->VecView();

    MPI_Barrier(MPI_COMM_WORLD);

    vec->VecAdd(prod);
   
    vec->VecView();

    MPI_Barrier(MPI_COMM_WORLD);

    vec->VecScale(c);

    double dot;

    dot = vec->VecDot(prod);

    if(world_rank == 0){printf("vecdot = %f\n", dot);}

    //Matrix Initialization

    parMatrixSparse<double,int> *Am = new parMatrixSparse<double,int>(vec,prod);

    if(world_rank == 0){printf("Matrix Initialized\n");}

    Am->SetDiagonal(vec);

    //Am->AddValueLocal(3,5,10.0);

    double x;

    x = Am->GetValue(0,0);

//    Am->MatView();

    if(world_rank == 0){printf("x = %f \n", x);}

    //Am->ConvertToCSR();
	
    //Am->FindColsToRecv();
    
   // Am->SetupDataTypes();

   // MPI_Barrier(MPI_COMM_WORLD);

   // Am->MatVecProd(vec,prod);

  //  if(world_rank == 0){printf("print SPMV results\n");}
 
    //prod->VecView();

//    parVector<float,int> *vec2 = new parVector<float,int>(MPI_COMM_WORLD, lower_b, upper_b);
//    vec2->SetTovalue(b);

    MPI_Finalize();

    return 0;
}

