//#include "parVectorMap.cc"
#include "../parMatrix/parMatrixSparse.cc"
#include <math.h>
#include <complex.h>
#include <sys/malloc.h>

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

    //MPI_Barrier(MPI_COMM_WORLD);
    
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);


    int probSize = 11;
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
    parVector<double,int> *prod2 = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);

    vec->SetTovalue(a); //1.0,1.0...1.0

    prod->SetTovalue(0.0);

    prod2->SetTovalue(0.0);

    MPI_Barrier(MPI_COMM_WORLD);

    vec->VecAdd(prod); //2.0,2.0...
   
    MPI_Barrier(MPI_COMM_WORLD);

    vec->VecScale(c); //4.0,4.0...

    double dot_local, dot;

    dot_local = vec->VecDot(prod);

    MPI_Allreduce(&dot_local, &dot, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    if(world_rank == 0){printf("vecdot = %f\n", dot);}

    MPI_Barrier(MPI_COMM_WORLD);

//    int as = vec->GetArraySize();

//    if(world_rank == 0){printf("vector array Size  = %d\n", as);}

    MPI_Barrier(MPI_COMM_WORLD);

    double *va = vec->GetArray();

    MPI_Barrier(MPI_COMM_WORLD);

    vec->VecScale(c); //8.0,8.0...

    vec->VecView();

    //Matrix Initialization

    parMatrixSparse<double,int> *Am = new parMatrixSparse<double,int>(vec,prod);
    parMatrixSparse<double,int> *Bm = new parMatrixSparse<double,int>(vec,prod);

    if(world_rank == 0){printf("Matrix Initialized\n");}

    MPI_Barrier(MPI_COMM_WORLD);

//    printf("\n\nprc %d: LOWER_X = %d, LOWER_Y = %d, UPPER_X = %d, UPPER_Y = %d \n", world_rank, Am->GetXLowerBound(), Am->GetYLowerBound(), Am->GetXUpperBound(), Am->GetYUpperBound());


    for(int j=6; j < probSize; j++){
        Am->SetValue(j-6,j,1);
        Am->SetValue(j, j-6,2);
    }

    for(int j=0; j < probSize; j++){
        Am->SetValue(j,j,1);
    }

    for(int j=0; j < probSize; j++){
        Bm->SetValue(j,j,0);
    }
/*
    for(int j=0; j < probSize - 5; j++){
        Am->SetValue(j,j+5,1.0);
    }
*/
    double x, y;

 //   x = Am->GetValue(10,11);

   // y = Am->GetValue(10,11);

    //printf("Prc %d: x = %f, y = %f \n", world_rank, x, y);

    MPI_Barrier(MPI_COMM_WORLD);

    Am->MatView();

    MPI_Barrier(MPI_COMM_WORLD);

    Am->ConvertToCSR();
	
    //Am->AXPY(Bm, 2.0);

    MPI_Barrier(MPI_COMM_WORLD);

    Am->MatView();

    MPI_Barrier(MPI_COMM_WORLD);

    Am->FindColsToRecv();
    
    MPI_Barrier(MPI_COMM_WORLD);

    Am->SetupDataTypes();

    MPI_Barrier(MPI_COMM_WORLD);

 //   Am->TestCommunication(vec,prod);

    Am->ELL_MatVecProd(vec,prod);

    MPI_Barrier(MPI_COMM_WORLD);

//    if(world_rank == 0){printf("print SPMV results\n");}
 
    prod->VecView();

    MPI_Barrier(MPI_COMM_WORLD);

    Am->CSR_MatVecProd(vec,prod2);

    MPI_Barrier(MPI_COMM_WORLD);

    prod2->VecView();

    MPI_Finalize();

    return 0;
}

