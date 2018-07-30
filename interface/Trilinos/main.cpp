//#include "parVectorMap.cc"
#include "../../parMatrix/parMatrixSparse.h"
#include <math.h>
#include <complex.h>
#include "../../utils/MPI_DataType.h"
#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

using Tpetra::CrsMatrix;
using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::arcp;
using Teuchos::tuple;
using Tpetra::Map;
using std::vector;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::tuple;
using std::cout;
using std::endl;

typedef double                			   ST;
typedef Teuchos::ScalarTraits<ST>          SCT;

int main(int argc, char** argv) {

	Teuchos::oblackholestream blackhole;

	Teuchos::GlobalMPISession mpisess (&argc, &argv, &std::cout);

    // Initialize the MPI environment
    //MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;

    int probSize = 11;

    double start, end;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

	RCP<const Teuchos::Comm<int> > comm =	
    	Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

	const int MyPID = comm->getRank ();

	Teuchos::oblackholestream blackHole;
  		std::ostream& out = std::cout;

  	RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

	RCP<const Map<> > map = rcp (new Map<> (probSize, 0, comm));


	int min = map->getMinGlobalIndex();
	int max = map->getMaxGlobalIndex();

    //MPI_Barrier(MPI_COMM_WORLD);
    
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);


    int span, lower_b, upper_b;

    //span = int(floor(double(probSize)/double(world_size)));
	span = int(ceil(double(probSize)/double(world_size)));

    printf("span = %d\n", span);
    if(world_rank == world_size - 1){
	lower_b = world_rank * span;
	upper_b = probSize - 1 + 1;
    }
    else{
	lower_b = world_rank * span;
	upper_b = (world_rank + 1) * span - 1 + 1;
    }

    printf("Proc. %d   Lower bound = %d   Upper bound = %d, map of Trilinos is min global index = %d, max global index = %d \n",world_rank, lower_b , upper_b, min, max) ;

    double a = 1.0, c = 2.0, b =2.0;

    parVector<double,int> *vec = new parVector<double,int>(MPI_COMM_WORLD, lower_b, upper_b);

    vec->SetTovalue(a);
    vec->VecView();

    parMatrixSparse<double,int> *Am = new parMatrixSparse<double,int>(vec,vec);
    Am->Loc_SetDiagonal(vec);

    double rnd;

    for(int i = 0; i < probSize; i++){
        for(int j = i - 3; j < i; j++){
            if(j >= 0){
            	rnd = 0.05;
            	Am->Loc_SetValue(i,j,rnd);
            }
        }
    }

    Am->LOC_MatView();

    RCP<CrsMatrix<ST> > K = rcp (new CrsMatrix<ST> (map, 0,  Tpetra::DynamicProfile));

    const int numMyElements = map->getNodeNumElements();
    
    ArrayView<const int> myGlobalElements = map->getNodeElementList ();

/*
    int dd = 3, m;
    int *ind = new int[dd];
    double *vv = new double[dd];
    for(m = 0; m < dd; m++){
    	ind[m] = m;
    	vv[m] = 1.0;
    }

    ArrayView<const int> cols(ind, dd);
    ArrayView<const double> vals(vv, dd);
*/

    printf("myProcId = %d: numMyElements = %d\n", MyPID, numMyElements);

    int i, j;
/*
    for(i = 0; i < numMyElements; i++){
    	//K->insertGlobalValues(myGlobalElements[i],tuple<int>(myGlobalElements[i]),tuple<double>(1.0));
    	K->insertGlobalValues(myGlobalElements[i],cols,vals);
    }
 */	
 
 	int col;
 	double val;

 	typename std::map<int ,double>::iterator it;

 	std::map<int,double> *dynloc;
 	dynloc = Am->GetDynMatLoc();
 	int container_size;
 	for(i = 0; i < numMyElements; i++){
 		for(it = dynloc[i].begin(); it != dynloc[i].end(); ++it){
 			col = it->first;
 			val = it->second;

	 		K->insertGlobalValues(myGlobalElements[i],tuple<int>(col),tuple<double>(val));
 		}
 	}

 	K->fillComplete ();

 	K->describe(*fos, Teuchos::VERB_EXTREME);

    //MPI_Finalize();

    return 0;
 }