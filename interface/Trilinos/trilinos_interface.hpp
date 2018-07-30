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

#if defined(__USE_COMPLEX__)
typedef std::complex<double>               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;

#else

typedef double                			   ST;
typedef Teuchos::ScalarTraits<ST>          SCT;

#endif

RCP<CrsMatrix<ST> > ConvertToTrilinosMat(parMatrixSparse<ST, int > *M){

	RCP<const Teuchos::Comm<int> > comm =	
    	Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

	parVectorMap<int> *pmap = M->GetYMap();

	int gsize = pmap->GetGlobalSize();

	RCP<const Map<> > map = rcp (new Map<> (gsize, 0, comm));

    RCP<CrsMatrix<ST> > K = rcp (new CrsMatrix<ST> (map, 0,  Tpetra::DynamicProfile));

    const int numMyElements = map->getNodeNumElements();
    
    ArrayView<const int> myGlobalElements = map->getNodeElementList ();

    int i, j;

    int col;
 	ST val;

 	typename std::map<int ,ST>::iterator it;

 	std::map<int,ST> *dynloc;
 	dynloc = M->GetDynMatLoc();
 	int container_size;
 	for(i = 0; i < numMyElements; i++){
 		for(it = dynloc[i].begin(); it != dynloc[i].end(); ++it){
 			col = it->first;
 			val = it->second;

	 		K->insertGlobalValues(myGlobalElements[i],tuple<int>(col),tuple<ST>(val));
 		}
 	}

 	K->fillComplete ();

 	return K;


}