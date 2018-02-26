#include <mpi.h> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <string> 
#include <vector>

#include "../parVector/parVector.h"
#include "MatrixCSR.h"

template<typename T, typename S>
class parMatrixSparse
{
	private:
		std::map<S, T> *dynamat_lloc, *dynmat_gloc;
		//size of local matrix
		S	ncols, nrows;
		
		MatrixCSR *CSR_lloc, *CSR_gloc;
		
		S	nnz_lloc, nnz_gloc;
		
		parVectorMap	*x_index_map;
		parVectorMap	*y_index_map;

		S	njloc;
		S	lower_x, lower_y, upper_x, upper_y;

		// mpi size and rank 
		int ProcID, nProcs;
		
		// pointers and buffers for setting up communications
		S	*VNumRecv, *VNumSend;
		S	maxRecv, maxSend;
		S	**Rbuffer, **Sbuffer;

		MPI_Datatype DTypeRecv , DTypeSend ;
			
	public:
		//constructor
		parMatrixSparse();
		parMatrixSparse(parVector *XVec, parVector *YVec);
		//deconstructor
		~parMatrixSparse();

		//get
		parVectorMap *GetXMap(){return x_index_map;};
		ParVectorMap *GetYMap(){return y_index_map;};

		S	GetXLowerBound();
		S	GetYlowerBound();

	        S       GetXUpperBound();
                S       GetYUpperBound();
		
		void	GetTrueLocalSize(S &rs, S &cs){
			rs = nrows;
			cs = njloc;
		};

		void	GetLocalSize(S &rs, S &cs){
			rs = nrows;
			cs = ncols;
		};

		std::map<S,T>	*GetDynMatGLobLoc(){return dynmat_lloc;};
		std::map<S,T>	*GetDynMatGlobLoc(){return dynmat_gloc;};

		MatrixCSR	*GetCSRLocLoc(){return CSR_lloc;};
		MatrixCSR	*GetCSRGlobLoc(){return CSR_gloc};

		//add
		void	AddValueLocal( S row, S col, T value);
  		void 	AddValuesLocal( S nindex, S *rows, S *cols, T *values);

		//Reader
		void	ReadExtMat();
		// convert from dyn to csr
		void	ConvertToCSR();

		void	FindColsToRecv();
		void	SetupDataTypes();

		void	TestCommunication(parVector	*XVec, parVec *YVec);

		//spmv
		void	MatVecProd(parVector *XVec, parVector *YVec);

};
