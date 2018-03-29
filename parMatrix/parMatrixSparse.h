#include <mpi.h> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <string> 
#include <vector>

#include "../parVector/parVector.cc"
//#include "MatrixCSR.h"
#include "MatrixCSR.h"


template<typename T, typename S>
class parMatrixSparse
{
	private:
		std::map<S, T> *dynmat_lloc, *dynmat_gloc;
		//size of local matrix
		S	ncols, nrows;
		
		MatrixCSR<T,S> *CSR_lloc, *CSR_gloc;
		
		S	nnz_lloc, nnz_gloc;
		
		parVectorMap<S>	*x_index_map;
		parVectorMap<S>	*y_index_map;

		S	njloc;
		S	lower_x, lower_y, upper_x, upper_y;

		// mpi size and rank 
		int ProcID, nProcs;
		
		// pointers and buffers for setting up communications
		S	*VNumRecv, *VNumSend;
		S	maxRecv, maxSend;
		S	**Rbuffer, **Sbuffer;

		MPI_Datatype *DTypeRecv , *DTypeSend ;
			
	public:
		//constructor
		parMatrixSparse();
		parMatrixSparse(parVector<T,S> *XVec, parVector<T,S> *YVec);
		//deconstructor
		~parMatrixSparse();

		//get
		parVectorMap<S> *GetXMap(){return x_index_map;};
		parVectorMap<S> *GetYMap(){return y_index_map;};

		S	GetXLowerBound();
		S	GetYLowerBound();

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

		MatrixCSR<T,S>	*GetCSRLocLoc(){return CSR_lloc;};
		MatrixCSR<T,S>	*GetCSRGlobLoc(){return CSR_gloc;};

		//add
		void	AddValueLocal( S row, S col, T value);
  		void 	AddValuesLocal( S nindex, S *rows, S *cols, T *values);

		//global add
		void	AddValue(S row, S col, T value);

  		//set
		void	SetValueLocal( S row, S col, T value);
  		void 	SetValuesLocal( S nindex, S *rows, S *cols, T *values);

		//global set
		void	SetValue(S row, S col, T value);

		//get
		T	GetLocalValue(S row, S col);
		T	GetValue(S row, S col);

		//show

		void	MatView();

		//set mat diagonal by vector given 
		void	SetDiagonal(parVector<T,S> *diag);

		//Mat Scala
		void	Scale(T scale);
		
		//AXPY
		void	AXPY(parMatrixSparse<T,S> *X, T scale);

		//AYPX
		void    AYPX(parMatrixSparse<T,S> *X, T scale);

		//Reader
		void	ReadExtMat();
		// convert from dyn to csr
		void	ConvertToCSR();

		void	FindColsToRecv();
		void	SetupDataTypes();

		void	TestCommunication(parVector<T,S> *XVec, parVector<T,S> *YVec);

		//spmv
		void	MatVecProd(parVector<T,S> *XVec, parVector<T,S> *YVec);

		//spgmm
		void	MatMatProd(parVector<T,S> *XVec, parVector<T,S> *YVec, parVector<T,S> *ZVec);


};
