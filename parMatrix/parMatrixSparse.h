#ifndef __PAR_MATRIX_SPARSE_H__
#define __PAR_MATRIX_SPARSE_H__

#include <mpi.h> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <string> 
#include <vector>
#include "../utils/MPI_DataType.h"
#include "../parVector/parVector.cc"
//#include "../utils/utils.h"
#include "MatrixCSR.h"

#include <complex>

template<typename T, typename S>
class parMatrixSparse
{
	private:
		std::map<S, T> *dynmat_lloc, *dynmat_gloc;

		std::map<S, T> *dynmat_loc;

		//size of local matrix
		S	ncols, nrows;
		
		MatrixCSR<T,S> *CSR_lloc, *CSR_gloc;
		
		S	nnz_lloc, nnz_gloc, nnz_loc;
		
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

	    S   GetXUpperBound();
        S   GetYUpperBound();
		
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
		T		GetLocalValue(S row, S col);
		T		GetValue(S row, S col);


		//combine gloc + lloc -> loc together
		void	glocPlusLloc();

		void	llocToGlocLoc();

		//show
		void	MatView();

		void	LOC_MatView();


  		//LOC set
		void	Loc_SetValueLocal( S row, S col, T value);
  		void 	Loc_SetValuesLocal( S nindex, S *rows, S *cols, T *values);

		//LOC global set
		void	Loc_SetValue(S row, S col, T value);

		//LOC get
		T		Loc_GetLocalValue(S row, S col);
		T		Loc_GetValue(S row, S col);

		//set mat diagonal by vector given 
		void	SetDiagonal(parVector<T,S> *diag);

		//Loc set mat diagonal by vector given 
		void	Loc_SetDiagonal(parVector<T,S> *diag);


		//Mat Scala
		void	MatScale(T scale);
		
		//AXPY
		void	MatAXPY(parMatrixSparse<T,S> *X, T scale);

		//AYPX
		void    MatAYPX(parMatrixSparse<T,S> *X, T scale);

		//Loc Mat Scala
		void	Loc_MatScale(T scale);
		
		//Loc AXPY
		void	Loc_MatAXPY(parMatrixSparse<T,S> *X, T scale);

		//Loc AYPX
		void    Loc_MatAYPX(parMatrixSparse<T,S> *X, T scale);


		void 	CSRMatView();

		//Reader
		void	ReadExtMat();

		//Writer
		void	WriteExtMat();

		// convert from dyn to csr
		void	ConvertToCSR();

		// Zeros all entries with keeping the previous matrix pattern
		void	ZeroEntries();

		// Loc: Zeros all entries with keeping the previous matrix pattern
		void	Loc_ZeroEntries();

		void	FindColsToRecv();

		void	SetupDataTypes();

		void	TestCommunication(parVector<T,S> *XVec, parVector<T,S> *YVec);

		//spmv
		void	CSR_MatVecProd(parVector<T,S> *XVec, parVector<T,S> *YVec);

		void	ELL_MatVecProd(parVector<T,S> *XVec, parVector<T,S> *YVec);
		//spgmm
		void	MatMatProd(parVector<T,S> *XVec, parVector<T,S> *YVec, parVector<T,S> *ZVec);
		
		//matrix multiple a special nilpotent matrix
		void	MA(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod);

		//special nilpotent matrix multiple another matrix
		void	AM(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod);

};


#endif