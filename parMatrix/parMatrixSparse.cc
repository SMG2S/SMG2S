#include "parMatrixSparse.h"

template<typename T, typename S>
parMatrixSparse<T,S>::parMatrixSparse()
{
	dynmat_lloc = NULL;
	dynmat_gloc = NULL;

	CSR_lloc = NULL; 
	CSR_gloc = NULL;
	
	nnz_lloc = 0;
	nnz_gloc = 0;

	ncols = 0;
	nrows = 0;
	njloc = 0;

	lower_x = 0;
	lower_y = 0;
	upper_x = 0;
	upper_y = 0;

	x_index_map = NULL;
	y_index_map = NULL;

	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	VNumRecv = NULL;
	VNumSend = NULL;
	Rbuffer = NULL;
	Sbuffer = NULL;
	DTypeRecv = NULL;
	DTypeSend = NULL;
}

