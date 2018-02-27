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

template<typename T, typename S>
parMatrixSparse<T,S>::parMatrixSparse(parVector *XVec, parVector *YVec)
{
	dynmat_lloc = NULL;
	dynmat_gloc = NULL;

	CSR_lloc = NULL;
	CSR_gloc = NULL;

	nnz_lloc = 0;
	nnz_gloc = 0;

	x_index_map = NULL;
	y_index_map = NULL;

	ncols = 0;
	nrows = 0;
	njloc = 0;

	lower_x = 0;
	lower_y = 0;

	upper_x = 0;
	upper_y = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	VNumRecv = NULL;
	VNumSend = NULL;
	Rbuffer = NULL;
	Sbuffer = NULL;

	DTypeRecv = NULL;
	DTypeSend = NULL;

	//get vector map for x and y direction
	x_index_map = XVec->GetVectorMap();
	x_index_map->AddUser();
	y_index_Map = YVec->GetVectorMap();
	y_index_map->AddUser();
	
	if(x_index_map != NULL && y_index_map != NULL){
		//get num of rows and cols in this mpi procs
		ncols = x_index_map->GetGlobalSize();
		nrows = y_index_map->GetLocalSize();
		njloc = x_index_map->GetLocalSize();
		//get upper and lower bounds
		lower_x = x_index_map->GetLowerBound();
		lower_y = y_index_map->GetLowerBound();
		upper_x = x_index_map->GetUpperBound();
		upper_y = y_index_map->GetUpperBound();
	}
}

template<typename T, typename S>
parMatrixSparse<T,S>::~parMatrixSparse()
{
	//if index map is defined
	if(x_index_map != NULL){
		x_index_map->DeleteUser();
		
//		if(x_index_map->GetUser() == 0){
			delete x_index_map;
//		}
	}


	if(y_index_map != NULL){
		y_index_map->DeleteUser();

//		if(x_index_map->GetUser() == 0){
			delete y_index_map;
//		}
	}

	//if dynmat has been defined
	if(dynmat_lloc != NULL){
		delete [] dynmat_lloc;
	}
	if(dynmat_gloc != NULL){
		delete [] dynmat_gloc;
	}
	if(CSR_lloc != NULL){
		delete [] CSR_lloc;
	}
	if(CSR_gloc != NULL){
		delete [] CSR_gloc;
	}
	if(VNumRecv != NULL){
		delete [] VNumRecv;
	}
	if(VNumSend != NULL){
		delete [] VNumSend;
	}
	if(Rbuffer != NULL){
		int i;
		for(i = 0; i < nProcs; i++){
			if (Rbuffer[i] != NULL){
				delete [] Rbuffer[i];
			}
		}
		delete [] Rbuffer;
	}

	if(Sbuffer != NULL){
		int i;
	        for(i = 0; i < nProcs; i++){
                        if (Sbuffer[i] != NULL){
                                delete [] Sbuffer[i];
                        }
                }
                delete [] Sbuffer;
	}


	if(DTypeRecv != NULL){
		int i;
		for(i = 0; i < nProcs; i++){
			if(DTypeSend[i] != MPI_DATATYPE_NULL){
				MPI_Type_free(&DTypeSend[i]);
			}
		}
		delete [] DTypeSend;
	}
}

template<typename T, typename S>
S parMatrixSparse<T,S>::GetXLowerBound(){
	if(x_index_map != NULL){
		return x_index_map->GetLowerBound();
	}
	else{
		return 0;
	}
}

template<typename T, typename S>
S parMatrixSparse<T,S>::GetXUpperBound(){
        if(x_index_map != NULL){
                return x_index_map->GetUpperBound();
        }
        else{
                return 0;
        }
}

template<typename T,typename S>
S parMatrixSparse<T,S>::GetYLowerBound(){
        if(y_index_map != NULL){
                return y_index_map->GetLowerBound();
        }
        else{
                return 0;
        }
}

template<typename T, typename S>
S parMatrixSparse<T,S>::GetYUpperBound(){
        if(y_index_map != NULL){
                return y_index_map->GetUpperBound();
        }
        else{
                return 0;
        }
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValueLocal(S row, S col, T value)
{
	std::map<S,T>::iterator it;
	//if location is inside of local area then add to local dynamic map
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col > 0)){
		if(dynmat_lloc == NULL){
			dynmat_lloc = new std::map<S,T> [nrows];
		}
		it = dynmat_lloc[row].find(col);
		if(it == dynmat_lloc[row].end()){
			dynmat_lloc[row][col] = value;
			nnz_lloc++;
		}
		else{
			it->second = it->second + value;
		}
	//if location is inside of local-global area
	}
	else if ((row < nrows && row >= 0) && (col >= upper_x || col < lower_x) && (col >= 0)){
		if(dynmat_gloc == NULL){
			dynmat_gloc = new std::map<S,T> [nrows];
		}
		it = dynmat_gloc[row].find(col);
		if(it == dynmat_gloc[row].end()){
			dynmat_gloc[row][col] = value;
			nnz_gloc++;
		}
		else{
			it->second = it->second + value;
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValuesLocal(S nindex, S *rows, S *cols, T *values)
{
	std::map<S,T>::iterator it;
	
	for( S i = 0; i < nindex; i++){
		AddValueLocal(rows[i],cols[i],values[i]);
	}
}


