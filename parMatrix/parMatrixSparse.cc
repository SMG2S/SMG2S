#include "parMatrixSparse.h"

template<typename T, typename S>
parMatrixSparse<T,S>::parMatrixSparse()
{
	dynmat_lloc = NULL;
	dynmat_gloc = NULL;
	dynmat_loc  = NULL;

	CSR_lloc = NULL; 
	CSR_gloc = NULL;
	
	nnz_lloc = 0;
	nnz_gloc = 0;
	nnz_loc = 0;

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
parMatrixSparse<T,S>::parMatrixSparse(parVector<T,S> *XVec, parVector<T,S> *YVec)
{
	dynmat_lloc = NULL;
	dynmat_gloc = NULL;
	dynmat_loc  = NULL;

	CSR_lloc = NULL;
	CSR_gloc = NULL;

	nnz_lloc = 0;
	nnz_gloc = 0;
	nnz_loc  = 0;

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
	x_index_map = XVec->GetVecMap();
	x_index_map->AddUser();
	y_index_map = YVec->GetVecMap();
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
	typename std::map<S,T>::iterator it;
	//if location is inside of local area then add to local dynamic map
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col >= 0)){
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
T parMatrixSparse<T,S>::GetLocalValue(S row, S col)
{
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col >= 0)){
		return dynmat_lloc[row][col]; 
	}
	else if ((row < nrows && row >= 0) && (col >= upper_x || col < lower_x) && (col >= 0)){
		return dynmat_gloc[row][col];
	}
	else return 0.0;
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValuesLocal(S nindex, S *rows, S *cols, T *values)
{
	typename std::map<S,T>::iterator it;
	
	for( S i = 0; i < nindex; i++){
		AddValueLocal(rows[i],cols[i],values[i]);
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AddValue(S row, S col, T value)
{

	if((row >= lower_y) && (row < upper_y) && (col < ncols)){
		AddValueLocal(y_index_map->Glob2Loc(row),col,value);
	}
}

//set
template<typename T,typename S>
void parMatrixSparse<T,S>::SetValueLocal( S row, S col, T value)
{
	typename std::map<S,T>::iterator it;
	//if location is inside of local area then add to local dynamic map
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col >= 0)){
		if(dynmat_lloc == NULL){
			dynmat_lloc = new std::map<S,T> [nrows];
		}
		it = dynmat_lloc[row].find(col);
		if(it == dynmat_lloc[row].end()){
			dynmat_lloc[row][col] = value;
			nnz_lloc++;
		}
		else{
	//		it->second = it->second + value;
			it->second = value;
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
//			it->second = it->second + value;
			it->second = value;
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::SetValuesLocal( S nindex, S *rows, S *cols, T *values)
{
	typename std::map<S,T>::iterator it;
	
	for( S i = 0; i < nindex; i++){
		SetValueLocal(rows[i],cols[i],values[i]);
	}
}

//global set
template<typename T,typename S>
void parMatrixSparse<T,S>::SetValue(S row, S col, T value)
{
	S local_row = y_index_map->Glob2Loc(row);

	if((row >= lower_y) && (row < upper_y) && (col < ncols)){
		SetValueLocal(local_row,col,value);
	}
}


template<typename T,typename S>
T parMatrixSparse<T,S>::GetValue(S row, S col)
{
	S local_row = y_index_map->Glob2Loc(row);
	if(local_row >= 0 && local_row < nrows){
		return GetLocalValue(local_row, col);
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::glocPlusLloc(){
	
	S i,j;
	T v;
	typename std::map<S,T>::iterator it;

	if(ProcID == 0) {std::cout << "Combine the block-diagonal part and non block-diagonal part of parallel matrix together " << std::endl;}

	if(dynmat_loc == NULL){
		dynmat_loc = new std::map<S, T> [nrows];
	}
	for(i = 0; i < nrows; i++){
		if((dynmat_gloc != NULL) && (dynmat_lloc != NULL)){
			dynmat_loc[i].insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
			dynmat_loc[i].insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());
		}
		else if ((dynmat_lloc != NULL)){
			dynmat_loc[i].insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
		}
		else if (dynmat_gloc != NULL){
			dynmat_loc[i].insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());
		}
		else{
			if(ProcID == 0) {std::cout << "Cannot execute the function glocPlusLloc since the given matrix is NULL " << std::endl;}
		}
	}
	nnz_loc = nnz_gloc + nnz_lloc;

}

template<typename T,typename S>
void parMatrixSparse<T,S>::llocToGlocLoc()
{
	typename std::map<S,T>::iterator it;

	S col;

	//if location is inside of local area then add to local dynamic map
	if(dynmat_loc != NULL){
		if(dynmat_lloc == NULL){
			dynmat_lloc = new std::map<S,T> [nrows];
		}
		if(dynmat_gloc == NULL){
			dynmat_gloc = new std::map<S,T> [nrows];
		}

		for(S i = 0; i < nrows; i++){		
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
				col = it->first;
				if(col < upper_x && col >= lower_x && col >= 0){
					dynmat_lloc[i][col] = it->second;
				}
				else if((col >= upper_x || col < lower_x) && (col >= 0)){
					dynmat_gloc[i][col] = it->second;
				}
			}
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::MatView(){
	
	S i, j;
	T v;
	typename std::map<S,T>::iterator it;

	if(ProcID == 0) {std::cout << "Parallel MatView: " << std::endl;}

	for (i = 0; i < nrows; i++){
		std::map<S,T> merge;
		std::cout << "row " << y_index_map->Loc2Glob(i) << ": ";

		if((dynmat_gloc != NULL) && (dynmat_lloc != NULL)){
			merge.insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
			merge.insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());
	
			for(it = merge.begin(); it != merge.end(); ++it){
				std::cout <<"("<<it->first << "," << it->second << "); ";
			}
			merge.clear();
		}
		else if ((dynmat_lloc != NULL)){
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); ++it){
				std::cout <<"("<<it->first << "," << it->second << "); ";
			}
		}
		else if (dynmat_gloc != NULL){
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); ++it){
				std::cout <<"("<<it->first << "," << it->second << "); ";
			}
		}

		std::cout << std::endl;
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::LOC_MatView(){
	
	S i, j;
	T v;
	typename std::map<S,T>::iterator it;

	if(ProcID == 0) {std::cout << "LOC MODE Parallel MatView: " << std::endl;}

	for (i = 0; i < nrows; i++){
		if(dynmat_loc != NULL){
			std::cout << "row " << y_index_map->Loc2Glob(i) << ": ";
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
				std::cout <<"("<<it->first << "," << it->second << "); ";
			}
		}
	std::cout << std::endl;
	}
}


//Loc set
template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetValueLocal( S row, S col, T value)
{
	typename std::map<S,T>::iterator it;

	if(dynmat_loc == NULL){
		dynmat_loc = new std::map<S,T> [nrows];
	}
	it = dynmat_loc[row].find(col);

	if(it == dynmat_loc[row].end()){
		dynmat_loc[row][col] = value;
		nnz_loc++;
	}
	else{
		it->second = value;
	}	
}

template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetValuesLocal( S nindex, S *rows, S *cols, T *values)
{
	typename std::map<S,T>::iterator it;
	
	for( S i = 0; i < nindex; i++){
		Loc_SetValueLocal(rows[i],cols[i],values[i]);
	}
}

//Loc global set
template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetValue(S row, S col, T value)
{
	S local_row = y_index_map->Glob2Loc(row);

	if (local_row >= 0 && local_row < nrows){
		Loc_SetValueLocal(local_row,col,value);
	}

}

template<typename T,typename S>
T parMatrixSparse<T,S>::Loc_GetLocalValue(S row, S col)
{
	if(dynmat_loc != NULL){
		return dynmat_loc[row][col]; 
	}
	else return 0.0;
}


template<typename T,typename S>
T parMatrixSparse<T,S>::Loc_GetValue(S row, S col)
{
	S local_row = y_index_map->Glob2Loc(row);
	if(local_row >= 0 && local_row < nrows){
		return Loc_GetLocalValue(local_row, col);
	}
}



template<typename T,typename S>
void parMatrixSparse<T,S>::SetDiagonal(parVector<T,S> *diag)
{
	S    lbound;
	S 	 ubound;

	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set diagonal for non-square matrix.");
		}
	}
	else{
		T *a = diag->GetArray();
		S local_size = diag->GetArraySize();

		lbound = diag->GetLowerBound();
		ubound = diag->GetUpperBound();

		std::cout << "local size  = " << local_size << "  lower bound = " << lbound << "  upper bound = " << ubound << std::endl;
		for(S i = 0; i < local_size; i++){
			SetValueLocal(i,diag->Loc2Glob(i),a[i]);
		}

	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_SetDiagonal(parVector<T,S> *diag)
{
	S    lbound;
	S 	 ubound;

	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set diagonal for non-square matrix.");
		}
	}
	else{
		T *a = diag->GetArray();
		S local_size = diag->GetArraySize();

		lbound = diag->GetLowerBound();
		ubound = diag->GetUpperBound();

		std::cout << "local size  = " << local_size << "  lower bound = " << lbound << "  upper bound = " << ubound << std::endl;
		for(S i = 0; i < local_size; i++){
			Loc_SetValueLocal(i,diag->Loc2Glob(i),a[i]);
		}

	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::ConvertToCSR()
{
	S 	count, i, j, k;
	T	v;
	typename std::map<S,T>::iterator it;
	
	if(dynmat_lloc != NULL){
		//allocate csr matrix


		CSR_lloc = new MatrixCSR<T,S>(nnz_lloc, nrows);

//		CSR_lloc = new MatrixCSR<T,S>(nrows);

		//convert local local to csr
		count = 0;
		CSR_lloc->rows[0] = 0;

		for(i = 0; i < nrows; i++){
			CSR_lloc->rows[i] = count;
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); it++){
				j = it->first;
				v = it->second;
				CSR_lloc->vals[count] = v;
				CSR_lloc->cols[count] = j;
				count++;
			}
		}
		CSR_lloc->rows[nrows] = nnz_lloc;
	}

	if(dynmat_gloc != NULL){


		CSR_gloc = new MatrixCSR<T,S>(nnz_gloc, nrows);
//		CSR_gloc = new MatrixCSR<T,S>(nrows);

		//convert global-local to CSR
		count = 0;
		CSR_gloc->rows[0] = 0;
		for(i = 0; i < nrows; i++){
			CSR_gloc->rows[i] = count;
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); it++){
				j = it->first;
				v = it->second;
				CSR_gloc->vals[count] = v;
				CSR_gloc->cols[count] = j;
				count++;
			}
		}
		CSR_gloc->rows[nrows] = nnz_gloc;
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::ReadExtMat()
{
	//Reader
	std::ifstream file("matrix.mat");
	std::string line;
	S	row, col;
	T	value;
	row = 0;
	col = 0;
	value = 0.0;
	S quit = 0;

	if(ProcID == 0) printf("start reading matrix !!!!\n");
	//Past no number content in the file
	while(std::getline(file,line)){
		row = 0;
		col = 0;
		value = 0.0;
		std::stringstream linestream(line);
		linestream >> row >> col >> value;
		if(row != 0 && col != 0 && value != 0.0){
			break;
		}
	}
		
	//read values and add them to matrix
	while(std::getline(file,line)){
		std::stringstream linestream(line);
		linestream >> row >> col >> value;
		row = row - 1;
		col = col - 1;
		if((row >= lower_y && row < upper_y) && (col < ncols)){
			AddValueLocal(y_index_map->Glob2Loc(row),col,value);
		}
	}
	if(ProcID == 0) printf("The matrix is readed !!!!\n");
}

template<typename T,typename S>
void parMatrixSparse<T,S>::WriteExtMat(){

}


template<typename T,typename S>
void parMatrixSparse<T,S>::FindColsToRecv()
{
	std::map<S,S> Rrows;
	std::map<S,S> Srows;
	typename std::map<S,S>::iterator vit;
	typename std::map<S,T>::iterator mit;

	S	i, j, k;
	S	count, count1;

	count = 0;
	count1 = 0;

	S	nRecv;

	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	//Initialise vector containing numer of entries to send and recv on each procss, number of vectors

	VNumRecv = new S [nProcs];
	for(i = 0; i < nProcs; i++){
		VNumRecv[i] = 0;
	}
	VNumSend = new S [nProcs];
	for(i = 0; i < nProcs; i++){
		VNumSend[i] = 0;
	}

	MPI_Request	*Rreqs, *Sreqs;
	MPI_Status	status, *Rstat, *Sstat;
	int		tag1, tag2;
	tag1 = 0;
	tag2 = 1;
	int		Rtag, Stag;
	Rtag = 0;
	Stag = 1;

	S		maxRecv, maxSend;

	if(dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(mit = dynmat_gloc[i].begin(); mit != dynmat_gloc[i].end(); mit++){
				j = mit->first;
				vit = Rrows.find(j);
				if(vit == Rrows.end()){
					Rrows[j] = x_index_map->GetOwner(j);
					VNumRecv[Rrows[j]] = VNumRecv[Rrows[j]] + 1;
					//printf("Proc %d: Rrows[%d] = %d ---- %d ----------\n", ProcID,j,Rrows[j],VNumRecv[Rrows[j]]);
					count++;
				}
			}
		}

		nRecv = count;
		//printf("Proc %d: nRecv = %d, VNumRecv[%d] = %d\n", ProcID,count, ProcID, VNumRecv[ProcID]);
	}
//	printf("Proc : %d, %d\n\n\n", ProcID, VNumRecv[0]);
/*
	for(i = 0; i < nProcs; i++){
		printf("Proc %d: -------VNumRecv[%d] = %d\n", ProcID, i,VNumRecv[i]);
	}
*/
	/*
	else{
		printf("fucking check 2\n\n\n\n");
		for(i = 0; i < nProcs; i++){
			VNumRecv[i] = 0;
		}
	}
	*/

	//MPI non-blocking requests and status
	Rreqs = new MPI_Request [nProcs - 1];
	Sreqs = new MPI_Request [nProcs - 1];
	Rstat = new MPI_Status [nProcs - 1];
	Sstat = new MPI_Status [nProcs - 1];

	count = 0;
	maxRecv = 0;
	maxSend = 0;
	

	for(i = 0; i < nProcs; i++){
		if(VNumRecv[i] >maxRecv){
			maxRecv = VNumRecv[i];
		}
		if(i != ProcID){
			MPI_Isend(&VNumRecv[i],1,MPI_INT,i,tag1,MPI_COMM_WORLD,&Sreqs[count]);
			MPI_Irecv(&VNumSend[i],1,MPI_INT,i,tag1,MPI_COMM_WORLD,&Rreqs[count]);
			count++;
		}
	}
	
	//wait for receives to finish
	MPI_Waitall(nProcs-1,Rreqs,Rstat);

/*
	for(i = 0; i < nProcs; i++){
		printf("Proc %d: +++++VNumSend[%d] = %d\n", ProcID, i,VNumSend[i]);
	}
*/

	//find max num to send
	for(i = 0; i < nProcs; i++){
		if(VNumSend[i] > maxSend){
			maxSend = VNumSend[i];
		}
	}

	//Initialisation of send and receive buffers
	Rbuffer = new int * [nProcs];
	Sbuffer = new int * [nProcs];
	
	for(i = 0; i < nProcs; i++){
		Rbuffer[i] = NULL;
	}
	for(i = 0; i < nProcs; i++){
		Sbuffer[i] = NULL;
	}

	//MPI NON-BLOCKING 
	count = 0;
	count1 = 0;
	
	for(i = 0; i < nProcs; i++){
		count = 0;
		if(ProcID != i){
			Sbuffer[i] = new S [VNumRecv[i]];
			for(vit = Rrows.begin(); vit != Rrows.end(); vit++){
				if(vit->second == i){
					Sbuffer[i][count] = vit->first;
					count++;
				}
			}
			MPI_Isend(Sbuffer[i],VNumRecv[i], MPI_INT,i,tag1,MPI_COMM_WORLD, &Sreqs[count1]);
			count++;
		}
		else{
			Sbuffer[i] = new S [1];
			Sbuffer[i][0] = 0;
		}
	}

	count1 = 0;
	
	for(i = 0; i < nProcs; i++){
		if(ProcID != i){
			Rbuffer[i] = new S [VNumSend[i]];
			MPI_Irecv(Rbuffer[i],VNumSend[i],MPI_INT,i,tag1,MPI_COMM_WORLD,&Rreqs[count1]);
		}
		else{
			Rbuffer[i] = new S [1];
			Rbuffer[i][0] = 0;	
		}
	}

	//wait for receives to finish
	MPI_Waitall(nProcs-1, Rreqs, Rstat);

	printf ("Found the cols to recv\n");
	delete [] Rreqs;	
	delete [] Sreqs;
	delete [] Rstat;
	delete [] Sstat;

}


template<typename T,typename S>
void parMatrixSparse<T,S>::SetupDataTypes()
{
	S	i,j,k;
	S	count, *blength, *displac;

	DTypeSend = new MPI_Datatype [nProcs];
	DTypeRecv = new MPI_Datatype [nProcs];

	for(i = 0; i < nProcs; i++){
		count = VNumSend[i];
//		printf("!!!!!!>>>Proc %d::: VNumSend[%d] = %d\n", ProcID, i, count);
		blength = new S [count];
		displac = new S [count];

		for(j = 0; j < count; j++){
//			printf("Rbuffer[i][j] = %d\n", Rbuffer[i][j]);
			blength[j] = 1;
			displac[j] = x_index_map->Glob2Loc(Rbuffer[i][j]);
//			printf("!!!!!!>>>Proc %d: DTypeSend %d //// Rbuffer[i][%d] = %d \n", ProcID, displac[j], j, Rbuffer[i][j]);
		}

		MPI_Type_indexed(count, blength, displac, MPI_DOUBLE, &DTypeSend[i]);
		MPI_Type_commit(&DTypeSend[i]);

		count = VNumRecv[i];
//		printf("VNumRecv[%d] = %d\n", i, count);
		blength = new S [count];
		displac = new S [count];

		for(j = 0; j < count; j++){
			blength[j] = 1;
			displac[j] = Sbuffer[i][j];
//			printf("@@@>>>Proc %d: DTypeRecv %d //// Sbuffer[i][%d] = %d\n", ProcID, displac[j], j,Sbuffer[i][j] );
		}

		MPI_Type_indexed(count, blength, displac, MPI_DOUBLE, &DTypeRecv[i]);
        MPI_Type_commit(&DTypeRecv[i]);
	}

	printf("Setup data types finished\n");	

	delete [] blength;
	delete [] displac;
}

template<typename T,typename S>
void parMatrixSparse<T,S>::TestCommunication(parVector<T,S> *XVec, parVector<T,S> *YVec)
{
	S	i,j,k,l,ng;
	S	sender, receiver;
	sender = 0;
	receiver = 1;
	
	MPI_Status Rstat;

	T	*rBuf, *sBuf;

	k = XVec->GetLocalSize();
	l = XVec->GetGlobalSize();

	//printf("##### k = %d, l = %d\n", k,l);

	rBuf = new T [l];
	for(i = 0; i < l; i++){
		rBuf[i] = 0;	
	}
	
	sBuf = XVec->GetArray();


	int number_amount;

	if(ProcID == sender){
		for(i = 0; i < k; i++){printf("=> sBuf[%d] = %f \n",i, sBuf[i] );}
		MPI_Send(sBuf, 1, DTypeSend[receiver], receiver, 1, MPI_COMM_WORLD);
		printf("sending done! \n");
	}
	else{
		for(i = 0; i < l; i++){printf("*> rBuf[%d] = %f \n",i, rBuf[i] );}
		MPI_Recv(&rBuf,1,DTypeRecv[sender],sender,1,MPI_COMM_WORLD,&Rstat);
		for(i = 0; i < l; i++){printf("@> rBuf[%d] = %f \n",i, rBuf[i] );}

//		MPI_Get_count(&Rstat, DTypeRecv[sender], &number_amount);
//		printf("Message source = %d, %d numbers are received\n", Rstat.MPI_SOURCE, number_amount);
	}

		printf("##### k = %d, l = %d\n", k,l);



	if(ProcID == sender){
		for( i = 0; i < k; i++){
			printf("sBuf[%d] = %f \n", i, sBuf[i]);
		}
	} else{
		for( i = 0; i < l; i++){
            printf("rBuf[%d] = %f \n", i, rBuf[i]);
		}
	}
/*
	if(ProcID == receiver){
		for( i = 0; i < l; i++){
            printf("rBuf[%d] = %f \n", i, rBuf[i]);
		}
	}	
	*/
}

template<typename T,typename S>
void parMatrixSparse<T,S>::CSR_MatVecProd(parVector<T,S> *XVec, parVector<T,S> *YVec){
	S	i,j,k,l;
	S	llength, glength;
	S	count; count = 0;
	S	count2; count2 = 0;
	S	tag1; tag1 = 0;
	T	v; v = 0;

	MPI_Request *Rreqs, *Sreqs;
	MPI_Status  *Rstat, *Sstat;

	T	*rBuf, *sBuf;

	//get local and global length
	llength = XVec->GetLocalSize();
	glength = XVec->GetGlobalSize();

	//setting recv and send buffers
	rBuf = new T [glength];

	for (i = 0; i < glength; i++){
		rBuf[i] = 0;
	}

	sBuf = XVec->GetArray();

    Rreqs = new MPI_Request [nProcs - 1];
    Sreqs = new MPI_Request [nProcs - 1];
    Rstat = new MPI_Status [nProcs - 1];
    Sstat = new MPI_Status [nProcs - 1];

	count = 0;
	for(i = 0; i < nProcs; i++){
		if(i != ProcID){
			if(DTypeSend[i] != MPI_DATATYPE_NULL){
				MPI_Isend(sBuf,1,DTypeSend[i],i,tag1,MPI_COMM_WORLD,&Sreqs[count]);
			}
			else{
				Sreqs[count] = MPI_REQUEST_NULL;
			}
			if(DTypeRecv[i] != MPI_DATATYPE_NULL){
				MPI_Irecv(rBuf,1,DTypeRecv[i],i,tag1,MPI_COMM_WORLD,&Rreqs[count]);
			}
			else{
				Rreqs[count] = MPI_REQUEST_NULL;
			}
			count++;
		}
	}


#ifndef _OPENMP
	//calculate local-local product
	if(CSR_lloc != NULL){
		for(i = 0; i < nrows; i++){
			v = 0.0;
			for(k = CSR_lloc->rows[i]; k < CSR_lloc->rows[i+1]; k++){
				j = x_index_map->Glob2Loc(CSR_lloc->cols[k]);
				v += CSR_lloc->vals[k]*sBuf[j];
			}	
			YVec->AddValueLocal(i,v);
//				printf("v[%d] = %f\n", i,v);
			
		}
	}

	MPI_Waitall(nProcs-1, Rreqs, Rstat);

/*
	for (i = 0; i < glength; i++){
		printf("rBuf[%d] = %f\n", i, rBuf[i]);
	}
*/
	//calculate local-global product
	if(CSR_gloc != NULL){
		for(i = 0; i < nrows; i++){
			v = 0;
			for(k = CSR_gloc->rows[i]; k < CSR_gloc->rows[i+1];k++){
				j = CSR_gloc->cols[k];
				v += CSR_gloc->vals[k]*rBuf[j];
//				printf("CSR_gloc->cols[%d] = %d\n", k, CSR_gloc->cols[k]);
//				printf("i = %d, v = %f\n",i, v);
			}	
			YVec->AddValueLocal(i,v);
//				printf("vv[%d] = %f\n", i, v);
			
		}
	}

#else

	#pragma omp parallel default (shared) private (i,j,k)
	{
	        if(CSR_lloc != NULL){
        	#pragma omp for schedule(static)
		#pragma vector aligned
		       	for(i = 0; i < nrows; i++){
				v = 0;
//			#pragma omp parallel for reduction(+:v)
                        	for(k = CSR_lloc->rows[i]; k < CSR_lloc->rows[i+1]; k++){
                                	j = x_index_map->Glob2Loc(CSR_lloc->cols[k]);
                                	v += CSR_lloc->vals[k]*sBuf[j];
                         
                        	}
				YVec->AddValueLocal(i,v);
                	}
        	}
	}

	#pragma omp barrier
	#pragma omp master
	{
		MPI_Waitall(nProcs-1, Rreqs, Rstat);
	}

	#pragma omp barrier

	#pragma omp parallel default (shared) private (i,j,k)
        {
	if(CSR_gloc != NULL){
        #pragma omp for schedule(static) 
	#pragma vector aligned
	        for(i = 0; i < nrows; i++){
			v = 0;
//			#pragma omp parallel for reduction(+:v)
                        for(k = CSR_gloc->rows[i]; k < CSR_gloc->rows[i+1];k++){
                                j = CSR_gloc->cols[k];
                                v += CSR_gloc->vals[k]*rBuf[j];
                	}        
		        YVec->AddValueLocal(i,v);
                        
                }
        }
	
	}

#endif

	delete [] rBuf;
	delete [] Rreqs;
	delete [] Sreqs;
	delete [] Rstat;
	delete [] Sstat;
}	


template<typename T,typename S>
void parMatrixSparse<T,S>::MatAXPY(parMatrixSparse<T,S> *X, T scale){

	typename std::map<S,T>::iterator it, itv, itvv;

	S i, k;
	T v;
	if(dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}

	if(dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}

	if(dynmat_lloc != NULL && X->dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
			merge.insert(X->dynmat_lloc[i].begin(),X->dynmat_lloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_lloc[i][k] = dynmat_lloc[i][k]+X->dynmat_lloc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_lloc == NULL && X->dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_lloc[i].begin(),X->dynmat_lloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_lloc[i][k] = X->dynmat_lloc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_gloc != NULL && X->dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());
			merge.insert(X->dynmat_gloc[i].begin(),X->dynmat_gloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_gloc[i][k] =dynmat_gloc[i][k]+ X->dynmat_gloc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_gloc == NULL && X->dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_gloc[i].begin(),X->dynmat_gloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_gloc[i][k] = X->dynmat_gloc[i][k];
			}
			merge.clear();
		}
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::MatScale(T scale){
	typename std::map<S,T>::iterator it;

	S i;
	T v;

	if(dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}

	if(dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::MatAYPX(parMatrixSparse<T,S> *X, T scale){

	typename std::map<S,T>::iterator it, itv, itvv;

	S i, k;
	T v;
	if(X->dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = X->dynmat_lloc[i].begin(); it != X->dynmat_lloc[i].end(); it++){
				k = it->first;
				X->dynmat_lloc[i][k] = X->dynmat_lloc[i][k]*scale;
			}
		}
	}

	if(X->dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = X->dynmat_gloc[i].begin(); it != X->dynmat_gloc[i].end(); it++){
				k = it->first;
				X->dynmat_gloc[i][k] =X->dynmat_gloc[i][k]*scale;
			}
		}
	}

	if(dynmat_lloc != NULL && X->dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_lloc[i].begin(),dynmat_lloc[i].end());
			merge.insert(X->dynmat_lloc[i].begin(),X->dynmat_lloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_lloc[i][k] = dynmat_lloc[i][k]+X->dynmat_lloc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_lloc == NULL && X->dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_lloc[i].begin(),X->dynmat_lloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_lloc[i][k] = X->dynmat_lloc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_gloc != NULL && X->dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_gloc[i].begin(),dynmat_gloc[i].end());
			merge.insert(X->dynmat_gloc[i].begin(),X->dynmat_gloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_gloc[i][k] =dynmat_gloc[i][k]+ X->dynmat_gloc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_gloc == NULL && X->dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_gloc[i].begin(),X->dynmat_gloc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_gloc[i][k] = X->dynmat_gloc[i][k];
			}
			merge.clear();
		}
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_MatAXPY(parMatrixSparse<T,S> *X, T scale){

	typename std::map<S,T>::iterator it, itv, itvv;

	S i, k;
	T v;
	if(dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}

	if(dynmat_loc != NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_loc[i].begin(),dynmat_loc[i].end());
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = dynmat_loc[i][k]+X->dynmat_loc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_loc == NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = X->dynmat_loc[i][k];
			}
			merge.clear();
		}
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_MatScale(T scale){
	typename std::map<S,T>::iterator it;

	S i;
	T v;

	if(dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
				it->second = it->second*scale;
			}
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_MatAYPX(parMatrixSparse<T,S> *X, T scale){

	typename std::map<S,T>::iterator it, itv, itvv;

	S i, k;
	T v;
	if(X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = X->dynmat_loc[i].begin(); it != X->dynmat_loc[i].end(); it++){
				k = it->first;
				X->dynmat_loc[i][k] = X->dynmat_loc[i][k]*scale;
			}
		}
	}

	if(dynmat_loc != NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(dynmat_loc[i].begin(),dynmat_loc[i].end());
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = dynmat_loc[i][k]+X->dynmat_loc[i][k];
			}
			merge.clear();
		}
	}

	if(dynmat_loc == NULL && X->dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			std::map<S,T> merge;
			merge.insert(X->dynmat_loc[i].begin(),X->dynmat_loc[i].end());
			for(it = merge.begin(); it != merge.end(); ++it){
				k = it->first;
				dynmat_loc[i][k] = X->dynmat_loc[i][k];
			}
			merge.clear();
		}
	}
}



template<typename T,typename S>
void parMatrixSparse<T,S>::ZeroEntries()
{
	typename std::map<S,T>::iterator it;

	S i;
	T v;

	if(dynmat_lloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_lloc[i].begin(); it != dynmat_lloc[i].end(); it++){
				it->second = 0;
			}
		}
	}

	if(dynmat_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_gloc[i].begin(); it != dynmat_gloc[i].end(); it++){
				it->second = 0;
			}
		}
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::Loc_ZeroEntries()
{
	typename std::map<S,T>::iterator it;

	S i;
	T v;

	if(dynmat_loc != NULL){
		for(i = 0; i < nrows; i++){
			for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); it++){
				it->second = 0;
			}
		}
	}
}


//matrix multiple a special nilpotent matrix
template<typename T,typename S>
void parMatrixSparse<T,S>::MA(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod)
{
	S i, j, k;
	T v;

	typename std::map<S,T>::iterator it;

	if(nilp.setup == true && y_index_map->GetGlobalSize() == nilp.matrix_size){
		if(ProcID == 0){
			std::cout << "INFO ]>: Nilpotent matrix diagPostion = " << nilp.diagPosition << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix nbOne = " << nilp.nbOne << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix size = " << nilp.matrix_size << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix setup = " << nilp.setup << std::endl;
		}
	}
	else if(nilp.setup == false){
		if(ProcID == 0){
			std::cout << "ERROR ]>: Nilpotent Matrix is not setup." << std::endl;
		}
		return;
	}
	else if(y_index_map->GetGlobalSize() != nilp.matrix_size){
		if(ProcID == 0){
			std::cout << "ERROR ]>: Nilpotent Matrix dimension do not equal to given matrix size." << std::endl;
		}
		return;
	}

	//use the given nilpotency matrix, MA operation will make elements of matrix right move diaPosition-1 offset.
	//And the positions of 0: pos = nbOne*integer - 1

	for(i = 0; i < nrows; i++){
		if(dynmat_loc == NULL) {
			return;
		}
		if(dynmat_loc != NULL && prod->dynmat_loc == NULL){
			prod->dynmat_loc = new std::map<S,T> [nrows];
		}
		
		for(it = dynmat_loc[i].begin(); it != dynmat_loc[i].end(); ++it){
			j = it->first + nilp.diagPosition - 1;
			k = (j+1)%(nilp.nbOne + 1);
			if(j < ncols && k != 0){
				prod->dynmat_loc[i][j] = it->second;
			}
		}
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::AM_SetUpDataTypes(Nilpotency<S> nilp){
	
	if(nilp.setup == true && y_index_map->GetGlobalSize() == nilp.matrix_size){
		if(ProcID == 0){
			std::cout << "INFO ]>: Nilpotent matrix diagPostion = " << nilp.diagPosition << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix nbOne = " << nilp.nbOne << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix size = " << nilp.matrix_size << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix setup = " << nilp.setup << std::endl;
		}
	}
	else if(nilp.setup == false){
		if(ProcID == 0){
			std::cout << "ERROR ]>: Nilpotent Matrix is not setup." << std::endl;
		}
		return;
	}
	else if(y_index_map->GetGlobalSize() != nilp.matrix_size){
		if(ProcID == 0){
			std::cout << "ERROR ]>: Nilpotent Matrix dimension do not equal to given matrix size." << std::endl;
		}
		return;
	}
}
//special nilpotent matrix multiple another matrix
template<typename T,typename S>
void parMatrixSparse<T,S>::AM(Nilpotency<S> nilp, parMatrixSparse<T,S> *prod)
{
	S i, j, k, p, q, c, loc;
	T v;

	typename std::map<S,T>::iterator it;

	if(nilp.setup == true && y_index_map->GetGlobalSize() == nilp.matrix_size){
		if(ProcID == 0){
			std::cout << "INFO ]>: Nilpotent matrix diagPostion = " << nilp.diagPosition << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix nbOne = " << nilp.nbOne << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix size = " << nilp.matrix_size << std::endl;
			std::cout << "INFO ]>: Nilpotent matrix setup = " << nilp.setup << std::endl;
		}
	}
	else if(nilp.setup == false){
		if(ProcID == 0){
			std::cout << "ERROR ]>: Nilpotent Matrix is not setup." << std::endl;
		}
		return;
	}
	else if(y_index_map->GetGlobalSize() != nilp.matrix_size){
		if(ProcID == 0){
			std::cout << "ERROR ]>: Nilpotent Matrix dimension do not equal to given matrix size." << std::endl;
		}
		return;
	}

	//use the given nilpotency matrix, AM operation will make elements of matrix up move diaPosition-1 offset.
	//And the positions of 0: pos = nbOne*integer - 1

	if((nilp.diagPosition - 1) > nrows){
		if(ProcID == 0){
			std::cout << "INFO ]>: Nilpotent matrix diagPostion - 1 is larger than local rows number," << std::endl << "         sending should be divided into block"  << std::endl;
		}
	}
	else{
		if(ProcID == 0){
			std::cout << "INFO ]>: The MPI sending and receiving row num = " << nilp.diagPosition - 1 << std::endl;
		}
	}

	MPI_Request	*Rreqs, *Sreqs, rtypereq, stypereq, *idxRreqs, *idxSreqs;

	MPI_Status	*status, *Rstat, *Sstat, typestat, *idxRstat, *idxSstat;

	int up, down;

	int			tag1[nilp.diagPosition - 1], tag2[nilp.diagPosition - 1], tagtype = nilp.diagPosition - 1 + nilp.diagPosition ;
	int			Rtag, Stag;
	S           count, count1;

	count = 0;
	count1 = 0;

	up = ProcID - 1; 
	down =ProcID + 1;
	S size[nilp.diagPosition - 1], rsize[nilp.diagPosition - 1];

	for(S a = 0; a < nilp.diagPosition - 1; a++){
		tag1[a] = a;
		tag2[a] = nilp.diagPosition - 1 + a;
	}

	if(ProcID != 0){
		for(S a = 0; a < nilp.diagPosition - 1; a++){
			size[a] = dynmat_loc[a].size();
		}
	}

	if(ProcID != nProcs - 1){
		for(S a = 0; a < nilp.diagPosition - 1; a++){
			rsize[a] = 0;
		}
	}



	if(ProcID != 0){
		MPI_Isend(size, nilp.diagPosition - 1, MPI_INT, up, tagtype, MPI_COMM_WORLD, &stypereq);
	}

	if(ProcID != nProcs - 1){
		MPI_Irecv(rsize,nilp.diagPosition - 1, MPI_INT, down, tagtype, MPI_COMM_WORLD, &rtypereq);
	}
	
	if(ProcID != nProcs - 1){
		MPI_Wait(&rtypereq,&typestat);
	}

/*
	if(ProcID != nProcs - 1){
		std::cout << "s-s-s--s-s-s-s-s-s-s----->" << rsize[0] << std::endl;
	}

*/
	Rtag = 0; Stag = 1;

	//MPI non-blocking requests and status
	
	Rreqs = new MPI_Request [nilp.diagPosition - 1];
	Sreqs = new MPI_Request [nilp.diagPosition - 1];
	Rstat = new MPI_Status [nilp.diagPosition - 1];
	Sstat = new MPI_Status [nilp.diagPosition - 1];

	idxRreqs = new MPI_Request [nilp.diagPosition - 1];
	idxSreqs = new MPI_Request [nilp.diagPosition - 1];
	idxRstat = new MPI_Status [nilp.diagPosition - 1];
	idxSstat = new MPI_Status [nilp.diagPosition - 1];

	if(dynmat_loc == NULL){
		return;
	}

	if(std::is_same<T,std::complex<double> >::value){
		if(ProcID == 0){
			std::cout << "INFO ]>: Using Complex Double values" << std::endl; 
		}
	}

	if(std::is_same<T,std::complex<float> >::value){
		if(ProcID == 0){
			std::cout << "INFO ]>: Using Complex Single values" << std::endl; 
		}
	}

	for(p = nilp.diagPosition - 1; p < nrows; p++){
		if(prod->dynmat_loc == NULL){
			prod->dynmat_loc = new std::map<S,T> [nrows];
		}

		i = p - nilp.diagPosition + 1;
		q = y_index_map->Loc2Glob(i);

		k = (q + 1)%(nilp.nbOne + 1);
		
		for(it = dynmat_loc[p].begin(); it != dynmat_loc[p].end(); ++it){
			
			j = it->first;

			if(i >= 0 && i < nrows && k != 0){
				prod->dynmat_loc[i][j] = it->second;
				//std::cout << "ProcID = " << ProcID << ": i = " << i << ", j = " << j << std::endl;
			}
		}
		
	}

	MPI_Datatype MPI_SCALAR = MPI_Scalar<T>();

	MPI_Barrier(MPI_COMM_WORLD);

#if defined (__USE_COMPLEX__) && defined(__USE_DOUBLE__)
//complex double
	double *sBuf, *rBuf;

	S *sIndx, *rIndx;
	S cnt = 0, cnt2 = 0;

	for(S b = 0; b < nilp.diagPosition -1; b++){

		if(ProcID != 0){
			loc = y_index_map->Loc2Glob(b);
			q = loc - nilp.diagPosition + 1 + b;
//					std::cout << "----->loc = " << loc << " divisblae = "<< (q + 1)%(nilp.nbOne + 1) << std::endl;
			if(q == 0){
				return;
			}
		}
		
		if(ProcID != 0){
			sBuf  = new double [2*size[b]];
			sIndx  = new S [size[b]];
			for(it = dynmat_loc[b].begin(); it != dynmat_loc[b].end(); ++it){
				sBuf[cnt] = (it->second).real();
				sBuf[cnt+1] = (it->second).imag();
				sIndx[cnt2] = it->first;
				cnt = cnt + 2;
				cnt2++;
			}

			MPI_Isend(sBuf, size[b], MPI_SCALAR, up, tag1[b], MPI_COMM_WORLD, &Sreqs[b]);
			MPI_Isend(sIndx, size[b], MPI_INT, up, tag2[b], MPI_COMM_WORLD, &idxSreqs[b]);
		}

		if(ProcID != nProcs - 1){
			rBuf = new double [2*rsize[b]];
			rIndx = new S [2*rsize[b]];
			for(S tt =0; tt < 2*rsize[b]; tt++){
				rBuf[tt] = 0.0;
			}

			for(S tt =0; tt < rsize[b]; tt++){
				rIndx[tt] = 0;
			}

			MPI_Irecv(rBuf,rsize[b], MPI_SCALAR, down, tag1[b], MPI_COMM_WORLD, &Rreqs[b]);	
			MPI_Irecv(rIndx,rsize[b], MPI_INT, down, tag2[b], MPI_COMM_WORLD, &idxRreqs[b]);	

		}

		if(ProcID != nProcs - 1){
			MPI_Wait(&Rreqs[b],&Rstat[b]);
			MPI_Wait(&idxRreqs[b],&idxRstat[b]);
		}

		if(ProcID != nProcs - 1){
				for(S tt = 0; tt < rsize[b]; tt++){
						//std::cout << " ===> rBuf[" << tt << "] = " << rBuf[tt] << "   rIndx[" << tt << "] = " << rIndx[tt] << std::endl;
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).real(rBuf[2*tt]);
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).imag(rBuf[2*tt + 1]);
				}
		
		}

	}

#elif defined (__USE_COMPLEX__)
//complex single

	float *sBuf, *rBuf;

	S *sIndx, *rIndx;
	S cnt = 0, cnt2 = 0;

	for(S b = 0; b < nilp.diagPosition -1; b++){

		if(ProcID != 0){
			loc = y_index_map->Loc2Glob(b);
			q = loc - nilp.diagPosition + 1 + b;
//					std::cout << "----->loc = " << loc << " divisblae = "<< (q + 1)%(nilp.nbOne + 1) << std::endl;
			if(q == 0){
				return;
			}
		}
		
		if(ProcID != 0){
			sBuf  = new float [2*size[b]];
			sIndx  = new S [size[b]];
			for(it = dynmat_loc[b].begin(); it != dynmat_loc[b].end(); ++it){
				sBuf[cnt] = (it->second).real();
				sBuf[cnt+1] = (it->second).imag();
				sIndx[cnt2] = it->first;
				cnt = cnt + 2;
				cnt2++;
			}

			MPI_Isend(sBuf, size[b], MPI_SCALAR, up, tag1[b], MPI_COMM_WORLD, &Sreqs[b]);
			MPI_Isend(sIndx, size[b], MPI_INT, up, tag2[b], MPI_COMM_WORLD, &idxSreqs[b]);
		}

		if(ProcID != nProcs - 1){
			rBuf = new float [2*rsize[b]];
			rIndx = new S [2*rsize[b]];
			for(S tt =0; tt < 2*rsize[b]; tt++){
				rBuf[tt] = 0.0;
			}

			for(S tt =0; tt < rsize[b]; tt++){
				rIndx[tt] = 0;
			}

			MPI_Irecv(rBuf,rsize[b], MPI_SCALAR, down, tag1[b], MPI_COMM_WORLD, &Rreqs[b]);	
			MPI_Irecv(rIndx,rsize[b], MPI_INT, down, tag2[b], MPI_COMM_WORLD, &idxRreqs[b]);	

		}

		if(ProcID != nProcs - 1){
			MPI_Wait(&Rreqs[b],&Rstat[b]);
			MPI_Wait(&idxRreqs[b],&idxRstat[b]);
		}

		if(ProcID != nProcs - 1){
				for(S tt = 0; tt < rsize[b]; tt++){
						//std::cout << " ===> rBuf[" << tt << "] = " << rBuf[tt] << "   rIndx[" << tt << "] = " << rIndx[tt] << std::endl;
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).real(rBuf[2*tt]);
						(prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]]).imag(rBuf[2*tt + 1]);
				}
		
		}

	}
}

#else
//real

	T *sBuf, *rBuf;
	S *sIndx, *rIndx;
	S cnt = 0, cnt2 = 0;

	for(S b = 0; b < nilp.diagPosition -1; b++){

		if(ProcID != 0){
			loc = y_index_map->Loc2Glob(b);
			q = loc - nilp.diagPosition + 1 + b;
//					std::cout << "----->loc = " << loc << " divisblae = "<< (q + 1)%(nilp.nbOne + 1) << std::endl;
			if(q == 0){
				return;
			}
		}
		
		if(ProcID != 0){
			sBuf  = new T [size[b]];
			sIndx  = new S [size[b]];
			for(it = dynmat_loc[b].begin(); it != dynmat_loc[b].end(); ++it){
				sBuf[cnt] = it->second;
				sIndx[cnt] = it->first;
				cnt++;
			}
			MPI_Isend(sBuf, size[b], MPI_SCALAR, up, tag1[b], MPI_COMM_WORLD, &Sreqs[b]);
			MPI_Isend(sIndx, size[b], MPI_INT, up, tag2[b], MPI_COMM_WORLD, &idxSreqs[b]);
		}

		if(ProcID != nProcs - 1){
			rBuf = new T [rsize[b]];
			rIndx = new S [rsize[b]];
			for(S tt =0; tt < rsize[b]; tt++){
				rBuf[tt] = 0.0;
				rIndx[tt] = 0;
			}
			MPI_Irecv(rBuf,rsize[b], MPI_SCALAR, down, tag1[b], MPI_COMM_WORLD, &Rreqs[b]);	
			MPI_Irecv(rIndx,rsize[b], MPI_INT, down, tag2[b], MPI_COMM_WORLD, &idxRreqs[b]);	

		}

		if(ProcID != nProcs - 1){
			MPI_Wait(&Rreqs[b],&Rstat[b]);
			MPI_Wait(&idxRreqs[b],&idxRstat[b]);
		}

		if(ProcID != nProcs - 1){
				for(S tt = 0; tt < rsize[b]; tt++){
						//std::cout << " ===> rBuf[" << tt << "] = " << rBuf[tt] << "   rIndx[" << tt << "] = " << rIndx[tt] << std::endl;
						prod->dynmat_loc[nrows - nilp.diagPosition + 1 + b][rIndx[tt]] = rBuf[tt];
				}
		
		}

	}

#endif
	
}
