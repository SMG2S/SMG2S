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

template<typename T,typename S>
void parMatrixSparse<T,S>::ConvertToCSR()
{
	S 	count, i, j, k;
	T	v;
	std::map<S,T>::iterator it;
	
	if(dynmat_lloc != NULL){
		//allocate csr matrix
		CSR_lloc = new MatrixCSR(nnz_lloc, nrows);

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
		CSR_gloc = new MatrixCSR(nnz_gloc, nrows);
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
}

template<typename T,typename S>
void parMatrixSparse<T,S>::FindColsToRecv()
{
	std::map<S,S> Rrows;
	std::map<S,S> Srows;
	std::map<S,S>::iterator vit;
	std::map<S,T>::iterator mit;

	S	i, j, k;
	S	count, count1;

	count = 0;
	count1 = 0;

	S	nRecv;

	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &nRrocs);

	//Initialise vector containing numer of entries to send and recv on each procss
	VNumRecv = new S [nProcs];
	for(i = 0; i < nProcs; i++){
		VNumRecv[i] = 0;
	}
	VNumRecv = new S [nProcs];
	for(i = 0; i < nProcs; i++){
		VNumSend[i] = 0;
	}

	MPI_Request	*Rreqs, *Sreqs;
	MPI_Status	statusm *Rstat, *Sstat;
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
					count++:
				}
			}
		}
		nRecv = count;
	}
	else{
		for(i = 0; i < nProcs; i++){
			VNumRecv[i] = 0;
		}
	}

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
			for(vit = Rows.begin(); vit  != Rrows.end(); vit++){
				if(vit->second == i){
					Sbuffer[i][count] = vit->first;
					count++;
				}
			}
			MPI_Isend(Sbuffer[i],VNumRecv[i], MPI_INT,i,tag1,MPI_COMM_WORLD, &Sreqs[count1]);
			count++:
		}
		else{
			Sbuffer[i] = new S [1];
			Sbuffer[i][0] = 0;
		}
	}

	count1 = 1;
	
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

	delete [] Rreqs;	
	delete [] Sreqs;
	delete [] Rstat;
	delete [] Sstat;

}

