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
parMatrixSparse<T,S>::parMatrixSparse(parVector<T,S> *XVec, parVector<T,S> *YVec)
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
void parMatrixSparse<T,S>::AddValuesLocal(S nindex, S *rows, S *cols, T *values)
{
	typename std::map<S,T>::iterator it;
	
	for( S i = 0; i < nindex; i++){
		AddValueLocal(rows[i],cols[i],values[i]);
	}
}

template<typename T,typename S>
T parMatrixSparse<T,S>::GetValue(S row, S col)
{
	if((row < nrows && row >= 0) && (col < upper_x && col >= lower_x && col > 0)){
		return dynmat_lloc[row][col];
	}
	else if ((row < nrows && row >= 0) && (col >= upper_x || col < lower_x) && (col >= 0)){
		return dynmat_gloc[row][col];
	}
}


template<typename T,typename S>
void parMatrixSparse<T,S>::MatView(){
	S i, j;
	T v;

	if(ProcID == 0) {std::cout << "Parallel MatView: " << std::endl;}

	for(i = 0; i < nrows; i++){
		std::cout << "row " << i << ":";
		for (j = 0; j < ncols; j++){
			if(j < upper_x && j >= lower_x){
				if(dynmat_lloc[i][j] != 0.0){
					std::cout << ", " << dynmat_lloc[i][j] << " "<< j;
				}
			}
			else if(j >= upper_x || j < lower_x){
				std::cout << ", " << j;
			}
		}


/*
		for (j = 0; j < ncols; j++){
			if(j < upper_x && j >= lower_x){
				std::cout << " "<< "("<< j << ","<< dynmat_lloc[i+1][j+1] << ")";
			}
			else if(j >= upper_x || j < lower_x){
				std::cout << " "<< "("<< j << "," << dynmat_gloc[i+1][j+1] << ")";
			}
		}
*/
	std::cout << std::endl;
	}
}

template<typename T,typename S>
void parMatrixSparse<T,S>::SetDiagonal(parVector<T,S> *diag)
{
	if (nrows != njloc ){
		if(ProcID == 0){
			printf("ERROR: cannot set digonal for non-square matrix.");
		}
	}
	else{
		T *array = diag->GetArray();

		for(S i = 0; i < nrows; i++){
			AddValueLocal(i,i,array[i]);
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
	typename std::map<S,S>::iterator vit;
	typename std::map<S,T>::iterator mit;

	S	i, j, k;
	S	count, count1;

	count = 0;
	count1 = 0;

	S	nRecv;

	MPI_Comm_rank(MPI_COMM_WORLD, &ProcID);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	//Initialise vector containing numer of entries to send and recv on each procss
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
					count++;
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
			for(vit = Rrows.begin(); vit  != Rrows.end(); vit++){
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


template<typename T,typename S>
void parMatrixSparse<T,S>::SetupDataTypes()
{
	S	i,j,k;
	S	count, *blength, *displac;

	DTypeSend = new MPI_Datatype [nProcs];
	DTypeRecv = new MPI_Datatype [nProcs];

	for(i = 0; i < nProcs; i++){
		count = VNumSend[i];
		blength = new S [count];
		displac = new S [count];

		for(j = 0; j < count; j++){
			blength[j] = 1;
			displac[j] = x_index_map->Glob2Loc(Rbuffer[i][j]);
		}

		MPI_Type_indexed(count, blength, displac, MPI_DOUBLE, &DTypeSend[i]);
		MPI_Type_commit(&DTypeSend[i]);

		count = VNumRecv[i];
		blength = new S [count];
		displac = new S [count];

		for(j = 0; j < count; j++){
			blength[j] = 1;
			displac[j] = Sbuffer[i][j];
		}

		MPI_Type_indexed(count, blength, displac, MPI_DOUBLE, &DTypeRecv[i]);
                MPI_Type_commit(&DTypeRecv[i]);
	}
	
	delete [] blength;
	delete [] displac;
}

template<typename T,typename S>
void parMatrixSparse<T,S>::TestCommunication(parVector<T,S> *XVec, parVector<T,S> *YVec)
{
	S	i,j,k,l,ng;
	S	sender, receiver;
	sender = 2;
	receiver = 1;
	
	MPI_Status Rstat;

	T	*rBuf, *sBuf;

	k = XVec->GetLocalSize();
	l = XVec->GetGlobalSize();
	rBuf = new T [1];
	for(i = 0; i < 1; i++){
		rBuf[i] = 0;	
	}
	
	sBuf = XVec->GetArray();

	if(ProcID == sender){
		MPI_Send(sBuf, 1, DTypeSend[receiver], receiver, 1, MPI_COMM_WORLD);
		printf("sending done! \n");
	}
	else{
		MPI_Recv(rBuf,1,DTypeRecv[sender],sender,1,MPI_COMM_WORLD,&Rstat);
		printf("receiving done! \n");
	}

	if(ProcID == sender){
		for( i = 0; i < k; i++){
			printf("sBuf[%d] = %f \n", i, sBuf[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(ProcID == receiver){
		for( i = 0; i < l; i++){
                        printf("rBuf[%d] = %f \n", i, rBuf[i]);
		}
	}	
}

template<typename T,typename S>
void parMatrixSparse<T,S>::MatVecProd(parVector<T,S> *XVec, parVector<T,S> *YVec){
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
			for(k = CSR_lloc->rows[i]; k < CSR_lloc->rows[i+1]; k++){
				j = x_index_map->Glob2Loc(CSR_lloc->cols[k]);
				v = CSR_lloc->vals[k]*sBuf[j];
				YVec->AddValueLocal(i,v);
			}
		}
	}

	MPI_Waitall(nProcs-1, Rreqs, Rstat);

	//calculate local-global product
	if(CSR_gloc != NULL){
		for(i = 0; i < nrows; i++){
			for(k = CSR_gloc->rows[i]; k < CSR_gloc->rows[i+1];k++){
				j = CSR_gloc->cols[k];
				v = CSR_gloc->vals[k]*rBuf[j];
				YVec->AddValueLocal(i,v);
			}
		}
	}
#else

	#pragma omp parallel default (shared) private (i,j,k,v)
	{
	        if(CSR_lloc != NULL){
        	#pragma omp for schedule(guided)
		       	for(i = 0; i < nrows; i++){
                        	for(k = CSR_lloc->rows[i]; k < CSR_lloc->rows[i+1]; k++){
                                	j = x_index_map->Glob2Loc(CSR_lloc->cols[k]);
                                	v = CSR_lloc->vals[k]*sBuf[j];
                                	YVec->AddValueLocal(i,v);
                        	}
                	}
        	}
	}

	#pragma omp barrier
	#pragma omp master
	{
		MPI_Waitall(nProcs-1, Rreqs, Rstat);
	}

	#pragma omp barrier
        if(CSR_gloc != NULL){
        #pragma omp for schedule(guided)
	        for(i = 0; i < nrows; i++){
                        for(k = CSR_gloc->rows[i]; k < CSR_gloc->rows[i+1];k++){
                                j = CSR_gloc->cols[k];
                                v = CSR_gloc->vals[k]*rBuf[j];
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




