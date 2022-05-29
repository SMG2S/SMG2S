/*
MIT License
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France
Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
									 Materials,  Forschungszentrum Juelich GmbH.
									 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef __PAR_VECTOR_H__
#define __PAR_VECTOR_H__

#include <parVector/parVectorMap.hpp>
#include <utils/utils.hpp>
#include <utils/MPI_DataType.hpp>

//!  @brief A class which defines a vector distributed across 1D MPI grid.
/*!
 * @ingroup group6 
  - This class can be constructed with user-provided of index ranges `[lower_bound, upper_bound)` of vector on each MPI proc.
  - This class can be constructed with a distribution scheme by a given parVectorMap object.

  @tparam T describes the scalar types of the entries of vector.   
  @tparam S type of integer to describes the dimension of vector to be generated. 
*/
template<typename T, typename S>
class parVector{
    private:
    //! An array store the local piece of a global vector on each MPI proc
    T	*array;
	//! The number of elements of vector stored on each MPI proc    
	S	local_size;
	//! Global size of this distributed vector	
	S       global_size;
    //! The working MPI Communicator	
	MPI_Comm  	comm;	
	//! A parVectorMap object which shows the distribution scheme of a vector
	parVectorMap<S> index_map;
	//! rank of each MPI procs within the working MPI communicator parVectorMap#comm	
	int 	MyPID;
	//! number of MPI procs within the working MPI communicator parVectorMap#comm
	int 	nProcs; 

    public:
	parVector();
    //! A constructor of `parVector`. 
    /*!
      * @param[in] ncomm the working MPI Communicator
      * @param[in] lbound the smallest index of a distributed vector on each MPI proc
      * @param[in] ubound `ubound-1 = `  the largest index of a distributed vector on each MPI proc 
      
      - parVector::array is also allocated. 
    */		
	parVector(MPI_Comm ncomm, S lbound, S ubound);
    //! A constructor of `parVectorMap`. 
    /*!
      * @param[in] map the distribution scheme determined by this object of type parVectorMap

      - parVector::array is also allocated. 
    */		
	parVector(parVectorMap<S> map);
	//! A destructor
	~parVector();

	//! Return parVector#index_map
	parVectorMap<S> GetVecMap(){return index_map;};

	//! Return the lower_bound on  each MPI proc 
	S GetLowerBound(){return index_map.GetLowerBound();};
	//! Return the upper_bound on  each MPI proc 
	S GetUpperBound(){return index_map.GetUpperBound();};
	//! Return parVector<S>#global_size			
	S GetGlobalSize(){return global_size;};
	//! Return parVector<S>#local_size		
	S GetLocalSize(){return local_size;};
	//! Return parVector<S>#MyPID
	S GetRank(){return index_map.GetRank();};

	//! For each MPI of rank `i`, get a value from the local vector stored on rank `i-1`.
    /*!
      * @param[in] offset get the value from local vector on rank `i-1` with index `upper_bound-1-offset`.

      - Attention, for the MPI rank `0`, the value is obtained from the MPI proc of rank `nProcs-1`. 
    */		
	T GetUpperNeighbor(S offset);
	//! For each MPI of rank `i`, get a value from the local vector stored on rank `i+1`.
    /*!
      * @param[in] offset get the value from local vector on rank `i+1` with index `offset`.

      - Attention, for the MPI rank `nProcs-1`, the value is obtained from the MPI proc of rank `0`. 
    */	
	T GetLowerNeighbor(S offset);
	//! Get a value of vector with a given global index
    /*!
      * @param[in] index the querying global index

      - Attention, this function is naturally distributed, so each MPI proc can only query the value within its range `[lower_bound, upper_bound)`. 
    */		
	T GetValue(S index)
	{
	    auto lindex = index_map.Glob2Loc(index);
	    return array[lindex];
	};
	//! Get a value of vector with a given local index on each MPI proc.
    /*!
      * @param[in] lindex the querying local index
    */	
	T GetValueLocal(S lindex)
	{
	    return array[lindex];
	};

	//! Get the pointer `*array` which stores the local piece of vector on each MPI proc
	T* GetArray(){return array;};
	//! Get working MPI communicator
	MPI_Comm GetComm(){return comm;};
	//! Convert a index of local vector on each MPI proc into its index in the global distributed vector
    /*!
      * @param[in] local_index the index local vector to be converted

      - Attention, this function is in distributed manner, each MPI proc can only convert the local index of vector stored on itself.
    */	
	S Loc2Glob(S local_index){return index_map.Loc2Glob(local_index);};
	//! Convert a index of global vector into its index in the local vector on each MPI proc.
    /*!
      * @param[in] global_index the index global vector to be converted

      - Attention, this function is in distributed manner, each MPI proc can only convert the global index of vector in the range `[lower_bound, upper_bound)`.
    */		
	S Glob2Loc(S global_index){return index_map.Glob2Loc(global_index);};

	//! Set all the entries of a vector to a same value
    /*!
      * @param[in] value the value to be set
    */		
	void SetToValue(T value);
	//! Set all the entries of a vector to zero	
	void SetToZero();
	//! Set a value on a local index of piece of distributed vector on each MPI proc 
    /*!
      * @param[in] row the local index
      * @param[in] value the scalar to be set
    */		
	void SetValueLocal(S row, T value);
	//! Set multiple values with multiple local indices of piece of distributed vector on each MPI proc
    /*!
      * @param[in] nindex number of values to be set
      * @param[in] rows an pointer stores all the values to be set
      * @param[in] values an pointer stores all the indices to be set
    */	
	void SetValuesLocal(S nindex, S *rows, T *values);
	//! Set a value with a global index distributed vector 
    /*!
      * @param[in] index the global index
      * @param[in] value the scalar to be set
    */		
	void SetValueGlobal(S index, T value);
	//! Set multiple values with multiple global indices of distributed vector
    /*!
      * @param[in] nindex number of values to be set
      * @param[in] rows an pointer stores all the values to be set
      * @param[in] values an pointer stores all the indices to be set
    */		
	void SetValuesGlobal(S nindex, S *rows, T *values);		
	//! add a value onto the value with a given global index of distributed vector 
    /*!
      * @param[in] row the global index
      * @param[in] value the scalar to be added
    */		
	void AddValueLocal(S row, T value);
	//! Add  multiple values with multiple global indices onto the related values of distributed vector
    /*!
      * @param[in] nindex number of values to be set
      * @param[in] rows an pointer stores all the values to be set
      * @param[in] values an pointer stores all the indices to be set
    */		
	void AddValuesLocal(S nindex, S *rows, T *values);
	//! Add with another vector
    /*!
      * @param[in] v the parVector object used to be added
    */		
	void VecAdd(parVector v);
	//! Scale all the elements of a vector with `scale`
    /*!
      * @param[in] scale the value used for scaling on the vector
    */		
	void VecScale(T scale);
	//! Compute the dot product with another vector `v`
    /*!
      * @param[in] v the parVector object used to perform a dot product

      - Attention, the results of dot product is stored rebundantly across all MPI procs 
    */		
	T    VecDot(parVector v);
	//! Read a vector from local text file.
    /*!
      * @param[in] spectrum a `std::string` indicates the path+filename of local text file
    */		
	void ReadExtVec(std::string spectrum);
    //! Write a parVector in real scalar to a local file
    /*!
      * @param[in] file_name a `std::string` indicates the path+filename to write into
    */      
    void writeToTxt(std::string file_name);
    //! Write a parVector in complex scalar to a local file
    /*!
      * @param[in] file_name a `std::string` indicates the path+filename to write into
    */        
    void writeToTxtCmplx(std::string file_name);

	//! Display the vector in a distributed manner		
	void VecView();
};

template<typename T,typename S>
parVector<T,S>::parVector(){}


template<typename T,typename S>
parVector<T,S>::parVector(MPI_Comm ncomm, S lbound, S ubound)
{
    MPI_Comm_dup(ncomm, &comm);
    index_map = parVectorMap<S>(ncomm, lbound, ubound);
    local_size = index_map.GetLocalSize();
    global_size = index_map.GetGlobalSize();
    array = new T[local_size];
    MPI_Comm_rank(ncomm, &MyPID);
    MPI_Comm_size(ncomm, &nProcs);
}

template<typename T,typename S>
parVector<T,S>::parVector(parVectorMap<S> map)
{
    index_map = map;
    comm = index_map.GetCurrentComm();
    local_size = index_map.GetLocalSize();
    global_size = index_map.GetGlobalSize();
    array = new T[local_size];
    MPI_Comm_rank(comm, &MyPID);
    MPI_Comm_size(comm, &nProcs);    
}

template<typename T,typename S>
parVector<T,S>::~parVector(){}


template<typename T, typename S>
void parVector<T,S>::SetToValue(T value)
{
    for(S i= 0; i < local_size; i++) {
    	array[i] = value;
    }
}

template<typename T, typename S>
void parVector<T,S>::SetToZero()
{
    T val = 0;
    SetToValue(val);
}

template<typename T, typename S>
void parVector<T,S>::VecView()
{
    int r;
    r = index_map.GetRank();
    if (r == 0){
	std::cout << "Parallel Vector View: " << std::endl << std::endl;
    }
    T *array = GetArray();
    S global;
    for(S i = 0; i < local_size; i++){
    	global = Loc2Glob(i);
	std::cout << "[" << global << "]: " << array[i] << std::endl;
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValueLocal(S row, T value)
{
    if (row < local_size){
	array[row] = value;
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValuesLocal(S nindex, S *rows, T *values)
{
    for(S i = 0; i < nindex; i++){
    	SetValueLocal(rows[i],values[i]);
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValueGlobal(S index, T value)
{
    int lower_bound = GetLowerBound();
    int upper_bound = GetUpperBound();
    if((index >= lower_bound) && (index < upper_bound)){
    	SetValueLocal(index_map.Glob2Loc(index), value);
    }
}

template<typename T, typename S>
void parVector<T,S>::SetValuesGlobal(S nindex, S *rows, T *values)
{
    for(S i = 0; i < nindex; i++){
	SetValueLocal(Glob2Loc(rows[i]),values[i]);
    }
}

template<typename T, typename S>
void parVector<T,S>::AddValueLocal(S row, T value)
{
    if (row < local_size){
	array[row] = array[row] + value;
    }
}

template<typename T, typename S>
void parVector<T,S>::AddValuesLocal(S nindex, S *rows, T *values)
{
    for(S i = 0; i < nindex; i++){
	AddValueLocal(rows[i],values[i]);
    }
}

template<typename T, typename S>
void parVector<T,S>::VecAdd(parVector<T,S> v)
{
    if(local_size != v.GetLocalSize()){
	std::cout << "vector size not coherant" << std::endl;
    }
    else{
	for(S i = 0; i < local_size; i++){
	    array[i] = array[i] + v.array[i];
	}
    }
}

template<typename T, typename S>
void parVector<T,S>::VecScale(T scale)
{
    for(S i = 0; i < local_size; i++){
	array[i] = scale*array[i];
    }
}

template<typename T, typename S>
T parVector<T,S>::VecDot(parVector v)
{
    T sum;
    for(S i = 0; i < local_size; i++){
	sum += array[i]*v.array[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, getMPI_Type<T>(), MPI_SUM, comm);
    return sum;
}


template<typename T, typename S>
void parVector<T, S>::ReadExtVec(std::string spectrum)
{
    std::ifstream file(spectrum);
    std::string line;

    int lower_bound = GetLowerBound();
    int upper_bound = GetUpperBound();

    int size1;
    size1 = sizeof(T) / sizeof(Base<T>);
    S idx;
    Base<T> in_vals[size1];
    T val;

    while (std::getline(file,line)) {
    	//std::cout << line << std::endl;
	idx = 0;
        for(int i = 0; i < size1; i++){
            in_vals[i] = 0.0;
        }

        std::stringstream linestream ( line ) ;
        linestream >> idx;

        for(int i = 0; i < size1; i++){
            linestream >> in_vals[i];
        }
        idx = idx - 1;

        for(int i = 0; i < size1; i++){
            reinterpret_cast<Base<T>(&)[size1]>(val)[i] = in_vals[i];
	}

	   if((idx >= lower_bound) && (idx < upper_bound)){
	       SetValueLocal(index_map.Glob2Loc(idx),val);
	   }		
    }
}

template<typename T, typename S>
void parVector<T, S>::writeToTxt(std::string file_name)
{
    std::string header;
    std::string data;

    if( (sizeof(T)/sizeof(Base<T>) != 1) ){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: for complex vector, please use writeToTxtCmplx" << std::endl;
        }               
    }

    //generate header
    header.append("%%SMG2S vector in real scalar\n");
    std::string dim_info = std::to_string(global_size) + " " + std::to_string(global_size) + "\n";
    header.append(dim_info);    

    for(auto i = 0; i < local_size; i++){
        data += std::to_string(index_map.Loc2Glob(i)+1) + " " + std::to_string(array[i]) + "\n";
    }

    MPI_File fh;
    MPI_Offset write_offset;
    MPI_Offset text_size;
    MPI_Offset *write_size_per_proc;
    MPI_Status sts;

    write_size_per_proc = (MPI_Offset *)malloc(sizeof(MPI_Offset) * nProcs);
    
    text_size = data.size();

    MPI_Allgather(&text_size, 1, MPI_OFFSET, write_size_per_proc, 1, MPI_OFFSET, comm);

    write_offset = header.size();

    for (auto i = 0; i < MyPID; ++i) {
        write_offset += write_size_per_proc[i];
    }


    MPI_File_delete (file_name.c_str(), MPI_INFO_NULL);
    MPI_File_open(comm, file_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);


    if (MyPID == 0) {
        MPI_File_write_at(fh, 0, header.c_str(), header.size(), MPI_CHAR, &sts);
    }

    MPI_File_write_at(fh, write_offset, data.c_str(), data.size(), MPI_CHAR, &sts);

    MPI_File_close(&fh);
    free(write_size_per_proc);
}

template<typename T, typename S>
void parVector<T, S>::writeToTxtCmplx(std::string file_name){
    std::string header;
    std::string data;

    if( (sizeof(T)/sizeof(Base<T>) != 2) ){
        try{
            throw 505;
        }catch(...){
            std::cout << "SMG2S]> Caught Exception: for real vector, please use writeToTxt" << std::endl;
        }               
    }

    //generate header
    header.append("%%SMG2S vector in complex scalar\n");
    std::string dim_info = std::to_string(global_size) + " " + std::to_string(global_size)  + " " + std::to_string(global_size) + "\n";
    header.append(dim_info);    

    for(auto i = 0; i < local_size; i++){
        data += std::to_string(index_map.Loc2Glob(i)+1) + " " + std::to_string(array[i].real()) + + " " + std::to_string(array[i].imag()) + "\n";
    }

    MPI_File fh;
    MPI_Offset write_offset;
    MPI_Offset text_size;
    MPI_Offset *write_size_per_proc;
    MPI_Status sts;

    write_size_per_proc = (MPI_Offset *)malloc(sizeof(MPI_Offset) * nProcs);
    
    text_size = data.size();

    MPI_Allgather(&text_size, 1, MPI_OFFSET, write_size_per_proc, 1, MPI_OFFSET, comm);

    write_offset = header.size();

    for (auto i = 0; i < MyPID; ++i) {
        write_offset += write_size_per_proc[i];
    }


    MPI_File_delete (file_name.c_str(), MPI_INFO_NULL);
    MPI_File_open(comm, file_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);


    if (MyPID == 0) {
        MPI_File_write_at(fh, 0, header.c_str(), header.size(), MPI_CHAR, &sts);
    }

    MPI_File_write_at(fh, write_offset, data.c_str(), data.size(), MPI_CHAR, &sts);

    MPI_File_close(&fh);
    free(write_size_per_proc);
}


template<typename T, typename S>
T parVector<T, S>::GetUpperNeighbor(S offset){
    MPI_Request	rtypereq, stypereq;
    MPI_Status	typestat;
    int up, down;

    if(MyPID != 0){
    	up = MyPID - 1;
    }else{
    	up = nProcs - 1;
    }

    if(MyPID != nProcs - 1){
    	down = MyPID + 1;
    }else{
    	down = 0;
    }

    T neighbor;

    try{
        if(offset < 1 || offset > local_size){
    	    throw 505;
        }
    }catch(...){
    	std::cout << "SMG2S]> Caught Exception: input parameter is invalid for GetUpperNeighbor" << std::endl;
    }

    MPI_Isend(&array[local_size-offset], 1, getMPI_Type<T>(), down, 1, comm, &stypereq);
    MPI_Irecv(&neighbor, 1, getMPI_Type<T>(), up, 1, comm, &rtypereq);

    MPI_Wait(&rtypereq,&typestat);

    return neighbor;
      
}

template<typename T, typename S>
T parVector<T, S>::GetLowerNeighbor(S offset){
    MPI_Request	rtypereq, stypereq;
    MPI_Status	typestat;
    int up, down;

    if(MyPID != 0){
    	up = MyPID - 1;
    }else{
    	up = nProcs - 1;
    }

    if(MyPID != nProcs - 1){
    	down = MyPID + 1;
    }else{
    	down = 0;
    }

    T neighbor;
    
    try{
        if(offset < 1 || offset > local_size){
    	    throw 505;
        }
    }catch(...){
    	std::cout << "SMG2S]> Caught Exception: input parameter is invalid for GetLowerNeighbor" << std::endl;
    }

    MPI_Isend(&array[offset-1], 1, getMPI_Type<T>(), up, 1, comm, &stypereq);
    MPI_Irecv(&neighbor, 1, getMPI_Type<T>(), down, 1, comm, &rtypereq);

    MPI_Wait(&rtypereq,&typestat);

    return neighbor;

}
#endif