#include <mpi.h> // Input/output
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <string>

#include <pVecMap.h>

template<typename T, typename S>
class pVec{
	private:
		T *array;
		S array_size;
		S local_size;

		pVecMap<T,S> *index_map;
 
	public:
	pVec();
	pVec(MPI_Comm ncomm, int lbound, int ubound);
	~pVec();

	pVecMap *GetVecMap(){return index_map;};
	
	S GetLowerBound();
	S GetUpperBound();
	S GetGlocalSize();
	S GetLocalSize();
	S GetArraySize();
	T *GetArray(){return array;};

	void AddValueLocal(S row, T value);
	void AddValuesLocal(S nindex, S *rows, T *values);
	void SetTovalue(T value);
	void SetToZero();

	void ReadExtVec();
	
	S Loc2Glob(S local_index);
	S Glob2Loc(S global_index);

	void RestoreArray(){};

};
