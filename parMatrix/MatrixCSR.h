template<typename T, typename S>
struct MatrixCSR
{
	S	nrows;
	S	nnz;
	S	*rows;
	S	*cols;
	T	*vals;

	MatrixCSR()
	{
		nnz = 0;
		rows = NULL;
		cols = NULL;
		vals = NULL;
	};

	MatrixCSR(S nnz_in, S nrows_in)
	{
		if (nnz_in !=0 && nrows_in != 0){
			nnz = nnz_in;
			rows = new S[nrows_in+1];
			cols = new S[nnz];
			vals = new T[nnz];
		}
		else{
			nnz = 0;
			rows = NULL;
			cols = NULL;
			vals = NULL;
		}	
	
	};

	~MatrixCSR()
	{
		if(nnz != 0){
			delete [] rows;
			delete [] cols;
			delete [] vals;
		}
	};

	void Free()
	{
		if(nnz != 0){
			delete [] rows;
                        delete [] cols;
                        delete [] vals;
			
			nnz = 0;
			rows = NULL;
			cols = NULL;
			vals = NULL;
		}
	};

};


