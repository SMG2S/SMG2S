#include <vector>

template<typename T, typename S>
struct MatrixCSR
{
	S	nrows;
	S	nnz;

	std::vector<S> rows;
	std::vector<S> cols;
	std::vector<T> vals;

	MatrixCSR()
	{
		nnz = 0;
		rows = NULL;
		cols = NULL;
		vals = NULL;
	};

	MatrixCSR(S nnz_in, S nrows_in)
	{
	
		nnz = nnz_in;
		rows.reserve(nrows_in+1);
		cols.reserve(nnz);
		vals.reserve(nnz);	
	
	};

	~MatrixCSR()
	{
		if(nnz != 0){
			rows.empty();
			cols.empty();
			vals.empty();
		}
	};

	void Free()
	{
		if(nnz != 0){
			std::vector<S>().swap(rows);
			std::vector<S>().swap(cols);
			std::vector<T>().swap(vals);
			
			nnz = 0;
			rows = NULL;
			cols = NULL;
			vals = NULL;
		}
	};

};


