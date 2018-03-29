#ifndef __MatrixCSR_H__
#define	__MatrixCSR_H__

#include <vector>
#include <iostream>

using namespace std;

template<typename T, typename S>
class MatrixCSR{
	public:

	//creation
	MatrixCSR(S n);
	MatrixCSR(S rows, S columns);

	MatrixCSR(const MatrixCSR<T,S> & m);
	MatrixCSR<T,S> & operator = (const MatrixCSR<T,S> & m);

	~MatrixCSR(void);

	//get / set

	S getRowCount(void) const;
	S getColumnCount(void) const;

	//values
	
	T get(S row, S col) const;
	MatrixCSR & set(T val, S row, S col);

	//operations
	vector<T> multiply(const vector<T> & x) const;
	vector<T> operator * (const vector<T> & x) const;

	MatrixCSR<T,S> multiply(const MatrixCSR<T,S> & m) const;
	MatrixCSR<T,S> operator * (const MatrixCSR<T,S> & m) const;

	MatrixCSR<T,S> add(const MatrixCSR<T,S> & m) const;
	MatrixCSR<T,S> operator + (const MatrixCSR<T,S> & m) const;

	MatrixCSR<T,S> subtract(const MatrixCSR<T,S> & m) const;
	MatrixCSR<T,S> operator - (const MatrixCSR<T,S> & m) const;

	//friend functions

	template<typename X, typename Y>
	friend bool operator == (const MatrixCSR<X, Y> & a, const MatrixCSR<X, Y> & b);

	template<typename X, typename Y>
	friend bool operator != (const MatrixCSR<X, Y> & a, const MatrixCSR<X, Y> & b);

	template<typename X, typename Y>
	friend ostream & operator << (ostream & os, const MatrixCSR<X, Y> & matrix);


	protected:

	S m, n;

	vector<T> * vals;
	vector<T> * rows, * cols;

	//help and validators

	void construct(S m, S n);
	void destruct(void);
	void deepCopy(const MatrixCSR<T,S> & m);
	void validateCoordinates(S row, S col) const;
	void insert(S index, S row, S col, T val);
	void remove(S index, S row);

};

#endif
