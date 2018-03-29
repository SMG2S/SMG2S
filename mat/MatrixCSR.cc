#include <vector>
#include <iostream>
#include "exceptions.h"
#include "MatrixCSR.h"

using namespace std;

//creation
template<typename T, typename S>
MatrixCSR<T,S>::MatrixCSR(S n)
{
	this->construct(n, n);
}

template<typename T, typename S>
MatrixCSR<T,S>::MatrixCSR(S rows, S columns)
{
	this->construct(rows, columns);
}

template<typename T, typename S>
MatrixCSR<T,S>::MatrixCSR(const MatrixCSR<T,S> & matrix)
{
	this->deepCopy(matrix);
}

template<typename T, typename S>
MatrixCSR<T,S> & MatrixCSR<T,S>::operator = (const MatrixCSR<T,S> & matrix)
{
	if (&matrix != this) {
		this->destruct();
		this->deepCopy(matrix);
	}

	return *this;
}

template<typename T, typename S>
void MatrixCSR<T,S>::deepCopy(const MatrixCSR<T,S> & matrix)
{
	this->m = matrix.m;
	this->n = matrix.n;
	this->rows = new vector<S>(*(matrix.rows));

	if (matrix.vals != NULL) {
		this->cols = new vector<S>(*(matrix.cols));
		this->vals = new vector<T>(*(matrix.vals));
	}
}

template<typename T, typename S>
MatrixCSR<T,S>::~MatrixCSR(void)
{
	this->destruct();
}

template<typename T, typename S>
void MatrixCSR<T,S>::construct(S rows, S columns)
{
	if (rows < 1 || columns < 1) {
		throw InvalidDimensionsException("Matrix dimensions cannot be zero or negative.");
	}

	this->m = rows;
	this->n = columns;

	this->vals = NULL;
	this->cols = NULL;
	this->rows = new vector<S>(rows + 1, 1);
}

template<typename T, typename S>
void MatrixCSR<T,S>::destruct(void)
{
	if (this->vals != NULL) {
		delete this->vals;
		delete this->cols;
	}

	delete this->rows;
}

//get set

template<typename T, typename S>
S MatrixCSR<T,S>::getRowCount(void) const
{
	return this->m;
}

template<typename T, typename S>
S MatrixCSR<T,S>::getColumnCount(void) const
{
	return this->n;
}

//values

template<typename T, typename S>
T MatrixCSR<T,S>::get(S row, S col) const
{
	this->validateCoordinates(row, col);

	S currCol;

	for (S pos = (*(this->rows))[row - 1] - 1; pos < (*(this->rows))[row] - 1; ++pos) {
		currCol = (*(this->cols))[pos];

		if (currCol == col) {
			return (*(this->vals))[pos];

		} else if (currCol > col) {
			break;
		}
	}

	return T();
}

template<typename T, typename S>
MatrixCSR<T,S> & MatrixCSR<T,S>::set(T val, S row, S col)
{
	this->validateCoordinates(row, col);

	S pos = (*(this->rows))[row - 1] - 1;
	S currCol = 0;

	for (; pos < (*(this->rows))[row] - 1; pos++) {
		currCol = (*(this->cols))[pos];

		if (currCol >= col) {
			break;
		}
	}

	if (currCol != col) {
		if (!(val == T())) {
			this->insert(pos, row, col, val);
		}

	} else if (val == T()) {
		this->remove(pos, row);

	} else {
		(*(this->vals))[pos] = val;
	}

	return *this;
}


//operations

template<typename T, typename S>
vector<T> MatrixCSR<T,S>::multiply(const vector<T> & x) const
{
	if (this->n != (S) x.size()) {
		throw InvalidDimensionsException("Cannot multiply: Matrix column count and vector size don't match.");
	}

	vector<T> result(this->m, T());

	if (this->vals != NULL) { // only if any value set
		for (S i = 0; i < this->m; i++) {
			T sum = T();
			for (S j = (*(this->rows))[i]; j < (*(this->rows))[i + 1]; j++) {
				sum = sum + (*(this->vals))[j - 1] * x[(*(this->cols))[j - 1] - 1];
			}

			result[i] = sum;
		}
	}

	return result;
}

template<typename T, typename S>
vector<T> MatrixCSR<T,S>::operator * (const vector<T> & x) const
{
	return this->multiply(x);
}

template<typename T, typename S>
MatrixCSR<T,S> MatrixCSR<T,S>::multiply(const MatrixCSR<T,S> & m) const
{
	if (this->n != m.m) {
		throw InvalidDimensionsException("Cannot multiply: Left matrix column count and right matrix row count don't match.");
	}

	MatrixCSR<T,S> result(this->m, m.n);

	T a;

	// TODO: more efficient?
	// @see http://www.math.tamu.edu/~srobertp/Courses/Math639_2014_Sp/CRSDescription/CRSStuff.pdf

	for (S i = 1; i <= this->m; i++) {
		for (S j = 1; j <= m.n; j++) {
			a = T();

			for (S k = 1; k <= this->n; k++) {
				a = a + this->get(i, k) * m.get(k, j);
			}

			result.set(a, i, j);
		}
	}

	return result;
}


template<typename T, typename S>
MatrixCSR<T,S> MatrixCSR<T,S>::operator * (const MatrixCSR<T,S> & m) const
{
	return this->multiply(m);
}

template<typename T, typename S>
MatrixCSR<T,S> MatrixCSR<T,S>::add(const MatrixCSR<T,S> & m) const
{
	if (this->m != m.m || this->n != m.n) {
		throw InvalidDimensionsException("Cannot add: matrices dimensions don't match.");
	}

	MatrixCSR<T,S> result(this->m, this->n);

	// TODO: more efficient?
	// @see http://www.math.tamu.edu/~srobertp/Courses/Math639_2014_Sp/CRSDescription/CRSStuff.pdf

	for (S i = 1; i <= this->m; i++) {
		for (S j = 1; j <= this->n; j++) {
			result.set(this->get(i, j) + m.get(i, j), i, j);
		}
	}

	return result;
}

template<typename T, typename S>
MatrixCSR<T,S> MatrixCSR<T,S>::operator + (const MatrixCSR<T,S> & m) const
{
	return this->add(m);
}

template<typename T, typename S>
MatrixCSR<T,S> MatrixCSR<T,S>::subtract(const MatrixCSR<T,S> & m) const
{
	if (this->m != m.m || this->n != m.n) {
		throw InvalidDimensionsException("Cannot subtract: matrices dimensions don't match.");
	}

	MatrixCSR<T,S> result(this->m, this->n);

	// TODO: more efficient?
	// @see http://www.math.tamu.edu/~srobertp/Courses/Math639_2014_Sp/CRSDescription/CRSStuff.pdf

	for (S i = 1; i <= this->m; i++) {
		for (S j = 1; j <= this->n; j++) {
			result.set(this->get(i, j) - m.get(i, j), i, j);
		}
	}

	return result;
}


template<typename T, typename S>
MatrixCSR<T,S> MatrixCSR<T,S>::operator - (const MatrixCSR<T,S> & m) const
{
	return this->subtract(m);
}


//helps and validators

template<typename T, typename S>
void MatrixCSR<T,S>::validateCoordinates(S row, S col) const
{
	if (row < 1 || col < 1 || row > this->m || col > this->n) {
		throw InvalidCoordinatesException("Coordinates out of range.");
	}
}


template<typename T, typename S>
void MatrixCSR<T,S>::insert(S index, S row, S col, T val)
{
	if (this->vals == NULL) {
		this->vals = new vector<T>(1, val);
		this->cols = new vector<S>(1, col);

	} else {
		this->vals->insert(this->vals->begin() + index, val);
		this->cols->insert(this->cols->begin() + index, col);
	}

	for (S i = row; i <= this->m; i++) {
		(*(this->rows))[i] += 1;
	}
}

template<typename T, typename S>
void MatrixCSR<T,S>::remove(S index, S row)
{
	this->vals->erase(this->vals->begin() + index);
	this->cols->erase(this->cols->begin() + index);

	for (S i = row; i <= this->m; i++) {
		(*(this->rows))[i] -= 1;
	}
}

template<typename T, typename S>
bool operator == (const MatrixCSR<T,S> & a, const MatrixCSR<T,S> & b)
{
	return ((a.vals == NULL && b.vals == NULL)
				|| (a.vals != NULL && b.vals != NULL && *(a.vals) == *(b.vals)))
			&& ((a.cols == NULL && b.cols == NULL)
				|| (a.cols != NULL && b.cols != NULL && *(a.cols) == *(b.cols)))
			&& *(a.rows) == *(b.rows);
}

template<typename T, typename S>
bool operator != (const MatrixCSR<T,S> & a, const MatrixCSR<T,S> & b)
{
	return !(a == b);
}

template<typename T, typename S>
ostream & operator << (ostream & os, const MatrixCSR<T,S> & matrix)
{
	for (S i = 1; i <= matrix.m; i++) {
		for (S j = 1; j <= matrix.n; j++) {
			if (j != 1) {
				os << " ";
			}

			os << matrix.get(i, j);
		}

		if (i < matrix.m) {
			os << endl;
		}
	}

	return os;
}


