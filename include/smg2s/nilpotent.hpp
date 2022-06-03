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

#ifndef __NILPOTENT_H__
#define __NILPOTENT_H__


#include "utils/utils.hpp"

/** @defgroup group3 Nilpotent
 *  This module relates to the nilpotent matrix used by SMG2S
 *  @{
 */
//!  @brief A class which determines the information of a nilpotent matrix.
/*!
  This class determines a special subset of <a href="https://en.wikipedia.org/wiki/Nilpotent">nilpotent matrix</a>,
  in which only one single off-diagonal on the upper-triangular part of a square matrix has non-zeros entries.
  These non-zeros entries are fixed to be `1`.

  This types of nilpotent matrix is determined by the four variables: Nilpotent<S>#probSize, Nilpotent<S>#degree, Nilpotent<S>#offset and Nilpotent<S>#indOfZeros.

  @tparam S type of integer to describes the dimension of matrices to be generated. 
*/
template<typename S>
class Nilpotent
{
	private:
		//! size of nilpotent matrix
    	/*!
      	size of nilpotent matrix (square matrix with number of rows = number of columns = probSize)
    	*/			
		S probSize;
		//! nilpotency of nilpotent matrix
  		/*!
      	Nipotency of the determined nilpotent matrix. For a nilpotent matrix `A`, it is the minimal integer `k` 
      		which makes \f$A^k = 0\f$.
    	*/		
		S degree;
		//! offset of off-diagonal
		/*!
		Offset of the the off-diagonal on the upper-triangular part of nilpotent matrix. 
		It is `0`-indexed, which means that the offset of main diagonal is `0`. 

		- Currently, the maximum degree supported is `80`
		*/
		S offset;
		//! indices of all zeros entries on the off-diagonal.
		/*!
		Indices of all zero entries on the only off-diagonal of a nilpotent matrix
		*/
		std::vector<S> indOfZeros; // a pair: index

	public:
		//! A constructor of class `Nilpotent`.
	    /*!
	      This is a constructor of class `Nilpotent`. With this constructor, the nilpotent has single off-diagonal
	      with index `1`. On this diagonal, the entries starts with `nbOne` number of continuous `1`, then one single
	      entry to be `0`, then `nbOne` number of continuous `1`. The pattern repeats until to the end.

	      In this pattern, the nilpotency is easy to be determined as \f$nbOne + 1\f$.
		  
		  - Currently, the maximum degree supported is `80`, there's risk that SMG2S fails with \f$nbOne > 80\f$

    	  * @param[in] nbOne number of continuous `1` on the single off-diagonal
      	  * @param[in] size size of nilpotent matrix
	    */		
		Nilpotent(S nbOne, S size);
		//! A constructor of class `Nilpotent`.
	    /*!
	      This is a constructor of class `Nilpotent`. With this constructor, the nilpotent has single off-diagonal
	      with index `diag`. On this diagonal, the entries starts with a number `k` of continuous `1`, in which `k`
	      is randomly generated as an integer between `0` and `nbOne`. Then a new `k` is re-generated,
	      for a continuous `k` of zero-entries. 

	      In this pattern, the nilpotency is not trivial and it can be computed by Nilpotent<S>#computeDegree.

		  - Currently, the maximum degree supported is `80`, there's risk that SMG2S fails with \f$nbOne - diagonal > 80\f$
	
    	  * @param[in] nbOne determines the maximum lengths of continous `1` and `0` on the off-diagonal. The lengths are randomly generated between `0` and `nbOne`, step by step.
    	  * @param[in] diag offset of the single off-diagonal of nilpotent matrix
      	  * @param[in] size size of nilpotent matrix
	    */		
		Nilpotent(S nbOne, S diag, S size);
		//! A constructor of class `Nilpotent`.
	    /*!
	      This is a constructor of class `Nilpotent`. It determines a nilpotent with user provided vector `nilpvec`, and
	      set it onto the off-diagonal indexing `1`.

	      In this pattern, the nilpotency is not trivial and it can be computed by Nilpotent<S>#computeDegree.

		  - Currently, the maximum degree supported is `80`, there's risk that SMG2S fails with user-given `nilpvec`.

    	  * @param[in] nilpvec user provided vector
      	  * @param[in] size size of nilpotent matrix
	    */			
		Nilpotent(std::vector<S> nilpvec, S size);
		//! A constructor of class `Nilpotent`.
	    /*!
	      This is a constructor of class `Nilpotent`. It determines a nilpotent with user provided vector `nilpvec`, and
	      set it onto the off-diagonal indexing `diag`.

	      In this pattern, the nilpotency is not trivial and it can be computed by Nilpotent<S>#computeDegree.

		  - Currently, the maximum degree supported is `80`, there's risk that SMG2S fails with user-given `nilpvec`.

    	  * @param[in] nilpvec user provided vector
    	  * @param[in] diag offset of the single off-diagonal of nilpotent matrix
      	  * @param[in] size size of nilpotent matrix
	    */			
		Nilpotent(std::vector<S> nilpvec, S diag, S size);
		//! A destructor of class `Nilpotent`. 
		~Nilpotent();
		//! compute the degree of a nilpotent with its off-diagonal entries stored in a vector `nilpvec`.
		/*!
	      This is a member function of class `Nilpotent`, which computes the degree of a nilpotent with its off-diagonal entries stored in a vector `nilpvec`.
		  
		  - Attention, the size of nilpotent matrix and the offset of diagonal should already be provided.

    	  * @param[in] nilpvec user provided vector
	    */		
		S computeDegree(std::vector<S> nilpvec);
		//! Get the size of nilpotent matrix
		/*!
			return Nilpotent<S>#probSize
		*/
		S getProbSize(){return probSize;};
		/*!
			return Nilpotent<S>#degree
		*/		
		//! Get the nilpotency of nilpotent matrix
		S getDegree(){return degree;};
		/*!
			return Nilpotent<S>#offset
		*/		
		//! Get off-diagonal offset of nilpotent matrix
		/*!
			return Nilpotent<S>#offset
		*/		
		S getOffset(){return offset;};
		//! Get indices of all zeros entries on the off-diagonal
		/*!
			return Nilpotent<S>#indOfZeros
		*/
		std::vector<S> getIndOfZeros(){return indOfZeros;};
		//! Display the information of a nilpotent matrix
		/*!
		    print Nilpotent<S>#probSize, Nilpotent<S>#degree, Nilpotent<S>#offset, and Nilpotent<S>#indOfZeros
		*/		
		void show();

};



template<typename S>
Nilpotent<S>::Nilpotent(S nbOne, S size){

    probSize = size;
    offset = 1; //in the default mode, the Nilpotent is generated with the diagonal of offset 1
    degree = nbOne + 1;

    if(nbOne > probSize - 1){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for Nilpotent, nbOne is too large" << std::endl;
    	}
    }

    for(auto i = 0; i < probSize - 1; i++){
    	if( (i + 1) % (nbOne + 1) == 0 ){
    		indOfZeros.push_back(i);
    	}
    }    	
}

template<typename S>
Nilpotent<S>::Nilpotent(S nbOne, S diag, S size){

    probSize = size;
    offset = diag; //in the default mode, the Nilpotent is generated with the diagonal of offset 1
    //degree = nbOne + 1;
    
    if(nbOne > probSize - diag){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for Nilpotent, nbOne is too large" << std::endl;
    	}
    }

    S length = probSize - diag;

    S cnt = 0;

    while (cnt < length){
    	S l1 =  1 + (rand() % nbOne);
    	S l2 =  1 + (rand() % nbOne);

    	for(auto i = l1 + cnt; i < l1 + l2 + cnt; i++){
    		if(i < length){
    			indOfZeros.push_back(i);
    		}
    	}

    	cnt += l1 + l2;
    }

	std::vector<S> nilpvec(length);
	std::fill(nilpvec.begin(), nilpvec.end(), 1);
	for(auto i = 0; i < indOfZeros.size(); i++){
		nilpvec[indOfZeros[i]] = 0;
	}

    degree = computeDegree(nilpvec);

}

template<typename S>
Nilpotent<S>::Nilpotent(std::vector<S> nilpvec, S size){
    probSize = size;
    offset = 1; 
    
    S nbOne = 0;
    S cnt = 0;

    if(nilpvec.size() < size - 1){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for Nilpotent, nilpvec size is not enough" << std::endl;
    	}    	
    }
    for(auto i = 0; i < size - 1; i++){
    	if (nilpvec[i] != 0){
    		cnt ++;
    	}else if(nilpvec[i] == 0){
    		indOfZeros.push_back(i);
    		if (cnt > nbOne){
    			nbOne = cnt;
    		}
    		cnt = 0;
    	}
    }

    degree = nbOne + 1;    
}
		
template<typename S>
Nilpotent<S>::Nilpotent(std::vector<S> nilpvec, S diag, S size){
    probSize = size;
    offset = diag; 

    if(nilpvec.size() < size - offset){
    	try{
    	    throw 505;
    	}catch(...){
    	    std::cout << "SMG2S]> Caught Exception: for Nilpotent, nilpvec size is not enough" << std::endl;
    	}  
    }
    for(auto i = 0; i < size - diag; i++){
    	if(nilpvec[i] == 0){
    		indOfZeros.push_back(i);
    	}
    }

    std::vector<S> vec(nilpvec.begin(), nilpvec.begin() + size - offset);
    degree = computeDegree(vec);   
}

template<typename S>
S Nilpotent<S>::computeDegree(std::vector<S> nilpvec){

	std::vector<S> vectmps(nilpvec.begin(), nilpvec.end());

	S diag_offset = offset;
	degree = 1;
	S maxdegree = 80;

	while (degree <= maxdegree){
		for(auto i = 0; i < probSize - offset - diag_offset; i++){
			vectmps[i] = nilpvec[i] * vectmps[i+offset];
		}

		degree = degree + 1;

		bool zeros = std::all_of(vectmps.begin(), vectmps.begin() + probSize-offset - diag_offset, [](S i) { return i==0; });
		if(zeros){
			break;
		}
		diag_offset = diag_offset + offset;
	}

	return degree;
}

template<typename S>
void Nilpotent<S>::show(){
	std::cout << "Nilpotent: " << "Matrix Size = " << probSize << ", diag_offset = " << offset << "; degree = " << degree << std::endl;
}

template<typename S>
Nilpotent<S>::~Nilpotent(){

}

/** @} */ // end of group3

#endif
