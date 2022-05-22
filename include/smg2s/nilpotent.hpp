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

//T scalar, S integer
template<typename S>
class Nilpotent
{
	private:
		S probSize;
		S degree;
		S offset;
		std::vector<S> indOfZeros; // a pair: index

	public:
		Nilpotent();
		Nilpotent(S nbOne, S size);
		Nilpotent(S nbOne, S diag, S size);
		Nilpotent(std::vector<S> nilpvec, S size);
		Nilpotent(std::vector<S> nilpvec, S diag, S size);
		~Nilpotent();

		S computeDegree(std::vector<S> nilpvec);
		S getProbSize(){return probSize;};
		S getDegree(){return degree;};
		S getOffset(){return offset;};
		std::vector<S> getIndOfZeros(){return indOfZeros;};
		void show();

};

template<typename S>
Nilpotent<S>::Nilpotent(){}

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

	for(auto i = 0; i < nilpvec.size();i++){
		std::cout << nilpvec[i] << ", ";
	}
	std::cout << std::endl;

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

#endif
