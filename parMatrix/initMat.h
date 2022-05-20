
/*

MIT License

Copyright (c) 2019 Xinzhe WU

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

#ifndef __INIT_MAT_H__
#define __INIT_MAT_H__

template<typename S>
struct initMat
{
	S diag_l;
	S diag_u;
	double scale;
	double sparsity;

	
	initMat(){
		diag_l = -3;
		diag_u = -2;
		scale = 1.0;
		sparsity = 0.9;
	};

	initMat(S diagl, S diagu){
		diag_l = diagl;
		diag_u = diagu;
		scale = 1.0;
		sparsity = 0.9;			
	};

	initMat(S diagl, S diagu, double Sparsity){
		diag_l = diagl;
		diag_u = diagu;
		scale = 1.0;
		sparsity = Sparsity;			
	};

	initMat(S diagl, S diagu, double Scale, double Sparsity){
		diag_l = diagl;
		diag_u = diagu;
		scale = Scale;
		sparsity = Sparsity;			
	};		


	void show(){
		std::cout << "Init Mat parameters:" << std::endl;
		std::cout << "        diag_l: " << diag_l << ", diag_u: " << diag_u << ", scale: " << scale << ", sparsity: " << sparsity << std::endl;
	};

};


#endif
