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