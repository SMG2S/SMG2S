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

#ifndef __INIT_MAT_H__
#define __INIT_MAT_H__

/** @defgroup group4 initMat
 *  This module relates to the information of initialisation of matrix to be used by SMG2S.
 *  @{
 */
//!  @brief A struct which stores the information for the initial input matrix of SMG2S.
/*!
  This struct defines a initial input matrix by the offset of lower diagonal `diag_l` and upper
  diagonal `diag_u`. The entries between diagonal `diag_l` and `diag_u` of lower-triangular part
  of initial matrix are set randomly with (0, 1). These 
  entries can be modified by `scale`. The `sparsity` determines the possibility of the entries
  of initial matrix to be `0`.

  @tparam S type of integer to describes the dimension of matrices to be generated. 
*/
template<typename S>
struct initMat
{
    /*!
      Offset of lower diagonal
    */	
	S diag_l;
  	/*!
      Offset of upper diagonal
    */		
	S diag_u;
 	/*!
      Number can be multiplied on the randomly generated entries of initial matrix
    */		
	double scale;
  	/*!
      Possibility of the entries of initial matrix to be `0`
    */	
	double sparsity;

	//! A constructor of struct `initMat`.
    /*!
      This is a constructor of struct `initMat`, in which the variables are set
      with default values: `diag_l=-3`, `diag_u=-2`, `scale=1.0` and `sparsity=0.9`.
    */	
	initMat(){
		diag_l = -3;
		diag_u = -2;
		scale = 1.0;
		sparsity = 0.9;
	};

	//! A constructor of struct `initMat`.
    /*!
      This is a constructor of struct `initMat` with default values `scale=1.0` and `sparsity=0.9`.

      * @param[in] diagl Offset of lower diagonal
      * @param[in] diagu Offset of upper diagonal
    */	
	initMat(S diagl, S diagu){
		diag_l = diagl;
		diag_u = diagu;
		scale = 1.0;
		sparsity = 0.9;			
	};

	//! A constructor of struct `initMat`.
    /*!
      This is a constructor of struct `initMat` with default value `scale=1.0`.

      * @param[in] diagl Offset of lower diagonal
      * @param[in] diagu Offset of upper diagonal
      * @param[in] Sparsity Possibility of the entries of initial matrix to be `0`
    */	
	initMat(S diagl, S diagu, double Sparsity){
		diag_l = diagl;
		diag_u = diagu;
		scale = 1.0;
		sparsity = Sparsity;			
	};

	//! A constructor of struct `initMat`.
    /*!
      This is a constructor of struct `initMat` without any default values.

      * @param[in] diagl Offset of lower diagonal
      * @param[in] diagu Offset of upper diagonal
      * @param[in] Scale Number can be multiplied on the randomly generated entries of initial matrix
      * @param[in] Sparsity Possibility of the entries of initial matrix to be `0`
    */	
	initMat(S diagl, S diagu, double Scale, double Sparsity){
		diag_l = diagl;
		diag_u = diagu;
		scale = Scale;
		sparsity = Sparsity;			
	};		

	//! A member function of struct `initMat` to display all its variables.
    /*!
      This is a member function of struct `initMat` to display all its variables.
    */	
	void show(){
		std::cout << "Init Mat parameters:" << std::endl;
		std::cout << "        diag_l: " << diag_l << ", diag_u: " << diag_u << ", scale: " << scale << ", sparsity: " << sparsity << std::endl;
	};

};

/** @} */ // end of group4

#endif