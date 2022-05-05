#include "../nilpotent/nilpotent.h"
#include <math.h>
#include <complex>
#include <cstdlib>
#include <string.h>
#include <string>
#include <typeinfo>  

int main(int argc, char** argv) 
{

	Nilpotent<int> nilp = Nilpotent<int>(2, 8);

	auto iz = nilp.getIndOfZeros();
/*
	for(auto i = 0; i < iz.size(); i++){
		std::cout << iz[i] << std::endl;
	}
*/
	auto nilp2 = Nilpotent<int>(10, 5, 40);
	auto iz2 = nilp2.getIndOfZeros();

	nilp2.show();
/*
	for(auto i = 0; i < iz2.size(); i++){
		std::cout << iz2[i] << std::endl;
	}
*/
	return 0;
}