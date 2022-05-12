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


	nilp.show();

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

	std::vector<int> v{ 1,1,1,1,0,0,1,1,0,0,0,1};

	auto nilp3 = Nilpotent<int>(v,12);
	nilp3.show();

	auto nilp4 = Nilpotent<int>(v,3,15);
	nilp4.show();

	auto nilp5 = Nilpotent<int>(v,4,15);
	nilp5.show();

	return 0;


}