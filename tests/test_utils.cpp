#include "../utils/utils.h"

#include <iostream>

int main(int argc, char **argv){

	std::cout << std::endl << "DEBUG >>> Test of Array Containing the information of Nilpotent matrix !!!" << std::endl; 

	std::map<int, int> nilp;

	nilp = nilpotent(50, 5);

	std::map<int, int>::iterator itr;

	for(itr = nilp.begin(); itr != nilp.end(); ++itr){
		printf("nilp[%d] = %d\n", itr->first, itr->second);
	}

	Nilpotency<int> a;

	printf("%d\n", a.diagPosition);

	a = Nilpotency<int>(2,2,10);
	
	printf("%d\n", a.diagPosition);
	
	return 0;
}