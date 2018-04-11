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

	Nilpotency<int> a, b, c;

	printf("%d\n", a.diagPosition);

	a = Nilpotency<int>(2,2,10,3);
	
	printf("%d\n", a.diagPosition);

	b.NilpType1(1,20);

	printf("%d\n", b.matrix_size);

	c.NilpType2(2,20);

	if(c.setup == true)
	{
		printf("%d\n", c.nilpotency);
	}	
	
	return 0;
}