#include "MatrixCSR.cc"

using namespace std;

int main(int argc, char ** argv){

	MatrixCSR<int,int> matrix(3);

	matrix.set(-5,2,3);

	std::cout << matrix << std::endl;

	int val;

	val = matrix.get(2, 3);
	std::cout << val << std::endl;

	vector<int> vec(3,3), result;
	result = matrix * vec;
	
	int i;
	for(i=0; i<vec.size();i++){
		std::cout << vec[i] <<" ";
	}
	std::cout << std::endl;

        for(i=0; i<result.size();i++){
                std::cout << result[i] <<" ";
        }
        std::cout << std::endl;


	return 0;
}
