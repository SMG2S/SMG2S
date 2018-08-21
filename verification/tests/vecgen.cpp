// writing on a text file
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#define PI 3.14159265

using namespace std;

int main () {

  	cout << "]> Starting to generate the eigenvalues\n";

  	ofstream myfile ("vector.txt");
	int N = 100,i;

	srand((unsigned)time(0));

	double random_value, random_value2;

	if (myfile.is_open()){
		myfile << "%%MatrixMarket matrix coordinate real general\n";
		myfile << N << " " << N << " "<< N << "\n";

		for (i=0;i<N/5;i++){
			random_value = 0.1*i;
			random_value2 = -0.2*random_value + 0.6;
			myfile << i+1 << " "<< random_value << " "<<  random_value2 << "\n";
		}	

		for (i=N/5;i<2*N/5;i++){
			random_value = 0.1*(i - N/5);
			random_value2 = -0.2*random_value + 0.3;
			myfile << i+1 << " "<< random_value << " "<<  random_value2 << "\n";
		}

		for (i=2*N/5;i<3*N/5;i++){
			random_value = 0.1*(i - 2*N/5);
			random_value2 = 0.0;
			myfile << i+1 << " "<< random_value << " "<<  random_value2 << "\n";
		}

		for (i=3*N/5;i<4*N/5;i++){
			random_value = 0.1*(i - 3*N/5);
			random_value2 = 0.2*random_value - 0.3;
			myfile << i+1 << " "<< random_value << " "<<  random_value2 << "\n";
		}

		for (i=4*N/5;i<N;i++){
			random_value = 0.1*(i - 4*N/5);
			random_value2 = 0.2*random_value - 0.6;
			myfile << i+1 << " "<< random_value << " "<<  random_value2 << "\n";
		}

  	cout << "]> Finsh to generate the eigenvalues\n";
	}
	else cout << "Unable to open file";

  	return 0;
}