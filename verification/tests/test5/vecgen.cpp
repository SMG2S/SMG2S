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

	if (myfile.is_open()){
		myfile << "%%MatrixMarket matrix coordinate real general\n";
		myfile << N << " " << N << " "<< N << "\n";

		for (i=0;i<N/4;i++){
			myfile << i+1 << " "<< 10 + 5*sin(4*180*i/N) << " "<<  4*cos(4*180*i/N) << "\n";
		}
	
        	for (i=N/4;i<N/2;i++){
                        myfile << i+1 << " "<< 10 + 5*sin(4*180*(i-N/4)/N) << " "<<  -4*cos(4*180*(i-N/4)/N) << "\n";
                }	

                for (i=N/2;i<3*N/4;i++){
                        myfile << i+1 << " "<< 10 + 5*sin(4*180*(i-N/2)/N)+0.1 << " "<<  4*cos(4*180*(i-N/2)/N)+0.02 << "\n";
                }

                for (i=3*N/4;i<N;i++){
                        myfile << i+1 << " "<< 10 + 5*sin(4*180*(i-3*N/4)/N)+0.1 << " "<<  -4*cos(4*180*(i-3*N/4)/N)-0.02 << "\n";
                }

  	cout << "]> Finsh to generate the eigenvalues\n";
	}
	else cout << "Unable to open file";

  	return 0;
}
