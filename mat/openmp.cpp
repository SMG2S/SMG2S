#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

#include <omp.h>

int main(int argc, char *argv[])
{


  //Below code will be run once for each processor (there are two)
    #pragma omp parallel 
  {
    cout << omp_get_thread_num() << endl; //this should output 1 and 0, in random order
  }


  //omp_set_nested(2);
  //The parallel example:
  vector <double> a(50000,0);

  clock_t start = omp_get_wtime();
#pragma omp parallel for  shared(a) 
  for (int i=0; i < 50000; i++)    
    {
      double StartVal=i;
#pragma simd
//#pragma vector aligned
//#pragma omp parallel for  shared(a, StartVal)
      for (int j=0; j<40000; ++j)
	a[i]=(StartVal + log(exp(exp((double) i)))); 
    } 

  cout<< "OpenMP Time: " << ( (double) ( omp_get_wtime() - start ) ) <<endl;

  //The serial example:
  start = omp_get_wtime();

  for (int i=0; i < 50000; i++)    
    {
      double StartVal=i;

      for (int j=0; j<40000; ++j)
	a[i]=(StartVal + log(exp(exp((double) i)))); 
    } 

  cout<< "Time: " << ( (double) ( omp_get_wtime() - start ) ) <<endl;

  return 0;
}
