#include "../parMatrix/parMatrixSparse.cc"
#include <math.h>
#include <complex.h>
#include <string.h>

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

int main(int argc, char** argv){

	int printEnds = 1, printMids = 1;

        int i,j,k;

        int span, lower_b, upper_b;

        const int tag1 = 0, tag2 = 1;

        int rank, size, intBuf1, intBuf2, number;

        double doubBuf1, start, finish, time;

        double a = 5;

        MPI_Status status;


        int commMethod = 0;
        int probSize = 10;

        int maxCount = 1;

	printf("HELLO\n");

        if (argc>=2){
                for (int i =0; i < argc; i++){
                        if (strcasecmp(argv[i],"-comm")){
                                commMethod = atoi(argv[i+1]) ;
                        }
                        if (strcasecmp(argv[i],"-size")==0){
                                probSize = atoi(argv[i+1]) ;
                        }
                        if (strcasecmp(argv[i],"-count")==0){
                                maxCount = atoi(argv[i+1]) ;
                        }

                }
        }

	return 0;
}

