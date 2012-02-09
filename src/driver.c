#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "par_mult.h"
#include "seq_mult.h"

//double timeval_diff(struct timeval *t1, struct timeval *t2);
int ROOT = 0;

/* 
   To run this driver:
    mpiexec -n 4 ./run parallel-or-sequential ubfactor
      where parallel-or-sequential is 1 for the parallel version, and 0 for sequential, and
            ubfactor specifies the unbalance factor which is passed to hMETIS
*/

int main(int argc, char** argv)
{
   char* matrix_file = "/tmp/matrix.txt";
   char* vector_file = "/tmp/vector.txt";
   char* par_results_file = "par_results.txt";
   char* seq_results_file = "seq_results.txt";

   FILE* matrix_fp = fopen(matrix_file, "r");
   FILE* vector_fp = fopen(vector_file, "r");
 
   int parallel; 
   struct timeval t0;
   struct timeval t1;   
   int rank = ROOT;

   int ubfactor = 25;
   if (argc == 2 || argc == 3)
   { 
      parallel = atoi(argv[1]);
   }
   else
   {
      printf("Usage: run flag ubfactor\n");
      printf("    run - 0 for sequential version, and 1 for parallel\n");
      printf("    ubfactor - option unbalance factor for partitioning, default is 25\n");
      exit(-1);
   }
   if (argc==3)
   {  
      ubfactor = atoi(argv[2]);
   }

   gettimeofday(&t0, NULL);

   if (parallel == 1)
   {
      FILE* par_results_fp = fopen(par_results_file, "w");
      rank = matrix_multiply(argc, argv, ubfactor, matrix_fp, vector_fp, par_results_fp);
   }
   else
   {
      FILE* seq_results_fp = fopen(seq_results_file, "w");
      matrix_multiply_sequential(matrix_fp, vector_fp, seq_results_fp);
   }
   gettimeofday(&t1, NULL);

   if (rank == ROOT)
   {
      double timediff = timeval_diff(&t1, &t0);
      printf("Time - total: %2.6f\n", timediff);
   }
   if (parallel == 1)
   {
      finalize();
   }
}


