#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "par_mult.h"

int* multiply_sequential(int* m, int* x, int size)
{
   int i;
   int j;
   int* y = (int*)malloc(sizeof(int)*size);
   if (y==NULL)
   { 
      printf("Error mallocing %i bytes\n", (sizeof(int)*size));
      exit(-1); 
   } 
   for (i=0; i<size; i++)
   {
      //printf("row %i\n", i);
      y[i] = 0;
      for (j=0; j<size; j++)
      {
        // printf("setting y[%i] = m[%i*%i+%i] * x[%i]", y[i], i, size, j, x[j]);
        // printf("m  %i  x %i", m[i*size+j], x[j]);
         y[i] += m[i*size+j] * x[j];
      }
   }  
   return y;
}

int matrix_multiply_sequential(FILE* matrix_file, FILE* vector_file, FILE* results_file)
{
   int size;
   struct timeval t0;
   struct timeval t1;

   int* m = read_matrix(&size, matrix_file);
   printf("Read matrix.\n");
   int* x = read_vector(size, vector_file);
   printf("Read vector.\n");

   gettimeofday(&t0, NULL);
   int* y = multiply_sequential(m,x,size);
   gettimeofday(&t1, NULL);

   double timediff = timeval_diff(&t1, &t0);
   printf("Time: %2.6f\n", timediff);

   printf("Done multiplying.\n");
   write_vector(y, size, results_file);

   return 0;
}

