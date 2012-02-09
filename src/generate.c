#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAX_VAL 1000
int negativity = 8;

/*
  To generate a matrix, run: 
       generate size sparsity
       where size is the number of rows/cols
             sparsity should be an integer greater than 0, it specifies how sparse the matrix should be
                      higher values of sparsity mean more sparse matrices
*/
void generate_matrix(int size, int sparsity, int diagonality, char* matrix_file) 
{
   srand(time(NULL));
   int i;
   int j;
   int nonzero;
   int negative;
   int val;
//   int diag;
//   int diff;
//   float spars = 1/(float)sparsity;

   FILE* fp = fopen(matrix_file, "w");   
   if (fp == NULL)
   {
      printf("Error opening file %s for writing.\n", matrix_file);
      exit(-1);
   }
   else
   {
      printf("Writing matrix to %s.\n", matrix_file);
   }
   fprintf(fp, "%d\n", size);

   for (i=0; i<size; i++)
   {
      nonzero = 1;
      for (j=0; j<size; j++)
      {
         //diff = abs(i-j);
         //if (diff == 0) diff = 1;
//         diag = abs(i-j)*diagonality;
//         printf("diff=%i, diag=%i\n", abs(i-j), diag);
//         if (diag == 0 || diag/sparsity == 0) diag = 1;
         //diag = (int)diff*spars;
         //if (diag == 0) diag = 1;
         //diff = 2*size - abs(i-j);
         //if (diff == 0) diff = 1;
         nonzero = (rand()%sparsity == 0);  


      //   printf("diff=%i, factor=%i\n", diff, diag);
         negative = (rand()%negativity == 0);

         if (nonzero == 1)
         {
            val = rand()%MAX_VAL;
         }
         else
         {
            val = 0;
         }

         if (i==j)  // force nonzero on diagonal
         {
            nonzero = 1;
            while (val == 0)
            {
               val = rand()%MAX_VAL;
            }
         }
         if (negative == 1)
         {
            val = val * (-1);
         }
         fprintf(fp, "%d ", val);
      }
      fprintf(fp, "\n");
   }

   fclose(fp);
}

void generate_vector(int size, char* vector_file)
{
   srand(time(NULL));
   int i;
   int val;

   FILE* fp = fopen(vector_file, "w");
   if (fp == NULL)
   {
      printf("Error opening file %s for writing.\n", vector_file);
      exit(-1);
   }
   else
   {
      printf("Writing vector to %s.\n", vector_file);
   }
   fprintf(fp, "%d\n", size);

   for (i=0; i<size; i++)
   {
      // dont want zeros
      val = 0;
      while (val == 0)
      {
         val = rand()%MAX_VAL;
      }      
      fprintf(fp, "%d\n", val);
   } 

   fclose(fp);
}

int main(int argc, char **argv)
{
   int size;
   int sparsity;
   int diagonality = 10;
 
   char* matrix_file = "/tmp/matrix.txt";
   char* vector_file = "/tmp/vector.txt";

   if (argc == 3 || argc == 4)
   { 
      size  = atoi(argv[1]); 
      sparsity = atoi(argv[2]);
    //  if (argc >= 5)
    //  {
    //     matrix_file = argv[3];
    //     vector_file = argv[4];
    //  }
      if (argc == 4)
      {
         diagonality = atoi(argv[3]);
      }
   } 
   else
   {
      printf("Usage: generate_matrix size sparsity diagonality filename\n  size - number of rows/columns\n");
      printf("  sparsity - integer > 0, higher sparsity indices result in sparser matrices\n");
      printf("  diagonality - 0 for uniform distribution, 1 for diagonal\n");
      printf("  filename - optional filename, 'matrix.txt' is default\n");
      exit(-1);
   }
   generate_matrix(size, sparsity, diagonality, matrix_file);
   generate_vector(size, vector_file);
   return 0;
}

