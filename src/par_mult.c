#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include <string.h>
#include "par_mult.h"

#define ROOT_PROCESS   0
#define WORLD          MPI_COMM_WORLD
#define TAG            1

int EIND_SIZE         = 2000000;
int INITIAL_SEND_SIZE = 10000000;
MPI_Comm MY_COMM;

double timeval_diff(struct timeval *t1, struct timeval *t0);

/* NOTES 
  size == nvtxs == #rows

  TODO -Double eind array if it gets full (currently exit)
       -Double send_buff array if it gets full (currently exit)
       -Change matrix* to matrix** so getting a row isn't so expensive         
       -Assign the partition with the most vertices to root processor, if partition is not very load balanced
         this will save a lot of communication time
*/

int matrix_multiply(int argc, char** argv, int ubfactor, FILE* matrix_fp, FILE* vector_fp, 
                    FILE* results_fp)
{
   int NUM_PROCESSORS;
   int rank;

   MPI_Init(&argc, &argv);
   MPI_Comm_size(WORLD, &NUM_PROCESSORS);
   MPI_Comm_rank(WORLD, &rank);
 
   MPI_Comm_split(WORLD, 0, 0, &MY_COMM);
// printf("After split, my rank in MY_COMM is %i\n", rank);
 
   run(ubfactor, NUM_PROCESSORS, rank, matrix_fp, vector_fp, results_fp);
 
   return rank;
}

int finalize()
{
   MPI_Finalize();
   return 1;
}

int** build_hypergraph_partitioning(int ubfactor, int nparts, int rank, int size, 
                                    int* matrix, int* vector, int* part, int*** row_array) 
{
   int eind[EIND_SIZE];
   int eptr[size+1];
   int nhedges;          //create_hypergraph fills this in
   int eind_length;      //create_hypergraph fills this in

#ifndef NO_TIMING
   struct timeval t0;
   struct timeval t1;
   gettimeofday(&t0, NULL);
#endif

   create_hypergraph(matrix, size, ubfactor, nparts, eind, eptr, &nhedges, &eind_length);

#ifndef NO_TIMING
   gettimeofday(&t1, NULL);
   double timediff = timeval_diff(&t1, &t0);
   printf("Time - hypergraph creation: %2.6f\n", timediff);
   gettimeofday(&t0, NULL);
#endif

   int edgecut = get_partitioning(eind, eptr, size, nhedges, eind_length, ubfactor, nparts, part);

#ifndef NO_TIMING   
   gettimeofday(&t1, NULL);
   timediff = timeval_diff(&t1, &t0);
   printf("Time - HMetis call: %2.6f\n", timediff);
#endif

   int** send_array = build_send_array(nparts, size, part, matrix, vector, edgecut, 
                                       eptr, eind, row_array);
   //print_send_array(nparts, send_array);

   return send_array;
}

int* multiply(int* recv_buff, int rank, int* send_buff_size)
{
   int i; int j;
   int size; int ind; 
   int num_rows; int num_cols;
   int row_num; int row; int col_ind;
   int num_cuts; int recv_ind;
   int* send_buff;
   
   size = recv_buff[2];
   int* x = (int*)malloc(sizeof(int)*size);
   check(x, "x in multiply");
   int* y = (int*)malloc(sizeof(int)*size);
   check(y, "y in multiply");

// memset(x, 0, size*sizeof(int));  // TODO remove me
// memset(y, 0, size*sizeof(int));  // TODO remove me
   num_rows = recv_buff[3];
   ind = 4;
   for (i=0; i<num_rows;i++)  // figure out which rows i'm responsible for
   {
      row_num = recv_buff[ind+1];
      x[row_num] = recv_buff[ind+2];
      ind += recv_buff[ind];
   }
   num_cuts = (recv_buff[0] - recv_buff[1])/2; 
   ind = recv_buff[1]+1;           
   for (i=0; i<num_cuts; i++)
   {
      x[recv_buff[ind]] = recv_buff[ind+1];
      ind += 2;
   }
  
   ind = 4;
   for (i=0; i<num_rows; i++) 
   {
      num_cols = (recv_buff[ind]-3)/2;
      row = recv_buff[ind+1];
      y[row] = 0;
      col_ind = ind+3;
      for (j=0; j<num_cols; j++)
      {
         y[row] += x[recv_buff[col_ind]] * recv_buff[col_ind+1]; // multiply!
         col_ind += 2;
      }
      ind += recv_buff[ind];
   }

   if (rank == ROOT_PROCESS)
   {
      *send_buff_size = size;
      printf("Done multiplying %i rows.\n", num_rows);
      return y;
   }
   else
   {
      *send_buff_size = num_rows*2+1;
      send_buff = (int*)malloc(sizeof(int)*(*send_buff_size));
      send_buff[0] = num_rows;

      ind = 1;
      recv_ind = 4;
      if (rank != ROOT_PROCESS)
      {
         for (i=0; i<num_rows; i++)
         {
            row = recv_buff[recv_ind+1];
            send_buff[ind++] = row;
            send_buff[ind++] = y[row];
            recv_ind += recv_buff[recv_ind];
         }     
      }
      printf("Done multiplying %i rows.\n", num_rows);
      return send_buff;
   }
}

void gather_results(int* y, int y_size, int rank, int nparts, int** row_array)
{
   int* results_buff;
   MPI_Status status; 
   int i; int j;
   int num_results;
   int recv_size;
   int y_ind;
   int res_ind;

   if (rank != ROOT_PROCESS) // send y-results back to root
   {
      MPI_Send(y, y_size, MPI_INT, ROOT_PROCESS, TAG, MY_COMM);
   }
   else
   {
      for (i=1; i<nparts; i++)
      {
         recv_size = (row_array[i][0]-1)*2+1;
         results_buff = (int*)malloc(sizeof(int)*recv_size);
   
//       printf("receving %i bytes from rank %i\n", recv_size, i);
         MPI_Recv(results_buff, recv_size, MPI_INT, i, TAG, MY_COMM, &status);
 
         num_results = results_buff[0];
         res_ind = 1;
         for (j=0; j<num_results; j++)
         {
            y_ind = results_buff[res_ind++];
            y[y_ind] = results_buff[res_ind++];
         }
         free(results_buff);
      }
   }
}


void run(int ubfactor, int nparts, int rank, FILE* matrix_fp, FILE* vector_fp, FILE* results_fp)
{ 
   int i;
   int size;
   int* matrix;
   int* vector;
   int* part;
   int** send_array;
   MPI_Status status; 
   int recv_buff_size;
   int* recv_buff;
   int** row_array;
  
#ifndef NO_TIMING
   struct timeval t0;
   struct timeval t1; 
#endif

   // parent process reads matrix, creates hypergraph, gets partitioning, and 
   //  builds arrays to distribute
   if (rank == ROOT_PROCESS)
   {
      matrix = read_matrix(&size, matrix_fp);
     //print_matrix(matrix, size);
      
      vector = read_vector(size, vector_fp);
      //print_vector(vector, size);

#ifndef NO_TIMING
      gettimeofday(&t0, NULL);
#endif 

      part = (int*)malloc(sizeof(int)*size);
      check(part, "part in run");
      send_array = build_hypergraph_partitioning(ubfactor, nparts, rank, size, matrix, 
                                                 vector, part, &row_array);
   }
   MPI_Barrier(MY_COMM); // so recvs dont sit there and wait while build_hypergraph is called

   if (rank==ROOT_PROCESS)
   {
      // for (i=0; i<nparts; i++) 
      //   printf("row_array[%i][0] = %i\n", i, row_array[i][0]);
     
      for (i=1; i<nparts; i++)
      {  
         MPI_Send(&send_array[i][0], 1, MPI_INT, i, TAG, MY_COMM); // send size
      }
      for (i=1; i<nparts; i++) 
      {
         MPI_Send(send_array[i], send_array[i][0], MPI_INT, i, TAG, MY_COMM); // send data of size
      }
   }
   else    
   {
      MPI_Recv(&recv_buff_size, 1, MPI_INT, ROOT_PROCESS, TAG, MY_COMM, &status);
      recv_buff = (int*)malloc(sizeof(int)*recv_buff_size);
      check(recv_buff, "recv_buff in run");
      
      MPI_Recv(recv_buff, recv_buff_size, MPI_INT, ROOT_PROCESS, TAG, MY_COMM, &status);
      //print_array(recv_buff, recv_buff_size);
   }

   if (rank == ROOT_PROCESS)
   {
      recv_buff = send_array[0]; 
   }
   int results_size;
   int *results = multiply(recv_buff, rank, &results_size);

   gather_results(results, results_size, rank, nparts, row_array);

    if (rank==ROOT_PROCESS)
    {

#ifndef NO_TIMING
      gettimeofday(&t1, NULL);
#endif
      write_vector(results, size, results_fp);

#ifndef NO_TIMING
      double timediff = timeval_diff(&t1, &t0);
      printf("Time - inner: %2.6f\n", timediff);
#endif      
      free(matrix); 
      free(vector);
   }
}

int get_partitioning(int* eind, int* eptr, int nvtxs, int nhedges, int eind_count, int ubfactor, 
                         int nparts, int* part)
{
   int* vwgts  = NULL;
   int* hewgts = NULL;
   int options[9];
   options[0] = 0;
   int edgecut;

   printf("calling HMETIS_PartRecursive with following params:\n");
   printf(" nvtxs:    %i\n nhedges:  %i\n nparts:   %i\n ubfactor: %i\n",
                   nvtxs, nhedges, nparts, ubfactor);

   HMETIS_PartRecursive(nvtxs, nhedges, vwgts, eptr, eind, hewgts, nparts, ubfactor, options,
                          part, &edgecut);
   printf("HMETIS_PartRecursive returned. Edgecut: %i\n", edgecut);

   return edgecut;
}


void print_partition(int nvtxs, int* part)
{
   int i;
   for (i=0; i<nvtxs; i++)
   {
      printf(" %i |", i);
   }
   printf("\n");
   for (i=0; i<nvtxs; i++)
   {
      printf(" %i |", part[i]);
   }
   printf("\n");
}

int* read_vector(int insize, FILE* fp)
{
   //FILE* fp = fopen(vector_file, "r");
   // first line is vector size
   int size;
   fscanf(fp, "%d", &size);
   if (size != insize)
   {
      printf("Vector size %i does not match matrix size %i. Terminating.\n", size, insize);
      exit(-1);
   }
   int* vector = (int*)malloc(sizeof(int)*size);
   check(vector, "vector in read_vector");

   int i;
   for (i=0; i<size; i++)
   {
      fscanf(fp, "%d", &vector[i]);
   }
   //fclose(fp);

   return vector;
}

int* read_matrix(int* outSize, FILE* fp)
{
   int size;
   //FILE* fp = fopen(matrix_file, "r");
   // first line is matrix size
   fscanf(fp, "%d", &size);
   printf("size: %d\n", size);
   
   // now save into matrix
   int* matrix = (int *)malloc(size*size*sizeof(int));
   check(matrix, "matrix in read_matrix");
   printf("Matrix malloced\n");
 
   int row;
   int col;
   for (row=0; row<size; row++)
   {
      for (col=0; col<size; col++)
      {
         fscanf(fp, "%i", &(matrix[size*row+col]));
      }
   }
   printf("Finished reading.\n");

   //fclose(fp);

   *outSize = size;
   return matrix;
}

void create_hypergraph(int* matrix, int size, int ubfactor, int nparts, int* eind, int* eptr, 
                           int* eptr_ind_out, int* total_eind_out)
{
   int row; int col;
   int eind_count = 0;
   int total_eind = 0;
   int eptr_ind   = 1;
   int eptr_count = 0;
   eptr[0] = 0;

   for(col=0; col<size; col++)
   {
      eind_count = 0;

      for(row=0; row<size; row++)
      {

         if (total_eind >= EIND_SIZE)
         {
            printf("Ran out of space in eind array\n");
            exit(-1);
         }

         if (matrix[size*row+col] != 0)
         {
            eind[eind_count+total_eind] = row;
            eind_count++;
         }
      }
      eptr[eptr_ind] = eind_count + eptr_count;
      eptr_count += eind_count;
      eptr_ind++;
      total_eind += eind_count;

      if (total_eind >= EIND_SIZE)
      {
         printf("Ran out of space in eind array\n");
         exit(-1);
      }
   }
   
   eptr_ind--; // hmetis wants eptr to be of size nhedges+1

   // advance eptr array to point to end of eind
   eptr[eptr_ind] = total_eind;
   // note that eptr_ind == nhedges and total_eind == eind_size
   *eptr_ind_out = eptr_ind;
   *total_eind_out = total_eind;
}

void print_vector(int* v, int size)
{
   int i;
   printf("[ ");
   for (i=0; i<size; i++)
   {
      printf("%d ", v[i]);
   }
   printf("]\n");
}

void print_matrix(int* m, int size)
{
   int i; int j;
   for(i=0; i<size; i++)
   {
      for(j=0; j<size; j++)
      {
         printf("%d ", m[size*i+j]);
      }
      printf("\n");
   }
}


void print_hypergraph(int* eind, int* eptr, int eind_size, int eptr_size, int eind_count)
{
   printf("\neptr: ");
   int i;
   for (i=0; i<=eptr_size; i++)
   {
      printf("%d ", eptr[i]);
   }

   printf("\neind: ");

   for (i=0; i<eind_count; i++)
   {
      printf("%d ", eind[i]);
   }
 
   printf("\n\n");
}

// Returns true if elem does not exist in array, false otherwise
int new_element(int* array, int alloced_size, int size, int elem)
{
   int i;
   for (i=0; i<size; i++)
   {
      if (array[i] == elem) return 0;
   }
   return 1;
}

void build_row_array(int size, int nparts, int** row_array, int* row_array_indices, 
                     int* part)
{
   int i;
   int which_part;
   int index;
   for (i=0; i<nparts; i++)
   {
      // init to 1, save size in send_array[i][0]
      row_array_indices[i] = 1; 
      // TODO assuming there is some sort of load balance this array can be much smaller than size
      row_array[i] = (int*)malloc(sizeof(int)*size);   
                                                       
      check(row_array[i], "row_array[i] in build_row_array");
   } 
   for (i=0; i<size; i++)
   { 
      which_part = part[i];
      index = row_array_indices[which_part];
      row_array[which_part][index] = i;
      row_array_indices[which_part] = index+1;
   }
   for (i=0; i<nparts; i++)
   {
      // store the size of send_array[i] in send_array[i][0]
      row_array[i][0] = row_array_indices[i]; 
   }
}

int** build_send_array(int nparts, int size, int* part, int* m, int* v, int edgecut, 
                       int* eptr, int* eind, int*** row_array_in)
{      
   int i;
   int ind; int val;
   int row; int col;
   int p; int size_per_row;
   // stores vertices that need to be sent to each partition
   int** row_array = (int**)malloc(sizeof(int)*nparts);  
   check_dbl(row_array, "row_array in build_send_array");
   int* row_array_indices = (int*)malloc(sizeof(int)*nparts);
   check(row_array_indices, "row_array_indices in build_send_array");
 
   // STEP 1: build row_array - stores which rows need to go to which processors //
   build_row_array(size, nparts, row_array, row_array_indices, part);

   // STEP 2: build structure to send to each processor, use sparse matrix rep //
   // total_size|row1size|row1|col1|val1|col2|val2|...|row2size|row2|col1|val2|.... ///
   int initial_size = INITIAL_SEND_SIZE;   
   int** send_buff = (int**)malloc(sizeof(int)*nparts);
   check_dbl(send_buff, "send_buff in build_send_array");
   int* row_buff;
   
   for (p=0; p<nparts; p++)
   {
      ind = 4;
      row_buff = (int*)malloc(sizeof(int)*initial_size);
      check(row_buff, "row_buff in build_send_array");

      int num_rows = row_array[p][0];
      row_buff[2] = size;
      row_buff[3] = (num_rows == 0) ? num_rows : num_rows-1;
      //printf("nparts=%i num_rows=%i\n", nparts, num_rows);
      for (i=1; i<num_rows; i++)
      {
         row = row_array[p][i];
         size_per_row = 3;     
         // store number of entries at ind
         row_buff[ind+1] = row;       // store row 
         row_buff[ind+2] = v[row];    // store v[row]
         for (col=0; col<size; col++)
         {
            if (ind+size_per_row+100 >= INITIAL_SEND_SIZE)
            {
               printf("Ran out of space in row_buff.\n");
               //MPI_Finalize();
               exit(-1);
            }
            val = m[row*size+col];
            if (val != 0)
            {
               row_buff[ind+size_per_row] = col;
               size_per_row++;
               row_buff[ind+size_per_row] = val;
               size_per_row++;;
            }
         }            
         row_buff[ind] = size_per_row;
         ind+=size_per_row;
         //printf("size_per_row=%i at row=%i\n", size_per_row, row);
      }
      row_buff[0] = ind;
      row_buff[1] = (ind>0) ? (ind-1) : ind; 
      send_buff[p] = &row_buff[0];  // send_buff** points to row_buff! 
   } 


   // STEP 3: build cut_nets - stores the processor who contains vi of a cut net i //
   //   check if all vertices ended up on the same processor

   //struct timeval t0;
   //struct timeval t1;
   //gettimeofday(&t0, NULL);

   int cut_count = 0;
   int eptr_ind;
   int vertexc;
   int vertex; 
   int are_different;
   int num_vertices;
   int new_size;
   int processor_owns_xi;
   // array of processors that contain vertices corresponding to hyperedge eptr_ind
   int* processors; 

   for (eptr_ind=0; eptr_ind<size; eptr_ind++) // for each hyperedge
   {
      are_different = 0; // false initially
      num_vertices = eptr[eptr_ind+1] - eptr[eptr_ind];
      processors = (int*)malloc(sizeof(int)*num_vertices);
      check(processors, "processors in build_send_array"); 

      for (vertexc=0; vertexc<num_vertices; vertexc++)
      {
         vertex = eind[eptr[eptr_ind]+vertexc];
         processors[vertexc] = part[vertex];

         if (vertexc > 0 && (processors[vertexc] != processors[vertexc-1]))
         {
            are_different = 1;
         }
      }
      if (are_different == 1) // cut net!
      {
         //printf("***net-cut on hyperedge %i ", eptr_ind); 
         processor_owns_xi = part[eptr_ind];
         //printf("processor %i at part[%i] should send x-%i to: ", processor_owns_xi, 
         //         eptr_ind, eptr_ind);
         cut_count++;
         
         processors = remove_dupes(processors, num_vertices, &new_size, processor_owns_xi);
         //for (i=0; i<new_size; i++) printf("%i ", processors[i]); printf("\n");
         
         // enter data into send_buff
         for (i=0; i<new_size; i++) 
         {
            int curr_proc = processors[i];
            int ind = send_buff[curr_proc][0];
            send_buff[curr_proc][ind++] = eptr_ind;  // save i
            send_buff[curr_proc][ind++] = v[eptr_ind]; 
            send_buff[curr_proc][0]+=2;
         }
      }
      free(processors); 
      if (cut_count == edgecut) break;
   }
   //printf("Step 3 complete.\n");
   //free(row_array);
   free(row_array_indices);
   *row_array_in = row_array;

   //gettimeofday(&t1, NULL);
   //double timediff = timeval_diff(&t1, &t0);
   //printf("Net-cut discovery time: %2.6f\n", timediff);
   return send_buff;
}

int* remove_dupes(int* array, int size, int* new_size, int skip)
{
   int i; int j;
   int* array_out = (int*)malloc(sizeof(int)*size);

   int size_out = 0;
   int exists; 
   for (i=0; i<size; i++)
   {
      exists = 0;   // false initially
      if (array[i] != skip)
      {
         for (j=0; j<size_out; j++)
         {
            if (array[i] == array_out[j]) 
            {
               exists = 1;
               break; // if it exists then stop checking
            }
         }
         if (exists == 0)
         {
            array_out[size_out] = array[i];
            size_out++;
         }
      }
   }
   *new_size = size_out;
   // free the memory that was passed in since we are returning a new pointer
   free(array); 
   return array_out;
}

void print_send_array(int nparts, int** send_array)
{
   printf("\nprint_send_array():\n");
   int i;
   int j;
   int size;
   // print send_array
   for (i=0; i<nparts; i++)
   {
      size = send_array[i][0];  // size of send_array[i]
      printf("Parition #%i: ", i);
      for (j=0; j<size; j++)
      {
         printf("%i ", send_array[i][j]);
      }
      printf(" size: %i\n", size);
   }
}

void print_array(int* array, int size)
{
   int i;
   for(i=0; i<size; i++)
   {
      printf("%i ", array[i]);
   }
}

void check(int* ptr, char* str)
{
   if (ptr == NULL)
   {
      printf("Error mallocing: ");
      printf("%s\n", str);
      MPI_Finalize();
      exit(-1);
   }
}

void check_dbl(int** ptr, char* str)
{
   if (ptr == NULL)
   {
      printf("Error mallocing: ");
      printf("%s\n", str);
       MPI_Finalize();
      exit(-1);
   }
}

void write_vector(int* v, int size, FILE* fp)
{
   int i;
   //FILE* fp = fopen(filename, "w");

   fprintf(fp, "%d\n", size);

   for (i=0; i<size; i++)
   {
      fprintf(fp, "%i\n", v[i]);
   }
   //fclose(fp);
}

double timeval_diff(struct timeval *t1, struct timeval *t0)
{
   return ((t1->tv_sec - t0->tv_sec) +
            (t1->tv_usec - t0->tv_usec)/1.0e6);
}

