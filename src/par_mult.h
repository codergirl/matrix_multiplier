
int matrix_multiply(int argc, char** argv, int ubfactor, FILE* matrix_fp, FILE* vector_fp, FILE* results_fp);
int finalize();

void run(int ubfactor, int nparts, int rank, FILE* matrix_fp, FILE* vector_fp, FILE* results_fp);
int** build_hypergraph_partitioning(int ubfactor, int nparts, int rank, int size, int* matrix, 
                                     int* vector, int* part, int*** row_array); 
int* read_matrix(int* outSize, FILE* matrix_file);
int* read_vector(int size, FILE* vector_file);
void create_hypergraph(int* matrix, int size, int ubfactor, int nparts, int* eind, int* eptr, 
                         int* nhedges, int* eind_length);
void print_send_array(int nparts, int** send_array);
int** build_send_array(int nparts, int size, int* part, int* matrix, int* vector, int edgecut, 
                         int* eptr, int* eind, int*** row_array);
int get_partitioning(int* eind, int* eptr, int nvtxs, int nhedges, int eind_count, int ubfactor, 
                         int nparts, int* part);
void print_partition(int nvtxs, int* part);
void print_hypergraph(int* eind, int* eptr, int eind_size, int eptr_size, int eind_count);
void print_matrix(int* m, int size);
void print_vector(int* v, int size);
void print_array(int* array, int size);
int* remove_dupes(int* array, int size, int* new_size, int skip);
int* multiply(int* buff, int rank, int* send_buff_size);
void check_dbl(int** ptr, char* str);
void check(int* ptr, char* str);
void write_vector(int* v, int size, FILE* fp);

double timeval_diff(struct timeval *t1, struct timeval *t0);


// external functions linked in with gcc
void HMETIS_PartKway(int nvtxs, int nhedges, int *vwgts, int *eptr, int *eind, int *hewgts, 
                      int nparts, int ubfactor, int *options, int *part, int *edgecut);
void HMETIS_PartRecursive(int nvtxs, int nhedges, int *vwgts, int *eptr, int *eind, int *hewgts, 
                      int nparts, int ubfactor, int *options, int *part, int *edgecut);


