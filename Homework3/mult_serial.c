#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>

double w_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1E-6 * tv.tv_usec;
}

// Preprocess the command line, read in matrix A, B from input files, allocate
// memory for buffers, i.e, f_stream_a, f_stream_b to cache them. Suppose A, B's
// size are n1 * n2, n2 * n3, then n1 ~ n3 will be stored at dim[0 ~ 2]
// Return value 0 means no error occurred during pre processing, otherwise a
// non-zero returns.
int setup(int argc, char* argv[], char** f_stream_a, char** f_stream_b, int* dim) {
    if (argc < 4) {
        printf("Invalid arguments!\n");
        printf("Usage: ./serial filea fileb filec\n");
        printf("filea, fileb and filec are file names for matrix A, B and C\n");
        return 1;
    }

    FILE* file_ptr_a;
    FILE* file_ptr_b;
    if (!(file_ptr_a = fopen(argv[1], "r"))) {
        printf("Can't open matrix file %s, Errno=%d\n", argv[1], errno);
        return 1;
    }
    if (!(file_ptr_b = fopen(argv[2], "r"))) {
        printf("Can't open matrix file %s, Errno=%d\n", argv[2], errno);
        return 1;
    }
    struct stat f_stat_a, f_stat_b;
    stat(argv[1], &f_stat_a);
    stat(argv[2], &f_stat_b);
    int f_size_a = f_stat_a.st_size;
    int f_size_b = f_stat_b.st_size;
    *f_stream_a = (char*)malloc(f_size_a);
    *f_stream_b = (char*)malloc(f_size_b);
    fread(*f_stream_a, sizeof(char), f_size_a, file_ptr_a);
    fread(*f_stream_b, sizeof(char), f_size_b, file_ptr_b);

    int n1 = ((int*)*f_stream_a)[0], m1 = ((int*)*f_stream_a)[1];
    int n2 = ((int*)*f_stream_b)[0], m2 = ((int*)*f_stream_b)[1];
    if (n1 <=0 || m1 <= 0 || n2 <= 0 || m2 <= 0 || m1 != n2) {
        printf("Matrix size error, %dx%d with %dx%d\n", n1, m1, n2, m2);
        return 1;
    }
    if (f_size_a < sizeof(int) * 2 + sizeof(double) * n1 * m1) {
        printf("Actual size of A mismatches with stated size\n");
        return 1;
    }
    if (f_size_b < sizeof(int) * 2 + sizeof(double) * n2 * m2) {
        printf("Actual size of B mismatches with stated size\n");
        return 1;
    }
    dim[0] = n1, dim[1] = m1, dim[2] = m2;
    fclose(file_ptr_a), fclose(file_ptr_b);
    return 0;
}

// Compute C = A * B. A is a n1 * n2 matrix. B is a n2 * n3 matrix.
void matmul(const double* A, const double* B, double* C, const int n1, const int n2, const int n3) {
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            C[i * n3 + j] = 0.0;
            for (int k = 0; k < n2; k++) {
                C[i * n3 + j] += A[i * n2 + k] * B[k * n3 + j];
            }
        }
    }
}

int main(int argc, char* argv[]) {
    // Buffers to cache matrix files of A, B
    char* f_stream_a;
    char* f_stream_b;
    // Preprocess the command line, read in files for
    // A, B and put their sizes in dim[].
    int dim[3];
    if (setup(argc, argv, &f_stream_a, &f_stream_b, dim)) {
        // Something error during pre processing
        exit(-1);
    }

    // Suppose A's size is n1 x n2, B's is n2 x n3.
    // n1 ~ n3 are read from input files.
    int n1 = dim[0];
    int n2 = dim[1];
    int n3 = dim[2];
    FILE* file_ptr;
    int f_size_c = (int)(sizeof(int) * 2 + sizeof(double) * n1 * n3);
    if (!(file_ptr = fopen(argv[3], "w"))) {
        printf("Can't open file %s, Errno=%d\n", argv[3], errno);
        exit(-1);
    }

    char* f_stream_c = (char *)malloc(f_size_c);
    ((int*)f_stream_c)[0] = n1;
    ((int*)f_stream_c)[1] = n3;
    double elapsed_time = w_time();
    matmul((double*)(f_stream_a + sizeof(int) * 2),
        (double*)(f_stream_b + sizeof(int) * 2),
        (double*)(f_stream_c + sizeof(int) * 2),
        n1, n2, n3);
    elapsed_time = w_time() - elapsed_time;

    printf("Serial algorithm: multiply a %d x %d with a %d x %d, use %.2f(s)\n",
        n1, n2, n2, n3, elapsed_time);
    fwrite(f_stream_c, sizeof(char), f_size_c, file_ptr);
    free(f_stream_a);
    free(f_stream_b);
    free(f_stream_c);
    fclose(file_ptr);
    return 0;
}