#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

// Preprocess the command line, read in matrix A, B from input files, allocate
// memory for buffers, i.e, f_stream_a, f_stream_b to cache them.
int setup(const int argc, char* argv[], char** f_stream_a, char** f_stream_b) {
    if (argc < 3) {
        printf("Invalid arguments!\n");
        printf("Usage: ./serial filea fileb\n");
        printf("filea, fileb are file names for matrix A, B to be compared\n");
        return 1;
    }

    FILE* file_ptr_a;
    FILE* file_ptr_b;
    if (!(file_ptr_a = fopen(argv[1], "r"))) {
        printf("Can't open matrix file %s, Errno = %d\n", argv[1], errno);
        return 1;
    }
    if (!(file_ptr_b = fopen(argv[2], "r"))) {
        printf("Can't open matrix file %s, Errno = %d\n", argv[2], errno);
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
    if (n1 <= 0 || m1 <= 0 || n2 <= 0 || m2 <= 0) {
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
    fclose(file_ptr_a);
    fclose(file_ptr_b);
    return 0;
}

void comp(const char* f_stream_a, const char* f_stream_b) {
    int n1 = ((int*)f_stream_a)[0], m1 = ((int*)f_stream_a)[1];
    int n2 = ((int*)f_stream_b)[0], m2 = ((int*)f_stream_b)[1];
    if (n1 != n2 || m1 != m2) {
        printf("Matrix size mismatch, %dx%d with %dx%d\n", n1, m1, n2, m2);
        return;
    }

    double norm = 0.0;
    double* A = (double*)(f_stream_a + sizeof(int) * 2);
    double* B = (double*)(f_stream_b + sizeof(int) * 2);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < m1; j++) {
            double diff = A[i * m1 + j] - B[i * m1 + j];
            norm += diff * diff;
        }
    }
    norm = sqrt(norm);
    if (norm > 0.000001) {
        printf("Matrix compare failed, with norm = %.8f\n", norm);
    }
    else {
        printf("Matrix compare succeeded, with norm = %.8f\n", norm);
    }
}

int main(int argc, char* argv[]) {
    // Buffers to cache matrix files of A, B
    char* f_stream_a;
    char* f_stream_b;
    // preprocess the command line, read in files for A, B.
    if (setup(argc, argv, &f_stream_a, &f_stream_b)) {
        // Something error during pre processing
        exit(-1);
    }
    comp(f_stream_a, f_stream_b);
    free(f_stream_a);
    free(f_stream_b);
    return 0;
}