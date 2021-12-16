#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

int main(int argc, char* argv[]) {
    FILE* file_ptr;
    if (argc < 2) {
        printf("Invalid arguments!\n");
        printf("Run the program as ./print filename\n");
        exit(-1);
    }

    if (!(file_ptr = fopen(argv[1], "r"))) {
        printf("Can't open file %s\n", argv[1]);
        exit(-1);
    }

    struct stat f_stat;
    stat(argv[1], &f_stat);
    int f_size = f_stat.st_size;
    char* f_stream = (char *)malloc(f_size);
    fread(f_stream, sizeof(char), f_size, file_ptr);
    int n = ((int*)f_stream)[0], m = ((int*)f_stream)[1];
    double* matrix = (double*)(f_stream + sizeof(int) * 2);

    if (n <= 0 || m <= 0) {
        printf("Matrix size error, %dx%d\n", n, m);
        exit(-1);
    }
    if (f_size < (sizeof(int) * 2 + sizeof(double) * n * m)) {
        printf("Actual size mismatches with stated size\n");
        exit(-1);
    }

    printf("             ---- %s: %d * %d Matrix -----\n", argv[1], n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%.4f    ", matrix[i * m + j]);
        }
        putchar('\n');
    }
    free(f_stream), fclose(file_ptr);
    return 0;
}