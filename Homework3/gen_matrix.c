#include <time.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        printf("Invalid arguments!\n");
        printf("Run the program as ./gen_matrix n m filename\n");
        printf("n, m are row/col number of the matrix, filename is the file to write\n");
        return 0;
    }

    char* filename = argv[3];
    int n = (int)strtol(argv[1], NULL, 10);
    int m = (int)strtol(argv[2], NULL, 10);
    int buf_size = (int)(sizeof(int) * 2 + sizeof(double) * n * m);
    double* matrix = (double*)malloc(buf_size);

    ((int*)matrix)[0] = n;
    ((int*)matrix)[1] = m;
    double* ptr = (double*)((int*)matrix + 2);
    srand(n * m * 0x19260817);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            ptr[i * m + j] = (double)random() / RAND_MAX;
        }
    }

    FILE* file_ptr;
    if (!(file_ptr = fopen(filename, "w"))) {
        printf("Can't open file %s\n", filename);
        exit(-1);
    }
    fwrite(matrix, sizeof(char), buf_size, file_ptr);
    fclose(file_ptr);
    return 0;
}