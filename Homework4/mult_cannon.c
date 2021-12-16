#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

/*
 * 定义节点间通信使用的 TAG
 * SHIFT_TAG 用于合并结果时的移动
 * SCATTER_A_TAG 和 SCATTER_B_TAG 分别在拆分任务和进行计算时使用
 * SCATTER_A_TAG 始终用于区分或处理矩阵 A
 * SCATTER_B_TAG 始终用于区分或处理矩阵 B
 */
#define SHIFT_TAG 77
#define SCATTER_A_TAG 88
#define SCATTER_B_TAG 99

/*
 * 将多个函数公用的变量提取为全局变量
 * 这种处理在减少调用传参时压栈开销的同时
 * 有利于编译器将变量优化为寄存器变量，从而提升性能
 * 三个矩阵的实际尺寸依次是 n x m, m x p, n x p
 */
int my_id, num_proc, n, m, p;
int rooted_proc, max_rows_a, max_cols_a, max_rows_b, max_cols_b;

/*
 * 交换两个 double* 类型数据
 * 一个更一般性的实现是采用 void* 类型的参数，使用 memcpy 完成交换
 * 考虑到本程序中只需要交换两个 double* 类型的数据，采用以下更简单高效的实现
 */
void swap(double** x, double** y) {
    double* temp = *x;
    *x = *y, *y = temp;
}

/*
 * 计算 (x, y) 向上循环移动 shift 个单位后的一维偏移
 * 采用作差对 rooted_proc 取模加 rooted_proc 再取模一次的方式
 * 可以保证取模结果处于区间 [0, rooted_proc) 而不会出现负数
 */
int shift_up(const int x, const int y, const int shift) {
    return ((x - shift) % rooted_proc + rooted_proc) % rooted_proc * rooted_proc + y;
}

/*
 * 计算 (x, y) 向左循环移动 shift 个单位后的一维偏移
 * 采用作差对 rooted_proc 取模加 rooted_proc 再取模一次的方式
 * 可以保证取模结果处于区间 [0, rooted_proc) 而不会出现负数
 */
int shift_left(const int x, const int y, const int shift) {
    return x * rooted_proc + ((y - shift) % rooted_proc + rooted_proc) % rooted_proc;
}

/*
 * 分发矩阵，通过 tag 对矩阵进行区分
 * tag 为 SCATTER_A_TAG 时处理 A 矩阵
 * tag 为 SCATTER_B_TAG 时处理 B 矩阵
 * row 和 col 分别对应矩阵的行数和列数
 */
void scatter_matrix(double* buf, const int row, const int col, const int tag, const double* buffer) {
    /*
     * 计算最大行数和最大列数
     * 缓存行数列数的乘积以提高程序性能
     */
    int max_rows = (row + rooted_proc - 1) / rooted_proc;
    int max_cols = (col + rooted_proc - 1) / rooted_proc;
    int total = max_rows * max_cols;
    if (!my_id) {
        /*
         * 只申请 buffer 一次，反复使用后释放
         * 避免多次动态申请内存造成程序性能下降
         */
        double* temp = (double*)malloc(sizeof(double) * total);
        for (int i = 0; i < rooted_proc; i++)  for (int j = 0; j < rooted_proc; j++) {
            /*
             * 遍历 i 和 j 并处理 (i, j) 节点对应的数据块
             * 对于不能整除时出现的空间冗余，将多余的位置赋 0 处理
             */
            for (int k = 0; k < max_rows; k++) for (int l = 0; l < max_cols; l++) {
                /*
                 * (i, j) 节点中位置 (k, l) 的元素对应了原矩阵中
                 * 位置 (i * max_rows_a + k, j * max_cols_b + l) 的元素
                 */
                temp[k * max_cols + l] = (i * max_rows + k < row && j * max_cols + l < col) ?
                        buffer[(i * max_rows + k) * col + j * max_cols + l] : 0;
            }
            /*
             * i 和 j 不全是 0 表明目标节点不是 0 号节点，此时需要使用 MPI 完成传送
             * 否则只需要由 0 号节点自身完成拷贝即可
             */
            if (i || j) {
                /*
                 * 对 A 矩阵的块 (i, j) 初始化时需要循环左移 i 列
                 * 对 B 矩阵的块 (i, j) 初始化是需要循环上移 j 列
                 */
                MPI_Send(temp, total, MPI_DOUBLE, tag == SCATTER_A_TAG ?
                        shift_left(i, j, i) : shift_up(i, j, j),
                        tag, MPI_COMM_WORLD);
            }
            else for (int k = 0; k < total; k++) {
                buf[k] = temp[k];
            }
        }
        free(temp);
    }
    else {
        /*
         * 非 0 节点不参与分发矩阵到其它节点
         * 只需要接收来自 0 号节点的矩阵即可
         */
        MPI_Status status;
        MPI_Recv(buf, total, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }
}

/*
 * 佳能算法的核心计算部分
 * 理论上可以使用 inline 关键字标记以使编译器内联此函数，提高性能
 * 然而在 mpicc 中开启 -O2 优化后 inline 关键字会诡异的导致 Runtime Error
 * 为使程序有更好的可读性，仍封装成函数
 */
void cannon(double* block_a, double* block_b, double* block_c
        , double* buffer_a, double* buffer_b, double* buffer_c) {
    /*
     * 初始化存储当前节点计算结果的块为 0
     * 在后续过程中不断对 block_c 进行累加
     */
    for (int i = 0; i < max_rows_a * max_cols_b; i++) {
        block_c[i] = 0;
    }
    for (int i = 0; i < rooted_proc; i++) {
        /*
         * 此处使用了类似 Copy On Write 的思想，维护了 block 和 buffer：
         *     其中的 block 是只读的，即传输过程和计算过程中都不会对其进行更新
         *     其中的 buffer 是只写的，即传输过程和计算过程中都不会对其进行读取
         * 上述性质允许了同时进行传输和计算，效率得以进一步提高
         */
        MPI_Status status[4];
        MPI_Request request[4];
        /*
         * 接下来四个函数调用对应了两个处理：
         *     将 A 矩阵中正在处理的块 (my_id / rooted_proc, my_id % rooted_proc) 循环左移一个单位
         *     将 B 矩阵中正在处理的块 (my_id / rooted_proc, my_id % rooted_proc) 循环上移一个单位
         * 开启传输后，在当前节点上串行完成矩阵乘法计算
         */
        MPI_Irecv(buffer_a, max_rows_a * max_cols_a, MPI_DOUBLE,
                shift_left(my_id / rooted_proc, my_id % rooted_proc, -1),
                SCATTER_A_TAG, MPI_COMM_WORLD, &request[0]);
        MPI_Irecv(buffer_b, max_rows_b * max_cols_b, MPI_DOUBLE,
                shift_up(my_id / rooted_proc, my_id % rooted_proc, -1),
                SCATTER_B_TAG, MPI_COMM_WORLD, &request[1]);
        MPI_Isend(block_a, max_rows_a * max_cols_a, MPI_DOUBLE,
                shift_left(my_id / rooted_proc, my_id % rooted_proc, 1),
                SCATTER_A_TAG, MPI_COMM_WORLD, &request[2]);
        MPI_Isend(block_b, max_rows_b * max_cols_b, MPI_DOUBLE,
                shift_up(my_id / rooted_proc, my_id % rooted_proc, 1),
                SCATTER_B_TAG, MPI_COMM_WORLD, &request[3]);
        for (int j = 0; j < max_rows_a; j++) for (int k = 0; k < max_cols_b; k++) {
            /*
             * 采用 temp 缓存，一方面保持 buffer 的只写性
             * 另一方面减少访存开销以提升效率
             */
            double temp = 0;
            for (int l = 0; l < max_cols_a; l++) {
                temp += block_a[j * max_cols_a + l] * block_b[l * max_cols_b + k];
            }
            buffer_c[j * max_cols_b + k] = temp;
        }
        /*
         * 传输和计算结束，开始交换 block 和 buffer 的角色
         * 对于结果 C 的交换不依赖于 A 和 B，因此可在同步前进行以提升效率
         * 对 A 和 B 的交换需要等待传输和计算完成，因此需要进行一次同步
         * 对 A 和 B 的交换无需遍历进行赋值，只需要交换两个指针即可
         */
        for (int j = 0; j < max_rows_a * max_cols_b; j++) {
            block_c[j] += buffer_c[j];
        }
        MPI_Waitall(4, request, status);
        swap(&block_a, &buffer_a), swap(&block_b, &buffer_b);
    }
}

/*
 * 计算结束，整合结果得到目标矩阵
 * 理论上可以使用 inline 关键字标记以使编译器内联此函数，提高性能
 * 然而在 mpicc 中开启 -O2 优化后 inline 关键字会诡异的导致 Runtime Error
 * 为使程序有更好的可读性，仍封装成函数
 */
void gather_matrix(double* source_c, double* block_c, double* buffer_c) {
    if (!my_id) {
        for (int i = 0; i < rooted_proc; i++) for (int j = 0; j < rooted_proc; j++) {
            /*
             * 遍历 i 和 j 并将其计算结果汇总到 0 号节点
             * i 和 j 不全是 0 表明目标节点不是 0 号节点，此时需要使用 MPI 完成传送
             * 否则只需要由 0 号节点自身完成拷贝即可
             */
            if (i || j){
                MPI_Status status;
                MPI_Recv(buffer_c, max_rows_a * max_cols_b, MPI_DOUBLE,
                         i * rooted_proc + j, SHIFT_TAG, MPI_COMM_WORLD, &status);
            }
            else for (int k = 0; k < max_rows_a * max_cols_b; k++) {
                buffer_c[k] = block_c[k];
            }
            /*
             * 整合 (i, j) 节点对应的数据块到合适的位置
             * 对于不能整除时出现的空间冗余，忽略多余的位置
             */
            for (int k = 0; k < max_rows_a; k++) for (int l = 0; l < max_cols_b; l++) {
                if (i * max_rows_a + k < n && j * max_cols_b + l < p) {
                    /*
                     * (i, j) 节点中位置 (k, l) 的元素对应了原矩阵中
                     * 位置 (i * max_rows_a + k, j * max_cols_b + l) 的元素
                     */
                    source_c[(i * max_rows_a + k) * p + j * max_cols_b + l] =
                            buffer_c[k * max_cols_b + l];
                }
            }
        }
    }
    else {
        /*
         * 非 0 节点不参与处理来自其它矩阵的节点
         * 只需要将自身的处理结果发送到 0 号节点即可
         */
        MPI_Send(block_c, max_rows_a * max_cols_b, MPI_DOUBLE,
                 0, SHIFT_TAG, MPI_COMM_WORLD);
    }
}

/*
 * 读文件过程，返回读取到的矩阵
 * 读取的矩阵的存储空间在此函数中完成申请
 * 需要外部在合适的时候释放申请的这部分内存
 * 通过 row 指针和 col 指针的方式，间接返回矩阵的大小信息
 */
double* read_file(const char* filename, int* row, int* col) {
    /*
     * 开文件并特判
     * 文件不存在直接报错并退出
     */
    FILE* file_ptr;
    if (!(file_ptr = fopen(filename, "r"))) {
        printf("Can't open file %s\n", filename);
        exit(-1);
    }
    /*
     * 获取文件大小并将文件一次性读入内存
     * 从文件前两个字获取矩阵的尺寸信息
     */
    struct stat f_stat;
    stat(filename, &f_stat);
    int f_size = f_stat.st_size;
    char* f_stream = (char*)malloc(f_size);
    fread(f_stream, sizeof(char), f_size, file_ptr);
    *row = ((int*)f_stream)[0], *col = ((int*)f_stream)[1];
    /*
     * 处理异常情况：
     * 矩阵的长宽必须是正数
     */
    if (*row <= 0 || *col <= 0) {
        printf("Matrix size error, %dx%d\n", *row, *col);
        exit(-1);
    }
    /*
     * 处理异常情况
     * 矩阵与其尺寸信息的总和不应超过文件大小
     */
    if (f_size < (sizeof(int) * 2 + sizeof(double) * *row * *col)) {
        printf("Actual size mismatches with stated size\n");
        exit(-1);
    }
    /*
     * 为读取出来的矩阵申请空间
     * 循环对读取结果进行赋值，返回结果
     */
    double* matrix_addr = (double*)(f_stream + sizeof(int) * 2);
    double* matrix = (double*)malloc(sizeof(double) * *row * *col);
    printf("       ---- %s: %d * %d Matrix -----\n", filename, *row, *col);
    for (int i = 0; i < *row; i++) for (int j = 0; j < *col; j++) {
        matrix[i * *col + j] = matrix_addr[i * *col + j];
    }
    free(f_stream), fclose(file_ptr);
    return matrix;
}

/*
 * 写文件过程，输出传入的矩阵
 * row 和 col 分别对应矩阵的行数和列数
 */
void write_file(const char* filename, const int row, const int col, const double* matrix) {
    /*
     * 申请存储空间，并记录矩阵信息
     * 存储到文件的矩阵，最开始的两个字分别表示其行数和列数
     * 接下来行数乘列数个字按行存储整个矩阵
     */
    int buf_size = (int)(sizeof(int) * 2 + sizeof(double) * row * col);
    double* buffer = (double*)malloc(buf_size);
    ((int*)buffer)[0] = row, ((int*)buffer)[1] = col;
    double* ptr = (double*)((int*)buffer + 2);
    printf("Result matrix: %d * %d\n", row, col);
    for (int i = 0; i < row; i++, putchar('\n')) for (int j = 0; j < col; j++) {
        /*
         * 循环遍历传入的矩阵
         * 输出矩阵信息并写入缓存
         */
        printf("%.4f  ", matrix[i * col + j]);
        ptr[i * col + j] = matrix[i * col + j];
    }
    FILE* file_ptr;
    /*
     * 开文件并特判
     * 文件不存在直接报错并退出，而不是只报错不退出
     * 否则会因空指针出现 Runtime Error
     */
    if (!(file_ptr = fopen(filename, "w"))) {
        printf("Can't open file %s\n", filename);
        exit(-1);
    }
    fwrite(buffer, sizeof(char), buf_size, file_ptr);
    free(buffer), fclose(file_ptr);
}

/*
 * 将运行时间(不考虑 IO 时间)输出到文件
 * 由于写结果的文件和标准输出均已经存在内容，另开一个文件写时间
 */
int write_time(const char* filename, const double time) {
    FILE* file_ptr;
    /*
     * 开文件并特判
     * 文件不存在直接报错并退出，而不是只报错不退出
     * 否则会因空指针出现 Runtime Error
     */
    if (!(file_ptr = fopen(filename, "w"))) {
        printf("Can't open file %s\n", filename);
        exit(-1);
    }
    fprintf(file_ptr, "Total time is %lf\n", time);
    fclose(file_ptr);
    return 0;
}

int main(int argc, char* argv[]) {
    /*
     * 初始化 MPI 环境并获取节点总数，自身 id
     * 在 cannon 算法中，节点总数必须为完全平方数
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    rooted_proc = (int)sqrt(num_proc);
    if (num_proc != rooted_proc * rooted_proc) {
        printf("Processor number must be a square!\n");
        exit(-1);
    }
    double* source_a = NULL; double* source_b = NULL; double* source_c = NULL;
    if (!my_id) {
        /*
         * 0 号节点读矩阵文件，并为计算结果申请空间
         * 读矩阵的同时完成 0 号节点中，对矩阵尺寸的更新
         */
        source_a = read_file("matrix_a.stdin", &n, &m);
        source_b = read_file("matrix_b.stdin", &m, &p);
        source_c = (double*)malloc(sizeof(double) * n * p);
    }
    /*
     * MPI_Bcast广播发送消息，将矩阵的大小信息广播到 comm 域
     * 由于每个节点都有收和发的操作，需要将 bcast 代码放在所有程序都能运行到的地方
     * 各个节点根据大小信息计算块的大小，并对各个块申请空间
     * 如果不能整除，则实际大小为向上取整后的大小
     */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_INT, 0, MPI_COMM_WORLD);
    max_rows_a = (n + rooted_proc - 1) / rooted_proc, max_cols_a = (m + rooted_proc - 1) / rooted_proc;
    max_rows_b = (m + rooted_proc - 1) / rooted_proc, max_cols_b = (p + rooted_proc - 1) / rooted_proc;
    double* block_a = (double*)malloc(sizeof(double) * max_rows_a * max_cols_a);
    double* buffer_a = (double*)malloc(sizeof(double) * max_rows_a * max_cols_a);
    double* block_b = (double*)malloc(sizeof(double) * max_rows_b * max_cols_b);
    double* buffer_b = (double*)malloc(sizeof(double) * max_rows_b * max_cols_b);
    double* block_c = (double*)malloc(sizeof(double) * max_rows_a * max_cols_b);
    double* buffer_c = (double*)malloc(sizeof(double) * max_rows_a * max_cols_b);
    /*
     * 空间申请完成，分别分发 A 矩阵和 B 矩阵到各个节点
     * 分发完成即表明 IO 完成，开始计时
     */
    scatter_matrix(block_a, n, m, SCATTER_A_TAG, source_a);
    scatter_matrix(block_b, m, p, SCATTER_B_TAG, source_b);
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    /*
     * 调用 cannon 算法核心核心计算部分完成计算
     * 完成计算后同步一次
     * 合并计算结果，得到目标矩阵
     * 完成计算后同步一次
     * 此时计算完成而 IO 即将开始，计时结束
     */
    cannon(block_a, block_b, block_c, buffer_a, buffer_b, buffer_c);
    MPI_Barrier(MPI_COMM_WORLD);
    gather_matrix(source_c, block_c, buffer_c);
    MPI_Barrier(MPI_COMM_WORLD);
    double finish_time = MPI_Wtime();
    /*
     * 开始输出，只有 0 号节点需要输出和释放三个 source 对应的空间
     * 全部节点都要释放其 block 和 buffer 对应的空间
     */
    if (my_id == 0) {
        write_time("run_time.stdout", finish_time - start_time);
        write_file("result.stdout", n, p, source_c);
        free(source_a), free(source_b), free(source_c);
    }
    free(buffer_a), free(buffer_b), free(buffer_c);
    free(block_a), free(block_b), free(block_c);
    MPI_Finalize(); return 0;
}