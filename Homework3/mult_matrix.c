#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>

/*
 * 定义矩阵的大小
 * 计算加速比时采用较大的居中，以便得到清晰稳定的区别
 */
#define LEN_N 2100
#define LEN_M 2100
#define LEN_P 2100

/*
 * 将计算最小值的过程封装
 * 由于比较简单，采用宏函数实现
 */
#define min(a, b) (((a) < (b)) ? (a) : (b))

/*
 * 用于计算的三个矩阵
 * 由于矩阵较大，从栈空间移到全局变量区
 */
double A[LEN_N][LEN_M], B[LEN_M][LEN_P], C[LEN_N][LEN_P];

/*
 * 定义线程的参数
 * 将传参封装成一个结构体
 */
struct threadArg {
    /*
     * tid：线程的 id
     * A_row：相乘行指针
     * C_row：结果行指针
     * num_threads：线程个数
     * len_m：A_row 的真实长度(而不是理论长度 LEN_M )
     * len_p：C_row 的真实长度(而不是理论长度 LEN_P )
     * 由于 B 已经移动到全局区，不需要通过参数传递
     */
    int tid;
    double* A_row;
    double* C_row;
    int num_threads;
    int len_m, len_p;
};

/*
 * 实际计算的过程
 * 完成结果行的计算
 */
void* worker(void* arg) {
    struct threadArg* my_arg = (struct threadArg*)arg;
    // printf("m = %d, p = %d\n", my_arg->len_m, my_arg->len_p);
    /*
     * 此处不需要开多少枚举多少，只需要用多少枚举多少，因此
     *     将理论长度 LEN_P 改为真实长度 len_p 枚举
     *     将理论长度 LEN_M 改为真实长度 len_m 枚举
     * 这种方式下(对于小矩阵)性能将大大提升
     *
     * 源码中的问题1：
     *     矩阵的尺寸在 LEN_N LEN_M LEN_P 处被写死，不利于提高计算效率
     */
    for (int i = my_arg->tid; i < my_arg->len_p; i += my_arg->num_threads) {
        /*
         * 动态申请的空间在堆中
         * 由于初值是不确定的，使用前需要先清空
         */
        my_arg->C_row[i] = 0.0;
        for (int j = 0; j < my_arg->len_m; j++) {
            /*
             * 串行的枚举完成计算即可
             * 将 B 中第 j 行第 i 个元素乘 A 对应的第 j 个元素
             * 累加到 C 对应的第 i 个元素
             */
            my_arg->C_row[i] += my_arg->A_row[j] * B[j][i];
        }
    }
    return NULL;
}

/*
 * 读文件过程，修改了 matrix 的内存使用
 * 具体的， matrix 每一行的存储长度是 width 以便使用 matrix[i] 访问
 * matrix 实际的大小则通过 n 和 m 作为结果输出
 */
int read_file(const char* filename, double* matrix, const int width, int* n, int* m) {
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
    *n = ((int*)f_stream)[0], *m = ((int*)f_stream)[1];
    /*
     * 处理异常情况：
     * 矩阵的长宽必须是正数
     */
    if (*n <= 0 || *m <= 0) {
        printf("Matrix size error, %dx%d\n", *n, *m);
        exit(-1);
    }
    /*
     * 处理异常情况
     * 矩阵与其尺寸信息的总和不应超过文件大小
     */
    if (f_size < (sizeof(int) * 2 + sizeof(double) * *n * *m)) {
        printf("Actual size mismatches with stated size\n");
        exit(-1);
    }

    /*
     * 此处原本是 matrix[i * *m + j] = matrix_addr[i * *m + j] 这样可以最小化矩阵占用的空间
     * 但矩阵已经以数组的形式开出，实际空间总在编译阶段确定，因此可以改用 matrix[i * width + j]
     * 使用 matrix[i * width + j] 的好处在于与 C 语言实际存储二维数组的方式一致
     * 因此可以方便的使用 matrix[i] 访问其第 i 行
     */
    double* matrix_addr = (double*)(f_stream + sizeof(int) * 2);
    printf("       ---- %s: %d * %d Matrix -----\n", filename, *n, *m);
    for (int i = 0; i < *n; i++) {
        for (int j = 0; j < *m; j++) {
            matrix[i * width + j] = matrix_addr[i * *m + j];
        }
    }
    free(f_stream), fclose(file_ptr);
    return 0;
}

/*
 * 写文件过程，由于修改了 matrix 的内存使用，需要同时修改输出方式
 * 具体的，存储到文件的矩阵每一行不存在空间冗余
 */
int write_file(const char* filename, const int n, const int m, const int width, const double* matrix) {
    /*
     * 分配存储空间，并记录矩阵信息
     * 存储到文件的矩阵，其最开始两个字分别表示其行数和列数
     * 接下来行数乘列数个字按行存储整个矩阵
     */
    int buf_size = (int)(sizeof(int) * 2 + sizeof(double) * n * m);
    double* buffer = (double*)malloc(buf_size);
    ((int*)buffer)[0] = n;
    ((int*)buffer)[1] = m;
    double* ptr = (double*)((int*)buffer + 2);
    printf("Result matrix: %d * %d\n", n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            /*
             * 输出矩阵信息并写入缓存
             * 读取时的处理为 matrix[i * *m + j] = matrix_addr[i * *m + j]
             * 写入时的处理为 ptr[i * m + j] = matrix[i * width + j]
             * 二者恰好相反
             */
            printf("%.4f  ", matrix[i * width + j]);
            ptr[i * m + j] = matrix[i * width + j];
        }
        putchar('\n');
    }

    FILE* file_ptr;
    /*
     * 开文件并特判
     * 文件不存在应直接报错并退出，而不是只报错不退出
     * 否则会因空指针出现 run time error
     */
    if (!(file_ptr = fopen(filename, "w"))) {
        printf("Can't open file %s\n", filename);
        exit(-1);
    }
    fwrite(buffer, sizeof(char), buf_size, file_ptr);
    fclose(file_ptr);
    return 0;
}

/*
 * 将运行时间(不考虑 IO 时间)输出到文件
 * 由于写结果的文件和标准输出均已经存在内容，另开一个文件写时间
 */
int write_time(const char* filename, const double time) {
    FILE* file_ptr;
    /*
     * 开文件并特判
     * 文件不存在应直接报错并退出，而不是只报错不退出
     * 否则会因空指针出现 run time error
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
     * 记录矩阵实际的尺寸
     * 三个矩阵的实际尺寸依次是 n x m, m x p, n x p
     */
    int n, m, p;
    /*
     * 初始化 MPI 环境
     * status 用于获取 MPI_TAG
     */
    MPI_Status status;
    int my_id, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    if (!my_id) {
        /*
         * 0 号节点读矩阵文件
         * 读矩阵的同时完成 0 号节点中，对矩阵尺寸的更新
         */
        read_file("matrix_a.stdin", (double*)A, LEN_M, &n, &m);
        read_file("matrix_b.stdin", (double*)B, LEN_P, &m, &p);
    }

    /*
     * MPI_Bcast广播发送消息，将 root 的 buffer 广播到 comm 域
     * 由于每个节点都有收和发的操作，需要将bcast代码放在所有程序都能运行到的地方
     */
    MPI_Bcast(B, LEN_M * LEN_P, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*
     * 除了矩阵 B 之外，广播矩阵的尺寸
     * 函数 void* worker(void* arg) 中将使用 m 和 p
     */
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&p, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (!my_id) {
        /*
         * IO 完成，在主节点(即 0 号节点)开始计时
         * 确定第一轮发送的数目，即取 MPI 节点数目和矩阵行数的较小者
         */
        double start_time = MPI_Wtime();
        int num_send = min(num_proc - 1, n);
        for (int i = 1; i < num_proc; i++) {
            /*
             * 源码中的问题2：
             *     循环终止条件有误，A[num_proc] 将在接下来的循环中计算
             *     因此 i < num_proc 不取等，否则将导致计算结果出错
             *
             * 向第 i 个节点发送矩阵的 i - 1 行
             * 需要注意，当节点多矩阵行数少的时候，需要向多余的节点发送结束消息
             * 否则这些节点将一直卡在 MPI_Recv 上，即死锁
             *
             * 源码中的问题3：
             *     没有特判节点数多矩阵行数少的情况，将导致死锁
             */
            if (i <= num_send) {
                MPI_Send(A[i - 1], m, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
            }
            else {
                int dev_null;
                MPI_Send(&dev_null, 0, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        for (int i = 1; i <= n; i++) {
            int sender = (i - 1) % (num_proc - 1) + 1;
            /*
             * 循环接收每一行，如果某一行计算完之后还有行
             * 就将这一行分配，继续计算，否则发送结束信号
             * 需要注意的是，原程序中 MPI_Send(A[num_send - 1], m, MPI_DOUBLE, sender, 99, MPI_COMM_WORLD); 有误
             * 在上一阶段的循环中 num_send - 1 及之前的部分全部计算完成，因此要从 num_send 开始计算
             *
             * 源码中的问题4：
             *     用于计算的行数错误，将导致计算结果出错
             */
            MPI_Recv(C[i - 1], p, MPI_DOUBLE, sender, 100, MPI_COMM_WORLD, &status);
            if (num_send < n) {
                MPI_Send(A[num_send], m, MPI_DOUBLE, sender, 99, MPI_COMM_WORLD);
                num_send++;
            }
            else {
                int dev_null;
                MPI_Send(&dev_null, 0, MPI_INT, sender, 0, MPI_COMM_WORLD);
            }
        }
        /*
         * 计算完成，停止计时(和串行程序一样，不计算 IO 时间)
         * 分别输出结果和时间信息到两个文件
         */
        double finish_time = MPI_Wtime();
        write_time("run_time.stdout", finish_time - start_time);
        write_file("result.stdout", n, p, LEN_P, C[0]);
    }
    else {
        /*
         * 根据处理器核心数设置线程数
         * 为每个线程分配存储线程号的空间
         * 为用于计算的行和存储结果的行分配空间
         * 为传给各个线程的参数分配空间
         */
        int num_threads = get_nprocs();
        pthread_t* thread_ids = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
        double* A_row = (double*)malloc(m * sizeof(double));
        double* C_row = (double*)malloc(p * sizeof(double));
        struct threadArg* thread_args = (struct threadArg *)malloc(num_threads * sizeof(struct threadArg));
        for (int i = 0; i < num_threads; i++) {
            /*
             * 初始化线程 ID 用于划分任务
             * 初始化总线程个数用于划分任务
             * 传入的参数 m 和 p 在此处设置
             */
            thread_args[i].tid = i;
            thread_args[i].len_m = m;
            thread_args[i].len_p = p;
            thread_args[i].A_row = A_row;
            thread_args[i].C_row = C_row;
            thread_args[i].num_threads = num_threads;
        }
        while (1) {
            /*
             * 循环接收消息，若收到结束信号(即 status.MPI_TAG == 0 时)则停止
             * 否则创建若干线程计算出 C_row 并返回
             */
            MPI_Recv(A_row, m, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == 0) {
                break;
            }
            /*
             * 先统一创建各个线程并开始执行
             * 再等待各个线程全部结束
             * 处理完成，将结果发送回 0 号节点
             */
            for (int i = 0; i < num_threads; i++) {
                pthread_create(thread_ids + i, NULL, worker, thread_args + i);
            }
            for (int i = 0; i < num_threads; i++) {
                pthread_join(thread_ids[i], NULL);
            }
            MPI_Send(C_row, p, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
        }
    }
    /*
     * 结束 MPI 环境并退出程序
     */
    MPI_Finalize();
    return 0;
}