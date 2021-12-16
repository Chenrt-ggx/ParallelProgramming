#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * 定义不同的 tag
 * 用以区分进程通信时的消息
 */
#define MULTI_LEN 600
#define MULTI_TYPE 300
#define OUTPUT_TYPE 200
#define ALLTOONE_TYPE 100

/* 输出错误信息 */
void merror(char* s) {
    printf("%s\n", s);
    exit(1);
}

/* 交换两个变量 */
void swap(int* x, int* y) {
    int t = *x;
    *x = *y, *y = t;
}

/* 串行快速排序算法 */
void quick_sort(int data[], int l, int r) {
	if (l > r) return;
	int mid = (l + r)  / 2;
	/*
	 * 对源程序进行了修改
	 * 通过三点取中的方式避免被基本有序的数据卡
	 */
	if (data[l] > data[r]) swap(data + l, data + r);
	if (data[r] < data[mid]) swap(data + r, data + mid);
	if (data[l] < data[mid]) swap(data + l, data + mid);
	int i = l, j = r, temp = data[l];
	/* 完成对基准元素的归位 */
	while (i - j) {
		while (data[j] >= temp && i < j) j--;
		while (data[i] <= temp && i < j) i++;
		if (i < j) swap(data + i, data + j);
	}
    data[l] = data[i], data[i] = temp, mid = j = i;
	/*
	 * 将直接在基准元素两侧递归改成两侧去重后递归
	 * 不进行此修改时会被全部相同的数据卡
	 */
	while (i >= l && data[i] == data[mid]) i--;
	while (j <= r && data[j] == data[mid]) j++;
	/* 在基准元素两侧递归快速排序 */
    quick_sort(data, l, i);
    quick_sort(data, j, r);
}

/* 合并两个有序数组 */
void merge(const int src[], int dst[], int l, int mid, int r) {
    /*
     * 修改了输入格式并重命名了变量
     * 原来的 s1 对应现在的 mid
     * 原来的 s1 + s2 对应现在的 r
     * 额外传入参数 l 使函数更具有一般性
     */
    int index = 0, l_ptr = l, r_ptr = mid;
    while (l_ptr < mid && r_ptr < r) {
        if (src[l_ptr] > src[r_ptr]) dst[index++] = src[r_ptr++];
        else dst[index++] = src[l_ptr++];
    }
    /*
     * 将原本的 for 循环按照 l_ptr < mid && r_ptr < r拆开
     * 削弱了循环内部分支的复杂程度
     * 有利于编译器产生更高效的汇编代码
     */
    while (l_ptr < mid) dst[index++] = src[l_ptr++];
    while (r_ptr < r) dst[index++] = src[r_ptr++];
}

/* 串行多路归并算法递归部分 */
void recursion_merge(int src[], int buf[], int index[], int l, int r) {
    if (r - l < 1) return;
    int mid = (l + r) / 2;
    /*
     * 分别递归两边，得到两个有序子数组
     * 再合并得到的两个有序子数组
     */
    recursion_merge(src, buf, index, l, mid);
    recursion_merge(src, buf, index, mid + 1, r);
    /*
     * 多路归并的关键在于“分治”
     * 对于 [1], [2], [3], [4], [5], [6], [7], [8] 这八个数组，合并的方式是：
     * [1], [2], [3], [4], [5], [6], [7], [8]
     * [1, 2], [3, 4], [5, 6], [7, 8]
     * [1, 2, 3, 4], [5, 6, 7, 8]
     * [1, 2, 3, 4, 5, 6, 7, 8]
     * 而不是：
     * [1, 2], [3], [4], [5], [6], [7], [8]
     * [1, 2, 3], [4], [5], [6], [7], [8]
     * [1, 2, 3, 4], [5], [6], [7], [8]
     * [1, 2, 3, 4, 5], [6], [7], [8]
     * [1, 2, 3, 4, 5, 6], [7], [8]
     * [1, 2, 3, 4, 5, 6, 7], [8]
     * [1, 2, 3, 4, 5, 6, 7, 8]
     * 前者的复杂度是 (n / p) * p * log(p) = n * log(p) 的
     * 后者的复杂度是 (n / p) * p * p = n * p 的
     */
    merge(src, buf, index[l], index[mid + 1], index[r + 1]);
    /*
     * 不采用原实现中的“滚动数组”
     * 两个数组一个只作为原数组使用，一个只作为 buffer 使用
     * 不会出现原实现中需要根据奇偶性判断最终结果在哪一个数组的情况
     * 因此在完成一次归并后，需要更新 src 数组对应的位置
     */
    for (int i = index[l]; i < index[r + 1]; i++) {
        src[i] = buf[i - index[l]];
    }
    /*
     * 此外，多路归并算法还可以通过堆实现，亦有 n * log(p) 的复杂度
     * 相比分治实现，堆实现有更小的函数调用开销，但空间连续性不如分治实现
     * 在不同的机器上运行，两种方式可能有不同的效果
     */
}

/* 串行多路归并算法递归前部分 */
void multi_merge(int src[], int buf[], const int length[], int total) {
    int* index = (int*)malloc((total + 1) * sizeof(int));
    index[0] = 0;
    /*
     * 输入的 length 意义是每一段子数组的长度
     * 为了在递归部分中快速得到每一段子数组的下标信息，需对 length 数组求和
     * 本质上是前缀和优化，实现 O(1) 得到下标
     */
    for (int i = 1; i <= total; i++) {
        index[i] = index[i - 1] + length[i - 1];
    }
    /*
     * 原实现本质上还是分治实现，不过其尾递归实现更像是“披着递归皮的循环”
     * 且其中存在较多特殊处理，如对 length 数组压缩和清零，对区间数的奇偶特判等
     * 考虑到其原实现不够清晰优雅，对此部分进行了重写
     */
    recursion_merge(src, buf, index, 0, total - 1);
    free(index);
}

/* 并行正则采样排序的主体 */
void psrs_main(int data_size) {
    int id, total, mutex;
    MPI_Status status[32 * 32 * 2];
    MPI_Request request[32 * 32 * 2];
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &total);

    /*
     * 申请空间并完成初始化 split 和 length
     * 此处采用一次性 malloc 出 data 和 buffer 的方式
     * 在归并时避免了重复 malloc 出 buffer 造成的性能下降
     *
     * split 是分段的数量
     * length 是此处理器处理的数据量
     * data 是长度为 data_size 的数组，存原始数据
     * buffer 是长度为 data_size 的数组，作为缓冲区
     */
    int split = total - 1;
    int length = data_size / total;
    int* data = (int*)malloc(data_size * sizeof(int) * 2);
    if (data == 0) merror("malloc memory for array error!");
    int* buffer = data + data_size;

    /*
     * 每个处理器上生成随机数用于排序
     * 原程序中输出是乱序的，需要进行处理，方法和第一次作业类似
     */
    if (id != 0) {
        MPI_Recv(&mutex, 1, MPI_CHAR, id - 1,
                 OUTPUT_TYPE, MPI_COMM_WORLD, status);
    }
    /*
     * 用自身 id 作为随机数种子
     * 确保每次生成的数据相同，便与调试
     */
    srand(id);
    /*
     * 需要注意的是，通过 MPI_Recv 和 MPI_Send 实现的串行输出只保证调用顺序
     * 而在 mpirun 中，多个处理器上的 io 顺序并不是调用顺序
     * 因此即使按顺序调用了 printf 输出顺序也有可能是错误的
     * 需要保证输出有序，可将数据发送到一个处理器，统一输出
     *
     * 此外，经过实验可以发现 mpirun 会倾向于优先满足短输出
     * 且 mpirun 的 io 是非抢占的，即在一个 printf 执行的过程中，它不会被打断
     *
     * 在保证调用顺序后，此处的输出往往是有序的
     * 这可能是因为各个处理器上的输出量差别不大，因此 mpirun 按照调用顺序输出
     */
    printf("----> This is node %d \n", id);
    printf("On node %d the input data is:\n", id);
    for (int i = 0; i < length; i++) {
        data[i] = (int)random() % 0xffff - 0x8000;
        printf("%d ", data[i]);
    }
    printf("\n");
    if (id != total - 1) {
        MPI_Send(&mutex, 1, MPI_CHAR, id + 1,
                 OUTPUT_TYPE, MPI_COMM_WORLD);
    }

    /*
     * 已清除部分没有必要的 MPI_Barrier 以提高性能
     * 每个处理器将自己的 n / P 个数据串行快速排序
     * 对应于算法 13.5 步骤 ( 1 )
     */
    quick_sort(data, 0, length - 1);
    MPI_Barrier(MPI_COMM_WORLD);

    // ------------------------------------------------------------------

    if (split) {
        /*
         * 处理器不止一个，此时有并行正则采样的过程
         * 申请数组 sample 存放自身选取的 split 个元素
         * 申请数组 index 存放此处理器的开始位点和结束位点
         */
        int* sample = (int*)malloc(total * split * sizeof(int));
        if (sample == 0) merror("malloc memory for sample error!");
        int* index = (int*)malloc(total * sizeof(int) * 2);
        if (index == 0) merror("malloc memory for index error!");

        /*
         * 每个处理器从排好序的序列中选取第 w, 2w, 3w, …, (P - 1)w
         * 共 P - 1 个数据作为代表元素，其中 w = n / P * P
         * 对应于算法 13.5 步骤 ( 2 )
         */
        for (int i = 0, n = length / total; i < split; i++)
            sample[i] = data[(i + 1) * n - 1];
        MPI_Barrier(MPI_COMM_WORLD);

        if (id == 0) {
            /*
             * 每个处理器将选好的代表元素送到处理器 P0 中
             * 对应于算法 13.5 步骤 ( 3 )
             */
            for (int i = 1, j = 0; i < total; i++) {
                /*
                 * MPI_Irecv 是非阻塞接受，可在没接收完的时候继续接受新的
                 * int MPI_Irecv(
                 *     void *buf,
                 *     int count,
                 *     MPI_Datatype datatype,
                 *     int source,
                 *     int tag,
                 *     MPI_Comm comm,
                 *     MPI_Request *request
                 * )
                 * 用 request 保存接受结构体，在 MPI_Waitall 中判断是否完成接收
                 */
                MPI_Irecv(sample + i * split, (int) (split * sizeof(int)), MPI_CHAR,
                          i, ALLTOONE_TYPE + i, MPI_COMM_WORLD, request + j++);
            }
            /* 等待全部数据完成传送 */
            MPI_Waitall(split, request, status);
            /* 第一次同步：发送完成，0 号处理器接收完成 */
            MPI_Barrier(MPI_COMM_WORLD);

            /*
             * 处理器 P0 将上一步送来的 P 段有序的数据序列做 P 路归并
             * 再选择排序后的第 P - 1, 2 * (P - 1), …, (P - 1) * (P - 1) 个共 P - 1 个主元
             * 对应于算法 13.5 步骤 ( 3 )
             *
             * 这里改为采用串行快速排序
             * 一方面快速排序也可以使数据有序，供接下来的过程使用
             * 另一方面，排序数据的个数是 P * P 的量级，相对 data_size 不大
             * 因此使用快速排序也有较高的效率
             */
            quick_sort(sample, 0, total * split - 1);
            /* 第二次同步：0 号处理器排序完成 */
            MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 1; i < total; i++)
                sample[i] = sample[i * split - 1];

            /*
             * 处理器 P0 将这 P - 1 个主元播送到所有处理器中
             * 对应于算法 13.5 步骤 ( 4 )
             */
            MPI_Bcast(sample, (int)(total * sizeof(int)), MPI_CHAR, 0, MPI_COMM_WORLD);
            /* 第三次同步：0 号处理器广播，其它处理器接收完成 */
            MPI_Barrier(MPI_COMM_WORLD);
        }
        else {
            /* 将当前处理器选出的代表元素发送到 0 号处理器 */
            MPI_Send(sample, (int)(split * sizeof(int)), MPI_CHAR, 0, ALLTOONE_TYPE + id, MPI_COMM_WORLD);
            /* 第一次同步：发送完成，0 号处理器接收完成 */
            MPI_Barrier(MPI_COMM_WORLD);
            /* 第二次同步：0 号处理器排序完成 */
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(sample, (int)(total * sizeof(int)), MPI_CHAR, 0, MPI_COMM_WORLD);
            /* 第三次同步：0 号处理器广播，其它处理器接收完成 */
            MPI_Barrier(MPI_COMM_WORLD);
        }

        /*
         * 每个处理器根据上步送来的 P - 1 个主元把自己的 n / P 个数据分成 P 段
         * 记为处理器 Pi 的第 j + 1 段，其中 i = 0, …, P - 1且 j = 0, …, P - 1
         * 对应于算法 13.5 步骤 ( 5 )
         */
        int current = 1;
        /* 对小于首个元素的主元，将它的开始位点和结束位点清零 */
        for (index[0] = 0; data[0] >= sample[current] && current < total; current++)
            index[2 * current] = index[2 * current - 1] = 0;
        /* 如果所有元素大于最大主元，直接将所有元素放在最后一段 */
        if (current == total) index[2 * current - 1] = length;
        for (int l = 0; current < total; current++) {
            /*
             * 原程序变量命名过于阴间
             * 已全部重命名，以便望文生义
             * 二分查找出比当前主元大的最小元素的位置存入 mid
             *
             * 一个更好的处理是让 data[l] <= sample[current]，data[r] > sample[current]
             * 在此基础上二分缩小 l 和 r 直到 r = l + 1，此时 r 即为所求
             * 按此方法处理无需额外的循环对 mid 进行修正
             */
            int r = length, mid = (l + r) / 2;
            while (data[mid] != sample[current] && l < r) {
                if (data[mid] > sample[current]) r = mid - 1;
                else l = mid + 1;
                mid = (l + r) / 2;
            }
            while (data[mid] <= sample[current] && mid < length) mid++;
            if (mid == length) {
                /*
                 * 当前段是最后一段，更新结束位点
                 * 后面的各段将不再有意义，将它们的下标开始位点和结束位点清零
                 * 最后更新 current 为 total 跳出循环
                 */
                index[2 * current - 1] = length;
                for (int i = current; i < total; i++)
                    index[2 * i] = index[2 * i + 1] = 0;
                current = total;
            }
            else {
                /*
                 * 上一段的结束位点是 mid
                 * 下一段的开始位点是 mid
                 */
                index[2 * current] = mid;
                index[2 * current - 1] = mid;
            }
            /* 更新二分开始的位置，进行新一轮二分 */
            l = mid;
        }
        /* 更新最后一段的结束位点 */
        if (current == total) index[2 * current - 1] = length;
        /* 各处理器完成分段 */
        MPI_Barrier(MPI_COMM_WORLD);

        // ------------------------------------------------------------------

        /*
         * 每个处理器发送它的第 i + 1 段给处理器 Pi
         * 从而使得第 i 个处理器含有所有处理器的第 i 段数据 (i = 0, …, P - 1)
         * 对应于算法 13.5 步骤 ( 6 )
         */
        for (int i = 0, j = 0; i < total; i++) {
            if (i == id) {
                /*
                 * 是本处理器，需要发送数据
                 * 此处复用了 sample 数组，存放各段的长度
                 * 自己的 sample[i] 直接计算即可，不需要通过 MPI_Recv 得到
                 */
                sample[i] = index[2 * i + 1] - index[2 * i];
                for (int k = 0; k < total; k++)
                    if (k != id) {
                        /* 向第 k 个处理器发送本机第 k + 1 段的长度 */
                        int l = index[2 * k + 1] - index[2 * k];
                        MPI_Send(&l, sizeof(int), MPI_CHAR, k, MULTI_LEN + id, MPI_COMM_WORLD);
                    }
            }
            else {
                /* 从第 i 个处理器获取其第 i + 1 段的长度，存入 sample[i] */
                MPI_Recv(sample + i, sizeof(int), MPI_CHAR, i,
                        MULTI_LEN + i, MPI_COMM_WORLD, status + j++);
            }
        }
        /* 发送完成，进行同步 */
        MPI_Barrier(MPI_COMM_WORLD);

        // ------------------------------------------------------------------

        /*
         * 每个处理器向其它处理器发送对应的数据
         * 第 i 个处理器向第 j 个处理器发送第 j 段
         * 用 length 维护当前处理器接受的数据的数量
         */
        length = 0;
        for (int i = 0, j = 0; i < total; i++) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (i == id) {
                /* 当前处理器无需给自己发送，直接复制即可 */
                for (int k = index[2 * i]; k < index[2 * i + 1]; k++)
                    buffer[length++] = data[k];
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (i == id) {
                for (int k = 0; k < total; k++)
                    if (k != id) {
                        /* 向其它处理器 k 发送对应的第 k 段 */
                        int count = index[2 * k + 1] - index[2 * k];
                        MPI_Send(data + index[2 * k], (int)(count * sizeof(int)),
                                 MPI_CHAR, k, MULTI_TYPE + id, MPI_COMM_WORLD);
                    }
            }
            else {
                /* 接收来自其它处理器的第 i 段，并维护 length */
                MPI_Recv(buffer + length, (int)(sample[i] * sizeof(int)), MPI_CHAR,
                         i , MULTI_TYPE + i, MPI_COMM_WORLD, status + j++);
                length = length + sample[i];
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        /*
         * 已清除部分没有必要的 MPI_Barrier 以提高性能
         * 每个处理器再通过 P 路归并排序将上一步的到的数据排序
         * 从而这 n 个数据便是有序的
         * 对应于算法 13.5 步骤 ( 7 )
         *
         * 此处的 data 是实质上的 buffer
         * 此处的 buffer 是实质上的 data
         * 归并的结果最终存在 buffer
         */
        multi_merge(buffer, data, sample, total);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // ------------------------------------------------------------------

    if (id != 0) {
        MPI_Recv(&mutex, 1, MPI_CHAR, id - 1,
                 OUTPUT_TYPE, MPI_COMM_WORLD, status);
    }
    /*
     * 此处同样存在输出顺序问题，通过 MPI_Recv 和 MPI_Send 保证 printf 的调用顺序并不能保证各处理器的输出顺序
     * 各个处理器上的输出量差别较大，会触发 mpirun 优先满足短输出的调度
     * 这显然不便于我们观察结果，例如可能出现 (#) ：
     *     某个处理器正在循环输出经过排序的数字
     *     另一个处理器输出 "On node %d the sorted data is : \n"
     *     此处理器继续循环输出经过排序的数字
     * 为此可以利用 mpirun 中 io 的非抢占性进行优化，具体的：
     *     申请一个缓冲区数组，将输出以字符串的形式存放在缓冲区
     *     一次性输出缓冲区，因为是单次输出，不会被打断
     * 虽然可能仍有输出不按顺序，但不再会出现上述 (#) 的情况，观察体验大大增加
     */
    char* s = (char*)malloc(16 * data_size);
    sprintf(s, "# On node %d the sorted data is : \n", id);
    for (int i = 0; i < length; i++) {
        if (i + 1 != length) sprintf(s, "%s %d ", s, buffer[i]);
        else sprintf(s, "%s %d\n", s, buffer[i]);
    }
    fputs(s, stdout), free(s);
    if (id != total - 1) {
        MPI_Send(&mutex, 1, MPI_CHAR, id + 1,
                 OUTPUT_TYPE, MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[]) {
    /* 当前算法需要保证数据总数时处理器的个数的倍数 */
    int id, data_size = 16 * 64;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (id == 0) printf("Data Size : %d\n\n", data_size);
    psrs_main(data_size);
    /* 输出换行前同步一次，避免非 0 号处理器还没有输出完而换行已输出 */
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) printf("\n");
    MPI_Finalize();
    return 0;
}