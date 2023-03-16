#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<mpi.h>

double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots){
    double* rebuiltCoord = (double*)malloc(sizeof(double) * n * k);
    int i;
    for(i=0; i<n*k; i++){
        rebuiltCoord[i] = 0;
    }

    // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
    for(i=0; i<n; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            double distance = 0;
            int pivoti = pivots[ki];
            int j;
            for(j=0; j<dim; j++){
                distance += pow(coord[pivoti*dim + j] - coord[i*dim + j] ,2);
            }
            rebuiltCoord[i*k + ki] = sqrt(distance);
        }
    }

    // Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
    double chebyshevSum = 0;
    for(i=0; i<n; i++){
        int j;
        for(j=0; j<n; j++){
            double chebyshev = 0;
            int ki;
            for(ki=0; ki<k; ki++){
                double dis = fabs(rebuiltCoord[i*k + ki] - rebuiltCoord[j*k + ki]);
                chebyshev = dis>chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
    }

    free(rebuiltCoord);

    return chebyshevSum;
}

//++++++++++++++将任意两点之间的距离求出来
double getDistance(const int dim, int a, int b, double* coord) {
    double distance = 0;
    for (int i = 0; i < dim; i++) {
        distance += pow(coord[a * dim + i] - coord[b * dim + i], 2);
    }
    return sqrt(distance);
}
//++++++++组合数计算
long C(int n,int k){
    long ans = 1;
    for(int i =1;i<=k;i++){
        ans = ans * (n-k+i)/i; // 注意一定要先乘再除
    }
    return ans;
}

//long C(int n,int k){
//    long sum = n;
//    int divide = 1;
//    for(int i=1;i<k;i++){
//        sum *= n-i;
//        divide *= i+1;
//    }
//    return sum/divide;
//}

//++++++++++++++++++新加的参数为N和my_rank
void Combination(int &myassignment,int ki, const int k, const int n, const int dim, const int M, double* coord, int* pivots,
                 double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots) {
    if (ki == k - 1) {
        int i;
        for (i = pivots[ki - 1] + 1; i < n; i++) {
            pivots[ki] = i;

            double distanceSum = SumDistance(k, n, dim, coord, pivots);
            maxDistanceSum[M] = distanceSum;
            minDistanceSum[M] = distanceSum;
            int kj;
            for (kj = 0; kj < k; kj++) {
                maxDisSumPivots[M * k + kj] = pivots[kj];
            }
            for (kj = 0; kj < k; kj++) {
                minDisSumPivots[M * k + kj] = pivots[kj];
            }
            // sort
            int a;
            for (a = M; a > 0; a--) {
                if (maxDistanceSum[a] > maxDistanceSum[a - 1]) {
                    double temp = maxDistanceSum[a];
                    maxDistanceSum[a] = maxDistanceSum[a - 1];
                    maxDistanceSum[a - 1] = temp;
                    int kj;
                    for (kj = 0; kj < k; kj++) {
                        int temp = maxDisSumPivots[a * k + kj];
                        maxDisSumPivots[a * k + kj] = maxDisSumPivots[(a - 1) * k + kj];
                        maxDisSumPivots[(a - 1) * k + kj] = temp;
                    }
                }
            }
            for (a = M; a > 0; a--) {
                if (minDistanceSum[a] < minDistanceSum[a - 1]) {
                    double temp = minDistanceSum[a];
                    minDistanceSum[a] = minDistanceSum[a - 1];
                    minDistanceSum[a - 1] = temp;
                    int kj;
                    for (kj = 0; kj < k; kj++) {
                        int temp = minDisSumPivots[a * k + kj];
                        minDisSumPivots[a * k + kj] = minDisSumPivots[(a - 1) * k + kj];
                        minDisSumPivots[(a - 1) * k + kj] = temp;
                    }
                }
            }
            myassignment--;
            if(myassignment==0){
                return;
            }
        }
        return;
    }
    if(ki==0){
        int i;
        for (i = pivots[0]; i < n; i++) {
            pivots[ki] = i;
            Combination(myassignment,ki + 1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);
            if(myassignment==0){
                return;
            }
        }
    }else{
        int i;
        for (i = pivots[ki - 1] + 1; i < n; i++) {
            pivots[ki] = i;
            Combination(myassignment,ki + 1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);
            if(myassignment==0){
                return;
            }
        }
    }


}

//++++++
void Combination1(const int k, const int n,const int dim,const int M,double* coord,int* pivots, double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots){
    double distanceSum = SumDistance(k, n, dim, coord, pivots);
    maxDistanceSum[M] = distanceSum;
    minDistanceSum[M] = distanceSum;
    int kj;
    for (kj = 0; kj < k; kj++) {
        maxDisSumPivots[M * k + kj] = pivots[kj];
    }
    for (kj = 0; kj < k; kj++) {
        minDisSumPivots[M * k + kj] = pivots[kj];
    }
    // sort
    int a;
    for (a = M; a > 0; a--) {
        if (maxDistanceSum[a] > maxDistanceSum[a - 1]) {
            double temp = maxDistanceSum[a];
            maxDistanceSum[a] = maxDistanceSum[a - 1];
            maxDistanceSum[a - 1] = temp;
            int kj;
            for (kj = 0; kj < k; kj++) {
                int temp = maxDisSumPivots[a * k + kj];
                maxDisSumPivots[a * k + kj] = maxDisSumPivots[(a - 1) * k + kj];
                maxDisSumPivots[(a - 1) * k + kj] = temp;
            }
        }
    }
    for (a = M; a > 0; a--) {
        if (minDistanceSum[a] < minDistanceSum[a - 1]) {
            double temp = minDistanceSum[a];
            minDistanceSum[a] = minDistanceSum[a - 1];
            minDistanceSum[a - 1] = temp;
            int kj;
            for (kj = 0; kj < k; kj++) {
                int temp = minDisSumPivots[a * k + kj];
                minDisSumPivots[a * k + kj] = minDisSumPivots[(a - 1) * k + kj];
                minDisSumPivots[(a - 1) * k + kj] = temp;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    // filename : input file namespace
    char* filename = (char*)"../uniformvector-2dim-5h.txt";
    if (argc == 2) {
        filename = argv[1];
    }
    else if (argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    // M : number of combinations to store
    const int M = 1000;
    // dim : dimension of metric space
    int dim;
    // n : number of points
    int n;
    // k : number of pivots
    int k;
    // Read parameter
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);
    fscanf(file, "%d", &n);
    fscanf(file, "%d", &k);
//    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Start timing
    struct timeval start;

    // Read Data
    double* coord = (double*)malloc(sizeof(double) * dim * n);
    int i;
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < dim; j++) {
            fscanf(file, "%lf", &coord[i * dim + j]);
        }
    }
    fclose(file);
    gettimeofday(&start, NULL);

    //++++++++++++++++++++++
    //N总进程数量
    //my_rank当前进程编号；
    int comm_sz;
    int my_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //++++++++++++++++++
    double* maxDistanceSum = (double*)malloc(sizeof(double) * 2 * (M + 1));
    for (i = 0; i < M; i++) {
        maxDistanceSum[i] = 0;
    }
    // maxDisSumPivots : the top M pivots combinations
    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * 2 * (M + 1));
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k; ki++) {
            maxDisSumPivots[i * k + ki] = 0;
        }
    }
    // minDistanceSum : the smallest M distance sum
    double* minDistanceSum = (double*)malloc(sizeof(double) * 2 * (M + 1));
    for (i = 0; i < M; i++) {
        minDistanceSum[i] = __DBL_MAX__;
    }
    // minDisSumPivots : the bottom M pivots combinations
    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * 2 * (M + 1));
    for (i = 0; i < (2 * M + 1); i++) {
        int ki;
        for (ki = 0; ki < k; ki++) {
            minDisSumPivots[i * k + ki] = 0;
        }
    }

//    //初始化距离表格
//    int local_n = n / comm_sz;
//    int local_a = my_rank * local_n;
//    int local_b = (my_rank + 1) * local_n;
//    double* reCoord = (double*)malloc(sizeof(double) * n * n),
//            * local_reCoord = (double*)malloc(sizeof(double) * local_n * n);
//    int endi = n * n, endj = local_n * n;
//    for (int i = 0; i < endj; i++)
//        local_reCoord[i] = 0;
//    for (i = 0; i < endi; i++)
//        reCoord[i] = 0;
//
//    for (int i = local_a, s = 0; i < local_b; i++, s++)
//    {
//        for (int j = i, t = i; j < n; j++, t++)
//        {
//            local_reCoord[s * n + t] = getDistance(dim,i,j,coord);
//        }
//    }
//    MPI_Allgather(local_reCoord, local_n * n, MPI_DOUBLE, reCoord, local_n * n, MPI_DOUBLE, MPI_COMM_WORLD);
//    for(int i = 0; i < n; i++)
//    {
//        for(int j = 0; j < n; j++)
//        {
//            reCoord[j * n + i] = reCoord[i * n + j];
//        }
//    }
//    //！！！！！！！！！！！！！！表格距离信息需要通信

    //++++++++++++去掉递归，均分进程任务
    int *myEndPivots = new int[k];
    int *pivots = new int[k];
    int myassignment = C(n,k)/comm_sz;
    int assignment = my_rank * myassignment + 1;
    if(my_rank==comm_sz-1){
        myassignment = myassignment + C(n,k)%comm_sz;
    }
//    cout<<"assignment = "<<assignment<<endl;
    int prvious = 0;
    for(int i=0;i<k;i++){
        int count=1;
        int reduce = C(n-prvious-count,k-1-i);
        while(assignment>reduce){
//            cout<<"reduce = "<<reduce<<endl;
            assignment -= reduce;
            count++;
            reduce = C(n-prvious-count,k-1-i);
        }
//        cout<<"reduce = "<<reduce<<endl;
        myEndPivots[i] = count+prvious;//下标从零开始
        pivots[i] = myEndPivots[i]-1;
//        cout<<"i = "<<i<<" myEndPivots[i] = "<<myEndPivots[i]<<endl;
        prvious = myEndPivots[i];
    }
//    for(int kj=0;kj<k;kj++){
//        cout<<myEndPivots[kj]<<" ";
//    }
//    cout<<endl;

    //开始计算
    for(int i = myEndPivots[0]-1;i<n-k+1;i++){
        pivots[0]=i;
        if(i != myEndPivots[0]-1){
            for(int kj = 1;kj < k;kj++){
                pivots[kj]=kj+i;
            }
        }
        bool flag = true;
        while(flag){
            flag = false;
            int start = pivots[k-1];
            for(int j = start;j < n;j++){
                //输出
                pivots[k-1] = j;
//                for(int kj = 0;kj < k;kj++){
//                    cout<<pivots[kj]<<" ";
//                }
//                cout<<endl;
                if(myassignment>0){
                    Combination1(k,n,dim,M,coord,pivots,maxDistanceSum,maxDisSumPivots,minDistanceSum,minDisSumPivots);
                    myassignment--;
                } else{
                    goto myend;
                }
            }
            for(int kj=k-1-1;kj>0;kj--){
                if(pivots[kj] != n-k+kj){
                    pivots[kj]++;
                    for(int ki = kj+1;ki<k;ki++){
                        pivots[ki] = pivots[kj]+ki-kj;
                    }
                    flag = true;
                    break;
                }
            }
        }
    }
    myend:
    //输出查看
//    printf("my_rank = %d ",my_rank);
//    printf("my_rank = %d,my_max = %.6lf,my_min = %.6lf\n",my_rank,maxDistanceSum[0],minDistanceSum[0]);
//    printf("max:");
//    for(int kj=0;kj<k;kj++){
//        printf("%d ",maxDisSumPivots[kj]);
//    }
//    printf("\nmin:");
//    for(int kj=0;kj<k;kj++){
//        printf("%d ",minDisSumPivots[kj]);
//    }



    bool currentSta = true;//值为true时值在数组的前半部分，否则在后半部分
//***************进行合并
    int remain = comm_sz, sum = my_rank, half, rm, M_1 = M + 1;
    double* maxDisSumTemp = (double*)malloc(sizeof(double) * 2 * (M + 1)),
            * minDisSumTemp = (double*)malloc(sizeof(double) * 2 * (M + 1));
    int* maxDisSumPivotsTemp = (int*)malloc(sizeof(int) * k * 2 * (M + 1)),
            * minDisSumPivotsTemp = (int*)malloc(sizeof(int) * k * 2 * (M + 1));
    bool currentStaTemp = true;
    const int k_m = k;
    int kj;
    int sta = 1;
    while (remain != 1) {
        half = remain / 2;
        rm = remain % 2;
        if (my_rank < half) {
            //printf("recv + %d\n", my_rank);
            MPI_Recv(maxDisSumTemp, 2 * (M + 1), MPI_DOUBLE, my_rank + half + rm, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(maxDisSumPivotsTemp, 2 * k * (M + 1), MPI_INT, my_rank + half + rm, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(minDisSumTemp, 2 * (M + 1), MPI_DOUBLE, my_rank + half + rm, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(minDisSumPivotsTemp, 2 * k * (M + 1), MPI_INT, my_rank + half + rm, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&sta, 1, MPI_INT, my_rank + half + rm, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            currentStaTemp = (bool)sta;
            int i, j;
            if (currentSta && currentStaTemp) {//数据都在前半部分
                i = 0; j = 0;
                for (int p = 0; p < M; p++)
                {
                    if (maxDistanceSum[i] >= maxDisSumTemp[j]) {
                        maxDistanceSum[M + 1 + p] = maxDistanceSum[i];
                        //处理PIVOTS
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[(M + 1 + p) * k + kj] = maxDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        maxDistanceSum[M + 1 + p] = maxDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[(M + 1 + p) * k + kj] = maxDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
                i = 0; j = 0;
                for (int p = 0; p < M; p++)
                {
                    if (minDistanceSum[i] <= minDisSumTemp[j]) {
                        minDistanceSum[M + 1 + p] = minDistanceSum[i];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[(M + 1 + p) * k + kj] = minDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        minDistanceSum[M + 1 + p] = minDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[(M + 1 + p) * k + kj] = minDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
            }
            else if (currentSta && !currentStaTemp) {//接收方数据在前，发送方在后
                i = 0; j = M + 1;
                for (int p = 0; p < M; p++)
                {
                    if (maxDistanceSum[i] >= maxDisSumTemp[j]) {
                        maxDistanceSum[M + 1 + p] = maxDistanceSum[i];
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[(M + 1 + p) * k + kj] = maxDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        maxDistanceSum[M + 1 + p] = maxDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[(M + 1 + p) * k + kj] = maxDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
                i = 0; j = M + 1;
                for (int p = 0; p < M; p++)
                {
                    if (minDistanceSum[i] <= minDisSumTemp[j]) {
                        minDistanceSum[M + 1 + p] = minDistanceSum[i];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[(M + 1 + p) * k + kj] = minDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        minDistanceSum[M + 1 + p] = minDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[(M + 1 + p) * k + kj] = minDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
            }
            else if (!currentSta && currentStaTemp) {//前面数据在后，后面在前
                i = M + 1; j = 0;
                for (int p = 0; p < M; p++)
                {
                    if (maxDistanceSum[i] >= maxDisSumTemp[j]) {
                        maxDistanceSum[p] = maxDistanceSum[i];
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[p * k + kj] = maxDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        maxDistanceSum[p] = maxDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[p * k + kj] = maxDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
                i = M + 1; j = 0;
                for (int p = 0; p < M; p++)
                {
                    if (minDistanceSum[i] <= minDisSumTemp[j]) {
                        minDistanceSum[p] = minDistanceSum[i];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[p * k + kj] = minDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        minDistanceSum[p] = minDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[p * k + kj] = minDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
            }
            else if (!currentSta && !currentStaTemp) {
                i = M + 1; j = M + 1;
                for (int p = 0; p < M; p++)
                {
                    if (maxDistanceSum[i] >= maxDisSumTemp[j]) {
                        maxDistanceSum[p] = maxDistanceSum[i];
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[p * k + kj] = maxDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        maxDistanceSum[p] = maxDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            maxDisSumPivots[p * k + kj] = maxDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
                i = M + 1; j = M + 1;
                for (int p = 0; p < M; p++)
                {
                    if (minDistanceSum[i] <= minDisSumTemp[j]) {
                        minDistanceSum[p] = minDistanceSum[i];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[p * k + kj] = minDisSumPivots[i * k + kj];
                        }
                        i++;
                    }
                    else {
                        minDistanceSum[p] = minDisSumTemp[j];
                        for (kj = 0; kj < k; kj++) {
                            minDisSumPivots[p * k + kj] = minDisSumPivotsTemp[j * k + kj];
                        }
                        j++;
                    }
                }
            }currentSta = !currentSta;
        }
        else if (my_rank >= half + rm && my_rank < remain) {
            if (currentSta) sta = 1; else sta = 0;
            MPI_Send(maxDistanceSum, 2 * (M + 1), MPI_DOUBLE, my_rank - half - rm, 0, MPI_COMM_WORLD);
            MPI_Send(maxDisSumPivots, 2 * k * (M + 1), MPI_INT, my_rank - half - rm, 1, MPI_COMM_WORLD);
            MPI_Send(minDistanceSum, 2 * (M + 1), MPI_DOUBLE, my_rank - half - rm, 2, MPI_COMM_WORLD);
            MPI_Send(minDisSumPivots, 2 * k * (M + 1), MPI_INT, my_rank - half - rm, 3, MPI_COMM_WORLD);
            MPI_Send(&sta, 1, MPI_INT, my_rank - half - rm, 4, MPI_COMM_WORLD);
            MPI_Finalize();
            //printf("send + %d\n", my_rank);
            return 0;
        }
        remain = half + rm;
    }

    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    //+++++++++++++++++++++++++++++++++++
    // Log
    int ki, first;
    printf("max : ");
    if (currentSta)
        first = 0;
    else
        first = M + 1;
    for (ki = 0; ki < k; ki++) {
        //printf("%d ", maxDisSumPivots[first + ki]);
    }
    printf("%lf\n", maxDistanceSum[first]);
    printf("min : ");
    for (ki = 0; ki < k; ki++) {
        //printf("%d ", minDisSumPivots[first + ki]);
    }
    printf("%lf\n", minDistanceSum[first]);
    MPI_Finalize();
    // Store the result
    FILE* out = fopen("result.txt", "w");
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", maxDisSumPivots[(first + i) * k + ki]);//后期可把for循环改为first作为变量
        }
        fprintf(out, "%d\n", maxDisSumPivots[(first + i) * k + k - 1]);
    }
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", minDisSumPivots[(first + i) * k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[(first + i) * k + k - 1]);
    }
    fclose(out);

    return 0;
    //printf("my_rank + %f\n" , minDistanceSum[0]);
    //printf("my_rank + %f\n" , minDistanceSum[M+1]);
    MPI_Finalize(); return 0;
}
