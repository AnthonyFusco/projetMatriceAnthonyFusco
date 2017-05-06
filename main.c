#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <mpi.h>

void mpi() {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    //omp_set_num_threads(4);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    int number;
    if(world_rank == 0) {
        /*for(int i = 1; i < world_size; i++) {
            printf("send %d\n", world_rank);
            number = i * 10;
            MPI_Send(&number, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }*/

// Create a buffer that will hold a subset of the random numbers
        // float *sub_rand_nums = malloc(sizeof(float) * elements_per_proc);

// Scatter the random numbers to all processes
        //MPI_Scatter(rand_nums, elements_per_proc, MPI_FLOAT, sub_rand_nums,
        //            elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    number = -1;
    if(world_rank != 0) {
        //printf("Goodbye from rank %d\n", world_rank);
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        printf("Process %d received number %d from process 0\n",
               world_rank, number);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}

typedef struct Matrix Matrix;
struct Matrix {
    long lines;
    long columns;
    int* tab;
};

void initMatrix(Matrix* matrix, long lines, long column) {
    long size = lines * column;
    matrix->lines = lines;
    matrix->columns = column;
    matrix->tab = malloc(sizeof(int) * size);
    for (int i = 0; i < size; ++i) {
        matrix->tab[i] = i;
    }
}

int get(Matrix* matrix, long line, long column) {
    assert(line < matrix->lines);
    assert(column < matrix->columns);
    return matrix->tab[column + line * matrix->columns];
}

void set(Matrix* matrix, long line, long column, int value) {
    assert(line < matrix->lines);
    assert(column < matrix->columns);
    matrix->tab[column + line * matrix->columns] = value;
}

long countElementOfFile(char* fileName) {
    int elem = 0;
    long nbElem = 0;
    FILE *file;
    file = fopen(fileName, "r");
    if (file) {
        while (fscanf(file, "%d", &elem) > 0) {
            nbElem++;
        }
        fclose(file);
    }
    return nbElem;
}

void readMatrixFromFile(Matrix* matrix, char* fileName) {
    long nbElement = countElementOfFile(fileName);
    long lines = (long) sqrt(nbElement);
    long columns = (long) sqrt(nbElement);
    initMatrix(matrix, lines, columns);
    int elem;
    FILE *file;
    file = fopen(fileName, "r");
    if (file) {
        int i = 0;
        while (fscanf(file, "%d", &elem) > 0) {
            matrix->tab[i] = elem;
            i++;
        }
        fclose(file);
    }
}

void printMatrix(Matrix* matrix) {
    for (int i = 0; i < matrix->lines; ++i) {
        for (int j = 0; j < matrix->columns; ++j) {
            printf("%d", get(matrix, i, j));
            if (j != matrix->columns - 1) {
                printf(" ");
            } else {
                printf("\n");
            }
        }
    }
    printf("\n");
}

void turnMatrix(Matrix* matrix, Matrix* turned) {
    turned->columns = matrix->lines;
    turned->lines = matrix->columns;
    for(int i = 0; i < matrix->lines; i++){
        for(int j = 0; j < matrix->columns; j++) {
            set(turned, j, i, get(matrix, i, j));
        }
    }
}

void dispatchMatrix(long n){

// Create a buffer that will hold a subset of the random numbers
   // float *sub_rand_nums = malloc(sizeof(float) * elements_per_proc);

// Scatter the random numbers to all processes
    //MPI_Scatter(rand_nums, elements_per_proc, MPI_FLOAT, sub_rand_nums,
    //            elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void initMPI() {
    MPI_Init(NULL, NULL);
    //omp_set_num_threads(4);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
}
/*
int* scatter(long* size, Matrix* matrix, int world_size) {
    MPI_Bcast(size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    int numberElement = (int) ((*size * *size) / world_size);
    int* buffer = malloc(sizeof(int) * numberElement);
    MPI_Scatter(matrix->tab, numberElement, MPI_INT, buffer,
                numberElement, MPI_INT, 0, MPI_COMM_WORLD);
    return buffer;
}
*/

Matrix scatter(long* size, Matrix* matrix, int world_size) {
    MPI_Bcast(size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    int numberElement = (int) ((*size * *size) / world_size);
    int* buffer = malloc(sizeof(int) * numberElement);
    MPI_Scatter(matrix->tab, numberElement, MPI_INT, buffer,
                numberElement, MPI_INT, 0, MPI_COMM_WORLD);
    Matrix result;
    result.lines = *size / world_size;
    result.columns = *size;
    result.tab = buffer;
    return result;
}

int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);
    //omp_set_num_threads(4);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    Matrix A;
    Matrix BTurned;
    long sizeA;
    long sizeB;

    if (world_rank == 0) {
        char* AFileName = argv[1];
        readMatrixFromFile(&A, AFileName);
        //printMatrix(&A);
        char *BFileName = argv[2];
        Matrix B;
        readMatrixFromFile(&B, BFileName);
        //printMatrix(&B);
        Matrix turned;
        turnMatrix(&B, &turned);
        //printMatrix(&turned);
        sizeA = A.lines;
        sizeB = turned.lines;
        BTurned.tab = turned.tab;
        BTurned.lines = turned.lines;
        BTurned.columns = turned.columns;
    }

    /*MPI_Bcast(&sizeA, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sizeB, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    int numberElementA = (int) ((sizeA * sizeA) / world_size);
    int* bufferA = malloc(sizeof(int) * numberElementA);
    MPI_Scatter(A.tab, numberElementA, MPI_INT, bufferA,
                numberElementA, MPI_INT, 0, MPI_COMM_WORLD);

    int numberElementB = (int) ((sizeB * sizeB) / world_size);
    int* bufferB = malloc(sizeof(int) * numberElementB);
    MPI_Scatter(BTurned.tab, numberElementB, MPI_INT, bufferB,
                numberElementB, MPI_INT, 0, MPI_COMM_WORLD);*/

    Matrix receivedA = scatter(&sizeA, &A, world_size);
    Matrix receivedB = scatter(&sizeB, &BTurned, world_size);

   /* Matrix receivedA;
    receivedA.lines = sizeA / world_size;
    receivedA.columns = sizeA;
    receivedA.tab = bufferA;*/
    printf("I am proc %d, %d\n", world_rank, (int) sizeA);
    printMatrix(&receivedA);

    /*Matrix receivedB;
    receivedB.lines = sizeB / world_size;
    receivedB.columns = sizeB;
    receivedB.tab = bufferB;*/
    printf("I am proc %d, %d\n", world_rank, (int) sizeB);
    printMatrix(&receivedB);


    //initb

    //calc

    //gather

    // Finalize the MPI environment.
    MPI_Finalize();
}


