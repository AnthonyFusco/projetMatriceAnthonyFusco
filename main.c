#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <mpi.h>

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
    turned->tab = malloc(sizeof(int) * matrix->lines * matrix->lines);
    turned->columns = matrix->lines;
    turned->lines = matrix->columns;
    #pragma omp parallel for
    for(int i = 0; i < matrix->lines; i++){
        #pragma omp parallel for
        for(int j = 0; j < matrix->columns; j++) {
            set(turned, j, i, get(matrix, i, j));
        }
    }
}

Matrix scatter(long size, Matrix* matrix, int world_size) {
    int numberElement = (int) ((size * size) / world_size);
    int* buffer = malloc(sizeof(int) * numberElement);
    MPI_Scatter(matrix->tab, numberElement, MPI_INT, buffer, numberElement, MPI_INT, 0, MPI_COMM_WORLD);
    Matrix result;
    result.lines = size / world_size;
    result.columns = size;
    result.tab = buffer;
    return result;
}

/**
 * Apply the matrix multiplication algorithm and reconstruct the matrix.
 * B lines and column are considered inverted.
 *
 * @param A
 * @param B
 * @param result
 * @param iteration The current segment of B being computed
 * @param rank The current process in which the computation occurs
 */
void multiplyMatrix(Matrix* A, Matrix* B, Matrix* result, int iteration, int rank) {
    int sum = 0;
    #pragma omp parallel for
    for (int i = 0; i < A->lines; ++i) {
        #pragma omp parallel for
        for (int j = 0; j < B->lines; ++j) {
            sum = 0;
            #pragma omp parallel for
            for (int k = 0; k < A->columns; ++k) {
                sum += get(A, i, k) * get(B, j, k);
            }
            //Make sure to replace the result in the right place in the matrix
            set(result, i, ((j + iteration * A->lines - (A->columns - (rank * A->lines)))
                            % A->columns + A->columns) % A->columns, sum);
        }
    }
}

void gather(Matrix* toSend, int size, Matrix* buffer) {
    MPI_Gather(toSend->tab, (int) (toSend->lines * size), MPI_INT, buffer->tab,
               (int) (toSend->lines * size), MPI_INT, 0, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);
    omp_set_num_threads(5);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    Matrix A;
    Matrix turned;
    long size; //number of element per line

    //Prepare the send -> read A and B, turn B
    if (world_rank == 0) {
        char* AFileName = argv[1];
        readMatrixFromFile(&A, AFileName);
        char *BFileName = argv[2];
        Matrix B;
        readMatrixFromFile(&B, BFileName);
        turnMatrix(&B, &turned);
        size = A.columns;
    }

    //Broadcast the column number of A, scatter a and B
    MPI_Bcast(&size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    Matrix receivedA = scatter(size, &A, world_size);
    Matrix receivedB = scatter(size, &turned, world_size);

    //Initialize and compute the first step of the result matrix
    Matrix afterCalcul;
    initMatrix(&afterCalcul, receivedA.lines, size);
    multiplyMatrix(&receivedA, &receivedB, &afterCalcul, 0, world_rank);

    //Send and receive parts of B from your neighbors,
    //then multiply it to your result until you received all of B.
    int prev;
    int next;
    for (int i = 1; i < world_size; ++i) {
        prev = (world_rank - 1) % world_size;
        if (world_rank == 0) {
            prev = world_size - 1;
        }
        MPI_Send(receivedB.tab, (int) (receivedB.lines * receivedB.columns), MPI_INT, prev, 0, MPI_COMM_WORLD);

        next = (world_rank + 1) % world_size;
        MPI_Recv(receivedB.tab, (int) (receivedB.lines * receivedB.columns),
                 MPI_INT, next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        multiplyMatrix(&receivedA, &receivedB, &afterCalcul, i, world_rank);
    }

    //Initialize, receive and print the final result.
    Matrix gathered;
    if (world_rank == 0) {
        initMatrix(&gathered, size, size);
    }

    gather(&afterCalcul, (int) size, &gathered);

    if (world_rank == 0) {
        printMatrix(&gathered);
    }

    MPI_Finalize();
}


