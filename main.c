#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <assert.h>

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

    // Print off a hello world message
    /*printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);*/

    int number;
    if(world_rank == 0) {
        for(int i = 1; i < world_size; i++) {
            printf("send %d\n", world_rank);
            number = i * 10;
            MPI_Send(&number, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
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

long countLinesOfFile(char* fileName) {
    int ch = 0;
    long lines = 0;
    FILE *file;
    file = fopen(fileName, "r");
    if (file) {
        while (!feof(file)) {
            ch = fgetc(file);
            if (ch == '\n') {
                lines++;
            }
        }
        fclose(file);
    }
    return lines;
}

long countColumnsOfFile(char* fileName) {
    long columns = 0;
    int ch = 0;
    int elem = 0;
    FILE *file;
    file = fopen(fileName, "r");
    if (file) {
        while (fscanf(file, "%d", &elem) > 0) {
            columns++;
            ch = fgetc(file);
            if (ch == '\n') {
                fclose(file);
                return columns;
            }
        }
        fclose(file);
    }
    return columns;
}

void readMatrixFromFile(Matrix* matrix, char* fileName) {
    long lines = countLinesOfFile(fileName);
    long columns = countColumnsOfFile(fileName);
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
}

int main(int argc, char** argv) {
    //mpi();
    char* AFileName = argv[1];

    Matrix A;
    readMatrixFromFile(&A, AFileName);
    //printMatrix(&A);

    char* BFileName = argv[2];

    Matrix B;
    readMatrixFromFile(&B, BFileName);
    //printMatrix(&B);
}


