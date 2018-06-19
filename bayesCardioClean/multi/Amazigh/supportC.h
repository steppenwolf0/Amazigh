#ifndef SUPPORT_H
#define SUPPORT_H

double** newDouble(int dim_1M1, int dim_2M2);
double* newDoubleV(int dim);
int* newIntV(int dim);
void freeMatrix(double** matrix, int dim);
double* newDoubleVInit(int dim, double constant);
void copyMatrix(double** matrix1, double** matrix2, int size1, int size2);
void copyVector(double* vector1, double* vector2, int size);
int** newInt(int dim_1, int dim_2);
int** newIntValue(int dim_1, int dim_2, int value);
void freeMatrixInt(int** matrix, int dim);
void copyMatrixInt(int** matrix1, int** matrix2, int size1, int size2);
//USAGE
//int size = 0;
//int** List;
//size = listInt2DExists(&List, 0, 5, size);
//printf("size %d \n", size);
//showMatrixInt(List, size, 2);
int listInt2D(int*** List, int index0, int index1, int size);
int listInt2DExists(int*** List, int index0, int index1, int size);
void copyVectorConst(double* vector1, double* vector2, int size, double constant);
#endif
