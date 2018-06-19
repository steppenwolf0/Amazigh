#ifndef MATH_MATRIX_H
#define MATH_MATRIX_H

double **idMatrix (int n);
double** multMatrix(double** Matrix1, double** Matrix2, int dim_1M1,int dim_2M1,int dim_1M2,int dim_2M2);
void swapRows(double** A, int r1, int r2,int m);
double** inverseMatrix(double** matrix, int dim);
double** addMatrix(double** Matrix1, double** Matrix2, int dim_1,int dim_2 );
double** multConstant(double** Matrix1, double constant, int dim_1,int dim_2);
double** subMatrix(double** Matrix1, double** Matrix2, int dim_1,int dim_2 );
double** transposeMatrix(double** Matrix, int dim_1,int dim_2);
double detMatrix( double** matrix, int n);
double* multVector(double** Matrix1, double* Vector, int dim_1,int dim_2);
double** cutMatrix(double** Matrix, int x0, int y0, int xn, int yn, int dim_1,int dim_2 );
double* getVector(double** Matrix,int index, int dim);
double frobeniusNorm(double** Matrix, int dim_1, int dim_2);
double matrixNorm(double** Matrix, int dim_1, int dim_2);
double traceMatrix(double** Matrix, int size);
double logDet(double **a, int n);
double** Cholesky(double** matrix, int size);
double det4(double** m);
#endif
