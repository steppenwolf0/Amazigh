#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H

double* idVector(int n, double constant);
double dotVector(double* vector1, double* vector2, int dim);
double* multConstantV(double* vector, double constant, int dim);
double* cross3(double* a, double* b);
double norm(double* vector, int dim);
double* addVector(double* Vector1, double* Vector2, int dim);
double* subVector(double* Vector1, double* Vector2, int dim);
double** transVector(double* vector, int dim);

double maxVector(double* Vector, int dim);
#endif
