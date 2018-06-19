#ifndef OPENFILE_H
#define OPENFILE_H

typedef struct matDouble
{ 
    int dim_1;
    int dim_2;
    double** matrix;
}matDouble;

double** openMatrix(const char* name);
struct matDouble openMatrixDouble(const char* name);
struct matDouble openMatrixDouble2(const char* name, int dim_1, int dim_2);
double* openVector(const char* name, int dim_1);

int** openMatrixInt(const char* name);
#endif
