#ifndef PLOT2D_H
#define PLOT2D_H
//dim size of the vector
//name of the file e.g "test.html"
void plot2D(double* x, double* y, int dim, char* title, const char* name);
//dim_1 how many graphs
//dim size of the vector
//title , title of the graph e.g "Trace 0"
//name of the file e.g "test.html"
void plot2DMat(double* x, double** y, int dim_1, int dim_2, const char* title, const char* name);
//dim_1 how many traces (in the same graph)
//dim size of the vector
//title , title of the graph e.g "Trace 0"
//name of the file e.g "test.html"
void plot2DComparison(double* x, double** y, int dim_1, int dim_2, const char* title, const char* name);

void plot2DComparisonNames(double* x, double** y, int dim_1, int dim_2, const char* title, const char* name,  char** variables);

void plot2DExtraDots(double* x, double* y, double** dots, int dim_2, int dim_3, const char* title, const char* name);

void plotContour(double* x, double* y, double** z, int dim_1, int dim_2, const char* title, const char* name);

void plotPoints(double* x, double* y, double* z, int dim_1, const char* title, const char* name);

void plotPointsExtraDots(double* x, double* y, double* z, double** dots, int dim_1, int dim_2, const char* title, const char* name);

//NEW
void plot2DComparisonNamesNew(double** y, int dim_1, int dim_2, const char* title, const char* name, char** variables, const char* Yaxis);
#endif
