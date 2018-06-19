#ifndef SOLVERS_H
#define SOLVERS_H

double* multVectorSparse(double* Avalues, int** Apositions, int Asize, double* x, int sizeVector);


double* PBCGSTABC(double** A, double* b, int size);

double* PBCGSTABCSparse(double* Avalues, int** Apositions, int Asize, double* b, int size);

void eulerODE(double* RATES, double* STATES, int sizeSTATES, double dt);

double* SolveSparsePCG(int** IJ,
	double* V, int length, double* B, double* x0, int size, double ptol);
#endif