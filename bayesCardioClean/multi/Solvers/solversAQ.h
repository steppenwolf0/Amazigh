#ifndef SOLVERSAQ_H
#define SOLVERSAQ_H
#include "problemDef.h"

Individual* RVGOMEAUCB(double** Ainv, double** xtrain, double* param, double* ytrain,
	int n_param, int n_x, int Dimension, double kappa, int kernelType, int evals);

Individual* RVGOMEAEI(double** Ainv, double** xtrain, double* param, double* ytrain, int n_param, int n_x, int Dimension,
	double bestSolution, int kernelType);


#endif
