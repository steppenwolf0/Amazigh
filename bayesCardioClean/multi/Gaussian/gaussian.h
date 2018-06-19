#ifndef GAUSSIAN_H
#define GAUSSIAN_H

double* randomNormalDistribution(int size);
double** KernelMV(double** xKernel, double** x2Kernel, int size_x, int size_x2, int Dimension, double* k_param, int type);
void shuffleIntV(int *vectorS, int n);
double cdf(double x);
double pdf(double x);
//
double marginalLikelihood(int n_points, double** xtrain, int Dimension, double* params, double* ytrain, int type);
double** getAinv(int n_points, double** xtrain, int Dimension, double* params, double* ytrain, int type);
double getMean(double* wi, int n_points, int Dimension, double* params, double** xtrain, double** points, int type);
double getVar(double** Ainv, int n_points, int Dimension, double* params, double** xtrain, double** points, int type);
//
double getCumulativeNormalValue(double x, double mean, double stdDev);
double averagePop(double* y, int size);
double standardDeviationPop(double* y, int size, double average, double averageFactor);
void moveY(double average, double standardDeviation, double* y, int size);
//
double acqFunctionEI2Kernels(double** AinvKexp, double** AinvKM52, double* wiKexp, double* wiKm52, int n_x,
	int number_of_parameters, double* paramKexp, double*paramKM52, double** xtrain,
	double** testpoint, double bestSolutionA);
double acqFunctionUCB2Kernels(double** AinvKexp, double** AinvKM52, double* wiKexp, double* wiKm52, int n_x,
	int number_of_parameters, double* paramKexp, double*paramKM52, double** xtrain,
	double** testpoint, double kappa);
//
double acqFunctionEI(double** Ainv, double* wi, int n_x, int dimension,
	double* params, double** xtrain, double** testpoint, double bestSolutionA, int type);
double acqFunctionUCB(double** Ainv, double* wi, int n_x,
	int dimension, double* params, double** xtrain,
	double** testpoint, double kappa, int type);
//
double getDerUCB(double** Ainv, double* ytrain, int n_points, int Dimension, double* params,
	double** xtrain, double** points, int type, int variableIndex, double kappa);
#endif
