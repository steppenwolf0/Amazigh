#ifndef SUPPORT_LBFGSB_H
#define SUPPORT_LBFGSB_H
double* solveKexpB(double** xtrainLBFGSA, double* ytrainLBFGSA, int n_pointsA, int DimensionA);
double* solveUCBB(double** Ainv, double* ytrainLBFGSA, int n_pointsA,
	int DimensionA, double* params, double** xtrainLBFGSA,
	double* randomPoint, double kappa, int type);
double* solveUCBBNumeric(double** Ainv, double* ytrainLBFGSA, int n_pointsA,
	int DimensionA, double* params, double** xtrainLBFGSA,
	double* randomPoint, double kappa, int type);
#endif