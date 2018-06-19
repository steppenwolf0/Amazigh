#ifndef FUNC_H
#define FUNC_H

double sphereFunctionProblemEvaluation(double *xtrain, int Dimension);
double griewankFunctionProblemEvaluation(double *xtrain, int Dimension);
double rosenbrockFunctionProblemEvaluation(double *xtrain, int Dimension);
double ackleyFunctionProblemEvaluation(double *xtrain, int Dimension);
double michalewiczFunctionProblemEvaluation(double *xtrain, int Dimension);

#endif