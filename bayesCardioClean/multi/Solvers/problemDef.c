/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "../Amazigh/mathVector.h"
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/openFile.h"
#include "../Amazigh/showFile.h"
#include "../Amazigh/printFile.h"
#include "../Amazigh/supportC.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "../Solvers/problemDef.h"
#include "../Gaussian/gaussian.h"
#include "./solversAQ.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798

double** xtrain;

double* ytrain;
int n_x;
int n_param;
double kappa;
int type;

int Dimension;
double bestSolutionA;

//double** AinvKexp;
//double** AinvKM52;
//double* paramKexp;
//double* paramKM52;
//double* wiKexp;
//double** transWiKexp;
//double* wiKM52;
//double** transWiKM52;

double** Ainv;
double* param;
double* wi;
double** transWi;

void initializeproblemUCB(double** _Ainv, double** _xtrain, double* _param, double* _ytrain, int _n_param, int _n_x,
	int _Dimension, double _Kappa, int _type)
{
	type = _type;
	Dimension = _Dimension;
	n_x = _n_x;
	n_param = _n_param;
	kappa = _Kappa;
	number_of_parameters = Dimension;

	Ainv = newDouble(n_x, n_x);
	copyMatrix(Ainv, _Ainv, n_x, n_x);

	xtrain = newDouble(n_x, Dimension);
	copyMatrix(xtrain, _xtrain, n_x, Dimension);

	ytrain = newDoubleV(n_x);
	copyVector(ytrain, _ytrain, n_x);

	param = newDoubleV(n_param);
	copyVector(param, _param, n_param);

	wi = multVector(Ainv, ytrain, n_x, n_x);

	transWi = transVector(wi, n_x);
	
	}

void initializeproblemEI(double** _Ainv, double** _xtrain, double* _param, double* _ytrain, int _n_param, int _n_x,
	int _Dimension, double _bestSolution)
{
	
	Dimension = _Dimension;
	n_x = _n_x;
	n_param = _n_param;
	bestSolutionA = _bestSolution;
	number_of_parameters = Dimension;

	Ainv = newDouble(n_x, n_x);
	copyMatrix(Ainv, _Ainv, n_x, n_x);

	xtrain = newDouble(n_x, Dimension);
	copyMatrix(xtrain, _xtrain, n_x, Dimension);

	ytrain = newDoubleV(n_x);
	copyVector(ytrain, _ytrain, n_x);

	param = newDoubleV(n_param);
	copyVector(param, _param, n_param);

	wi = multVector(Ainv, ytrain, n_x, n_x);

	transWi = transVector(wi, n_x);

}

//void initializeproblem2A(double** _AinvKexp, double** _AinvKM52, double** _xtrain, double* _paramKexp, 
//	double* _paramKM52, double* _ytrain, int _n_param, int _n_x, int _Dimension, double _bestSolution, double power)
//{
//	power = power;
//	Dimension = _Dimension;
//	n_x = _n_x;
//	n_param = _n_param;
//	bestSolutionA = _bestSolution;
//	number_of_parameters = Dimension;
//
//	AinvKexp = newDouble(n_x, n_x);
//	copyMatrix(AinvKexp, _AinvKexp, n_x, n_x);
//	AinvKM52 = newDouble(n_x, n_x);
//	copyMatrix(AinvKM52, _AinvKM52, n_x, n_x);
//
//	xtrain = newDouble(n_x, Dimension);
//	copyMatrix(xtrain, _xtrain, n_x, Dimension);
//	
//	ytrain = newDoubleV(n_x);
//	copyVector(ytrain, _ytrain, n_x);
//	/*double meanPop = averagePop(ytrain, n_x);
//	double averageFactor = 1.0 / n_x;
//	double standardPop = standardDeviationPop(ytrain, n_x, meanPop, averageFactor);
//	moveY(meanPop, standardPop, ytrain, n_x);*/
//		
//	paramKexp = newDoubleV(n_param);
//	copyVector(paramKexp, _paramKexp, n_param);
//
//	paramKM52 = newDoubleV(n_param);
//	copyVector(paramKM52, _paramKM52, n_param);
//
//	wiKexp = multVector(AinvKexp, ytrain, n_x, n_x);
//	wiKM52 = multVector(AinvKM52, ytrain, n_x, n_x);
//
//	transWiKexp = transVector(wiKexp, n_x);
//	transWiKM52 = transVector(wiKM52, n_x);
//}
//
//void initializeproblem0A(double** _AinvKexp, double** _AinvKM52, double** _xtrain, double* _paramKexp, double* _paramKM52, 
//	double* _ytrain, int _n_param, int _n_x, int _Dimension, double _Kappa)
//{
//
//	Dimension = _Dimension;
//	n_x = _n_x;
//	n_param = _n_param;
//	kappa = _Kappa;
//	number_of_parameters = Dimension;
//
//	AinvKexp = newDouble(n_x, n_x);
//	copyMatrix(AinvKexp, _AinvKexp, n_x, n_x);
//	AinvKM52 = newDouble(n_x, n_x);
//	copyMatrix(AinvKM52, _AinvKM52, n_x, n_x);
//
//	xtrain = newDouble(n_x, Dimension);
//	copyMatrix(xtrain, _xtrain, n_x, Dimension);
//
//	ytrain = newDoubleV(n_x);
//	copyVector(ytrain, _ytrain, n_x);
//
//	paramKexp = newDoubleV(n_param);
//	copyVector(paramKexp, _paramKexp, n_param);
//
//	paramKM52 = newDoubleV(n_param);
//	copyVector(paramKM52, _paramKM52, n_param);
//
//	wiKexp = multVector(AinvKexp, ytrain, n_x, n_x);
//	wiKM52 = multVector(AinvKM52, ytrain, n_x, n_x);
//
//	transWiKexp = transVector(wiKexp, n_x);
//	transWiKM52 = transVector(wiKM52, n_x);
//
//}

void FunctionProblemEvaluation( Individual *ind, int number_of_parameters  )
{
	//UCB
	if (problem_Solving == 0)
	{
		int i = 0;
		double** testpoint = newDouble(1, number_of_parameters);
		for (i = 0; i < number_of_parameters; i++)
		{
			double val = ind->parameters[i];
			testpoint[0][i] = val;
		}
			
		double result = acqFunctionUCB(Ainv, wi, n_x,
			Dimension, param, xtrain,
			testpoint, kappa, type);


		ind->objective_value = result;
		ind->constrained_value = 0;

		freeMatrix(testpoint, 1);
	}
	//EI
	if (problem_Solving == 2)
	{
		int i = 0;
		double** testpoint = newDouble(1, number_of_parameters);
		for (i = 0; i < number_of_parameters; i++)
		{
			double val = ind->parameters[i];
			testpoint[0][i] = val;
		}

		double result= acqFunctionEI(Ainv,  wi,  n_x,
			Dimension, param,  xtrain, testpoint, bestSolutionA,  type);
		

		ind->objective_value = result;

		ind->constrained_value = 0;

		freeMatrix(testpoint, 1);
	}
	////EI 2 Kernels
	//if (problem_Solving == 20)
	//{
	//	int i = 0;
	//	double** testpoint = newDouble(1, number_of_parameters);
	//	for (i = 0; i < number_of_parameters; i++)
	//	{
	//		double val = ind->parameters[i];
	//		testpoint[0][i] = val;
	//	}


	//	double result = acqFunctionEI2Kernels(AinvKexp, AinvKM52, ytrain, n_x, number_of_parameters
	//		, paramKexp, paramKM52, xtrain, testpoint, bestSolutionA);
	//
	//	
	//	ind->objective_value = result;

	//	ind->constrained_value = 0;

	//	freeMatrix(testpoint, 1);
	//}
	////UCB 2 Kernels
	//if (problem_Solving == 100)
	//{
	//	int i = 0;
	//	double** testpoint = newDouble(1, number_of_parameters);
	//	for (i = 0; i < number_of_parameters; i++)
	//	{
	//		double val = ind->parameters[i];
	//		testpoint[0][i] = val;
	//	}

	//	double result = acqFunctionUCB2Kernels(AinvKexp, AinvKM52, ytrain, n_x, number_of_parameters,
	//		paramKexp, paramKM52, xtrain, testpoint, kappa);

	//	ind->objective_value = result;
	//	ind->constrained_value = 0;

	//	freeMatrix(testpoint, 1);
	//}
}

void FunctionPartialProblemEvaluation( Individual* ind, int number_of_touched_parameters, int *touched_parameters_indices, 
double *touched_parameters, double *parameters_before, double objective_value_before, 
double constraint_value_before )
{
	/*ind->objective_value = 0;
	ind->constrained_value = 0;*/
}

double functionProblemLowerRangeBound( int dimension )
{
	return(0);
}

double functionProblemUpperRangeBound( int dimension )
{

	return(1);
}

void copyIndividual(Individual *ind2, Individual *ind1, int number_of_parameters)
{
	int j;
	for (j = 0; j < number_of_parameters; j++)
		ind2->parameters[j] = ind1->parameters[j];

	ind2->objective_value = ind1->objective_value;
	ind2->constrained_value = ind1->constrained_value;
}

void initializeIndividuals(int maximum_number_of_populations)
{
	populations = (Individual ***)malloc(maximum_number_of_populations*sizeof(Individual **));
	selections = (Individual ***)malloc(maximum_number_of_populations*sizeof(Individual **));

	
}

void initializePopulations(int population_index, int number_of_parameters, int* population_sizes, int* selection_sizes)
{
	int j;
	
	populations[population_index] = (Individual **)malloc(population_sizes[population_index] * sizeof(Individual*));
	for (j = 0; j < population_sizes[population_index]; j++)
		populations[population_index][j] = initializeIndividual(number_of_parameters);

	selections[population_index] = (Individual **)malloc(selection_sizes[population_index] * sizeof(Individual*));
	for (j = 0; j < selection_sizes[population_index]; j++)
		//selections[population_index][j].parameters = (double *)Malloc(number_of_parameters*sizeof(double));
		selections[population_index][j] = initializeIndividual(number_of_parameters);
}

Individual* initializeIndividual(int number_of_parameters)
{
	Individual* ind = (Individual *)malloc(sizeof(Individual));
	ind->parameters = (double *)malloc(number_of_parameters*sizeof(double));

	int i;
	for (i = 0; i < number_of_parameters; i++)
	{
		ind->parameters[i] = 0.0;
	}
	ind->constrained_value = 0.0;
	ind->objective_value = 0.0;



	return ind;
}

void freeIndividual(Individual* ind)
{
	if (ind!=NULL)
	{
		if(ind->parameters!=NULL)
		free(ind->parameters);
		free(ind);
	}
	
	
}

void closeProblem()
{

	if (Ainv != NULL)
		freeMatrix(Ainv, n_x);

	if (xtrain != NULL)
		freeMatrix(xtrain, n_x);

	if (param != NULL)
		free(param);

	if (ytrain != NULL)
		free(ytrain);

	if (transWi != NULL)
		freeMatrix(transWi, 1);

	if (wi != NULL)
		free(wi);

}

//void closeProblemA()
//{
//
//
//	if (AinvKexp != NULL)
//		freeMatrix(AinvKexp, n_x);
//	
//	if (AinvKM52 != NULL)
//		freeMatrix(AinvKM52, n_x);
//	
//
//
//	if (xtrain != NULL)
//		freeMatrix(xtrain, n_x);
//	
//
//
//	if (paramKexp != NULL)
//		free(paramKexp);
//
//	if (paramKM52 != NULL)
//		free(paramKM52);
//
//	if (ytrain != NULL)
//		free(ytrain);
//
//	if (transWiKexp != NULL)
//		freeMatrix(transWiKexp, 1);
//	
//	if (transWiKM52 != NULL)
//		freeMatrix(transWiKM52, 1);
//	
//	
//	if (wiKexp != NULL)
//		free(wiKexp);
//
//	if(wiKM52 != NULL )
//		free(wiKM52);
//}
