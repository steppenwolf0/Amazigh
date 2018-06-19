#include "../Solvers/solversAQ.h"
#include "../Amazigh/printFile.h"
#include "../Amazigh/openFile.h"
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/plot2D.h"
#include "../Amazigh/mathVector.h"
#include "../Amazigh/supportC.h"
#include "../Amazigh/showFile.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../Solvers/functions.h"
#include "../Gaussian/gaussian.h"
#ifdef _WIN32
#include "../Solvers/util/ToolsWIN.h"
#else
#include "../Solvers/util/Tools.h"
#endif
//#include<vld.h>
/**
* The main function:
* - interpret parameters on the command line
* - run the algorithm with the interpreted parameters
*/

#include "../Solvers/LBFGSB/driver1.h"
#include "../Solvers/CMAES/src/cmaes_interface.h"
#include "../Solvers/CMAES/src/boundary_transformation.h"

double** xtrainLBFGS;
double* ytrainLBFGS;

int Dimension = 2;	//2,5,10,20
int n_points = 5;
int evalsEA = 30000;	//30,000 75,000 150,000 300,000
int test = 0;	//Sphere,Ackley,Grienwank,Rosenbrock
int method = 2;	//Gradient,CMAES,GOMEA
int iterationsMax = 50; //100,150,300,500

double funcEval(double* xtrain)
{
	/*double* scaleX = newDoubleV(Dimension);
	for (int i = 0; i < Dimension; i++)
	{
	scaleX[i] = xtrain[i] * 200 - 100;
	}
	double tempY = cec14_eotp_problems(scaleX, Dimension, 2);
	free(scaleX);*/
	double tempY = 1e20;
	if (test == 0)
		tempY = sphereFunctionProblemEvaluation(xtrain, Dimension);
	if (test == 1)
		tempY = ackleyFunctionProblemEvaluation(xtrain, Dimension);
	if (test == 2)
		tempY = griewankFunctionProblemEvaluation(xtrain, Dimension);
	if (test == 3)
		tempY = rosenbrockFunctionProblemEvaluation(xtrain, Dimension);
	return tempY;
}


void generateInitialPoints(double** x, double*y, int n_points, int Dimension)
{
	initializeRandomNumberGenerator();


	int number = n_points;
	int dimensions = Dimension;
	// Generate slotted random points.
	double factor = 1.0 / number;

	double *points = newDoubleV(number * dimensions);
	double *start = points;

	for (int d = 0; d < dimensions; ++d)
	{
		// Generate points.
		for (int i = 0; i < number; ++i)
		{
			// Generate point.
			double value = (i + randomRealUniform01()) * factor;

			// Generate random index.
			int j = randomInt(i + 1);

			// Swap them.
			start[i] = start[j];
			start[j] = value;
		}

		// Advance.
		start += number;
	}

	// Assign points.
	for (int i = 0; i < number; ++i)
	{
		// Determine location.
		double *xDummy = newDoubleV(dimensions);

		for (int d = dimensions - 1; d >= 0; --d)
		{
			xDummy[d] = points[d * number + i];
			x[i][d] = xDummy[d];
		}

		// Determine value.
		double value = funcEval(x[i]);
		y[i] = value;


	}




}

double* randomVector(int size)
{
	initializeRandomNumberGenerator();
	double* x = newDoubleV(size);

	for (int d = 0; d < size; d++)
	{
		x[d] = randomRealUniform01();
	}

	return x;
}

double* gradientBasedNumeric(double** AinvKexp, double* ytrain, double* paramsKexp,
	double** xtrain, int n_points, int Dimension)
{
	double* wi = multVector(AinvKexp, ytrain, n_points, n_points);
	double* finalPoint = newDoubleV(Dimension);
	double bestFx = 1e20;
	for (int tries = 0; tries < 100; tries++)
	{
		double* randomPoint = randomVector(Dimension);

		showVector(randomPoint, Dimension);

		double* resultLBFGS = solveUCBBNumeric(AinvKexp, ytrain, n_points,
			Dimension, paramsKexp, xtrain,
			randomPoint, 1.96, 1);

		double** resultLBFGSTemp = newDouble(1, Dimension);
		for (int j = 0; j < Dimension; j++)
		{
			resultLBFGSTemp[0][j] = resultLBFGS[j];
		}

		double fx = acqFunctionUCB(AinvKexp, wi,
			n_points, Dimension, paramsKexp, xtrain, resultLBFGSTemp, 1.96, 1);
		printf("fx  \t %0.8f \n", fx);
		if (fx < bestFx)
		{
			bestFx = fx;
			copyVector(finalPoint, resultLBFGS, Dimension);
		}
		free(randomPoint);
		free(resultLBFGS);
		freeMatrix(resultLBFGSTemp, 1);
		printf("%d ", tries);
	}
	printf("\n");

	free(wi);

	return finalPoint;
}

double* gradientBased(double** AinvKexp, double* ytrain, double* paramsKexp,
	double** xtrain, int n_points, int Dimension)
{
	double* wi = multVector(AinvKexp, ytrain, n_points, n_points);
	double* finalPoint = newDoubleV(Dimension);
	double bestFx = 1e20;
	for (int tries = 0; tries < 1; tries++)
	{
		double* randomPoint = randomVector(Dimension);

		showVector(randomPoint, Dimension);

		double* resultLBFGS = solveUCBB(AinvKexp, ytrain, n_points,
			Dimension, paramsKexp, xtrain,
			randomPoint, 1.96, 1);

		double** resultLBFGSTemp = newDouble(1, Dimension);
		for (int j = 0; j < Dimension; j++)
		{
			resultLBFGSTemp[0][j] = resultLBFGS[j];
		}

		double fx = acqFunctionUCB(AinvKexp, wi,
			n_points, Dimension, paramsKexp, xtrain, resultLBFGSTemp, 1.96, 1);
		printf("fx  \t %0.8f \n", fx);
		if (fx < bestFx)
		{
			bestFx = fx;
			copyVector(finalPoint, resultLBFGS, Dimension);
		}
		free(randomPoint);
		free(resultLBFGS);
		freeMatrix(resultLBFGSTemp, 1);
		printf("%d ", tries);
	}
	printf("\n");

	free(wi);

	return finalPoint;
}

double* gomeaBased(double** AinvKexp, double** xtrain, double* paramsKexp, double* ytrain,
	int Dimension, int n_points, double kappa, int kernelType, int evals)
{
	Individual* test = NULL;
	//if (bestIterEI < 300)
	//{			
	//	printf("Expected Improvement \n");
	//	//test = RVGOMEA3A(AinvKexp, AinvKM52, xtrain, paramsKexp, paramsKM52, ytrainLBFGS, Dimension + 1, n_points, Dimension, bestValue, 0);
	//	test = RVGOMEAEI(AinvKexp, xtrainLBFGS, paramsKexp, ytrainLBFGS, Dimension + 1, n_points, Dimension, bestValue,1);
	//}
	////if (bestIterUCB < 5)
	////{
	//	freeIndividual(test);
	//	printf("UCB \n");
	//	//test = RVGOMEA1A(AinvKexp, AinvKM52, xtrain, paramsKexp, paramsKM52, ytrain, Dimension + 1, n_points, Dimension, 1.96);
	test = RVGOMEAUCB(AinvKexp, xtrain, paramsKexp, ytrain, Dimension + 1, n_points, Dimension, kappa, kernelType, evals);
	////}

	double* points = newDoubleV(Dimension);
	for (int i = 0; i < Dimension; i++)
	{
		points[i] = test->parameters[i];
	}
	freeIndividual(test);

	return points;
}

double fitfun(double const *x, int N, double** AinvKexp, double* ytrain, double* paramsKexp,
	double** xtrain) { /* function "cigtab" */

	int i = 0;
	int typeKernel = 1;
	double** testpoint = newDouble(1, N);
	for (i = 0; i < N; i++)
	{
		double val = x[i];
		testpoint[0][i] = val;
	}

	double* wi = multVector(AinvKexp, ytrain, n_points, n_points);

	double result = acqFunctionUCB(AinvKexp, wi, n_points,
		Dimension, paramsKexp, xtrain,
		testpoint, 1.96, typeKernel);

	free(wi);
	freeMatrix(testpoint, 1);
	return result;
}

double* cmaesBased(double** AinvKexp, double* ytrain, double* paramKexp, double** xtrain,
	int dimensionSize, int evals, int thread)
{
	char* tempName[20];
	sprintf(tempName, "cmaes_initials%d.par", thread);

	FILE *f = fopen(tempName, "w");
	//function number 1 2 5 6 8 23
	//	restarts 0
	//	N 10
	//	initialX 1:
	//0.5e0
	//	typicalX 1 :
	//	0.0
	//	initialStandardDeviations  1 :
	//	0.3
	//	stopMaxFunEvals   100000
	//	stopMaxIter       1e299
	//	stopTolFun 1e-12
	//	stopTolFunHist 1e-13
	//	stopTolX 1e-11
	//	stopTolUpXFactor 1e3
	//	seed 0
	//	maxTimeFractionForEigendecompostion 1
	fprintf(f, "function number 1 2 5 6 8 23\n");
	fprintf(f, "restarts 0\n");
	fprintf(f, "N %d\n", dimensionSize);
	fprintf(f, "initialX 1:\n");
	fprintf(f, "0.5e0\n");
	fprintf(f, "typicalX 1 :\n");
	fprintf(f, "0.0\n");
	fprintf(f, "initialStandardDeviations  1 :\n");
	fprintf(f, "0.3\n");
	fprintf(f, "stopMaxFunEvals   %d\n", evals);
	fprintf(f, "stopMaxIter       1e299\n");
	fprintf(f, "stopTolFun 1e-12\n");
	fprintf(f, "stopTolFunHist 1e-13\n");
	fprintf(f, "stopTolX 1e-11\n");
	fprintf(f, "stopTolUpXFactor 1e3\n");
	fprintf(f, "seed 0\n");
	fprintf(f, "maxTimeFractionForEigendecompostion 1\n");
	fclose(f);

	cmaes_t evo; /* an CMA-ES type struct or "object" */
	cmaes_boundary_transformation_t boundaries;
	double *arFunvals, *x_in_bounds, *const*pop;
	double lowerBounds[] = { 0 }; /* last number is recycled for all remaining coordinates */
	double upperBounds[] = { 1 };
	int nb_bounds = 1; /* numbers used from lower and upperBounds */
	unsigned long dimension;
	int i;

	/* initialize boundaries, be sure that initialSigma is smaller than upper minus lower bound */
	cmaes_boundary_transformation_init(&boundaries, lowerBounds, upperBounds, nb_bounds);
	/* Initialize everything into the struct evo, 0 means default */
	arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, tempName);
	dimension = (unsigned long)cmaes_Get(&evo, "dimension");

	printf("%s\n", cmaes_SayHello(&evo));

	x_in_bounds = cmaes_NewDouble(dimension); /* calloc another vector */
	cmaes_ReadSignals(&evo, "cmaes_signals.par");  /* write header and initial values */

												   /* Iterate until stop criterion holds */
	while (!cmaes_TestForTermination(&evo))
	{
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

											/* transform into bounds and evaluate the new search points */
		for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {
			cmaes_boundary_transformation(&boundaries, pop[i], x_in_bounds, dimension);
			/* this loop can be omitted if is_feasible is invariably true */
			//while (!is_feasible(x_in_bounds, dimension)) { /* is_feasible needs to be user-defined, in case, and can change/repair x */
			//	cmaes_ReSampleSingle(&evo, i);
			//	cmaes_boundary_transformation(&boundaries, pop[i], x_in_bounds, dimension);
			//}
			arFunvals[i] = fitfun(x_in_bounds, dimension,
				AinvKexp, ytrain, paramKexp, xtrain); /* evaluate */
		}

		/* update the search distribution used for cmaes_SampleDistribution() */
		cmaes_UpdateDistribution(&evo, arFunvals);  /* assumes that pop[i] has not been modified */

													/* read instructions for printing output or changing termination conditions */
		cmaes_ReadSignals(&evo, "cmaes_signals.par");
		fflush(stdout); /* useful in MinGW */
	}
	printf("Stop:\n%s\n", cmaes_TestForTermination(&evo)); /* print termination reason */
	cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

															/* get best estimator for the optimum, xmean */
	cmaes_boundary_transformation(&boundaries,
		(double const *)cmaes_GetPtr(&evo, "xmean"), /* "xbestever" might be used as well */
		x_in_bounds, dimension);

	/* do something with final solution x_in_bounds */

	/* ... */

	/* and finally release memory */
	cmaes_exit(&evo); /* release memory */
	cmaes_boundary_transformation_exit(&boundaries); /* release memory */
													 //free(x_in_bounds);

	return x_in_bounds;
}

double* moveBoundaries(double** xtrain, int n_points, int Dimension, double* point)
{

	if (n_points > 25)
	{
		double* vector = newDoubleV(Dimension);
		int initial = (n_points - 10);
		int points = 0;
		for (int i = initial; i < n_points; i++)
		{
			for (int j = 0; j < Dimension; j++)
			{
				vector[j] = vector[j] + xtrain[i][j];
			}
			points++;
		}
		for (int j = 0; j < Dimension; j++)
		{
			vector[j] = vector[j] / points;
			if (vector[j] <= 0.000001 || vector[j] >= 0.999999)
			{
				vector[j] = randomRealUniform01();
			}
			else
			{
				vector[j] = point[j];
			}
		}
		return vector;
	}
	else
	{
		return point;
	}

}

void code2DE(int runName)
{

	double* paramsKexp = newDoubleVInit(Dimension + 1, 1);
	double* paramsBestKexp = newDoubleVInit(Dimension + 1, 1);

	double** xtrain = newDouble(n_points, Dimension);
	double* ytrain = newDoubleV(n_points);

	generateInitialPoints(xtrain, ytrain, n_points, Dimension);

	double bestValue = 1e10;
	int bestIterEI = 0;
	int bestIterUCB = 10;
	int systemACQ = 0;

	int iterations = 0;




	char* tempName[20];
	sprintf(tempName, "results_%d.txt", runName);

	FILE *f = fopen(tempName, "w");
	for (iterations = 0; iterations < iterationsMax; iterations++)
	{

		double GAMMA = 0.5;
		double kappa = sqrt(2 * log(pow(iterations + 1, (Dimension / 2 + 2))* PI*PI / (3 * GAMMA)));


		int i, j;

		printf("LBFGS \n");
		xtrainLBFGS = newDouble(n_points, Dimension);
		ytrainLBFGS = newDoubleV(n_points);

		copyMatrix(xtrainLBFGS, xtrain, n_points, Dimension);
		copyVector(ytrainLBFGS, ytrain, n_points);

		double meanPop = averagePop(ytrainLBFGS, n_points);
		double averageFactor = 1.0 / n_points;
		double standardPop = standardDeviationPop(ytrain, n_points, meanPop, averageFactor);
		moveY(meanPop, standardPop, ytrainLBFGS, n_points);

		double* thetaKexp = solveKexpB(xtrainLBFGS, ytrainLBFGS, n_points, Dimension);


		copyVector(paramsKexp, thetaKexp, Dimension + 1);

		free(thetaKexp);

		//first system********************************************************************************************
		//second system********************************************************************************************
		double marginal_likelihoodKexp = marginalLikelihood(n_points, xtrain, Dimension, paramsKexp, ytrainLBFGS, 1);
		printf("ML new %.8f\t", marginal_likelihoodKexp);

		double marginal_likelihoodBestKexp = marginalLikelihood(n_points, xtrain, Dimension, paramsBestKexp, ytrainLBFGS, 1);
		printf("ML best %.8f\n", marginal_likelihoodBestKexp);

		double** AinvKexp;
		if (fabs(marginal_likelihoodBestKexp)<fabs(marginal_likelihoodKexp))
		{
			for (i = 0; i < Dimension + 1; i++)
			{
				paramsKexp[i] = paramsBestKexp[i];
			}
			AinvKexp = getAinv(n_points, xtrain, Dimension, paramsBestKexp, ytrainLBFGS, 1);
		}
		else
		{
			for (i = 0; i < Dimension + 1; i++)
			{
				paramsBestKexp[i] = paramsKexp[i];
			}
			AinvKexp = getAinv(n_points, xtrain, Dimension, paramsKexp, ytrainLBFGS, 1);
		}
		//second system********************************************************************************************
		double* finalPoint;
		if (method == 0)
			finalPoint = gradientBased(AinvKexp, ytrain, paramsKexp, xtrain, n_points, Dimension);
		if (method == 1)
			finalPoint = cmaesBased(AinvKexp, ytrain, paramsKexp, xtrain, Dimension, evalsEA, runName);
		if (method == 2)
		{
			finalPoint = gomeaBased(AinvKexp, xtrain, paramsKexp, ytrain, Dimension, n_points, 1.96, 1, evalsEA);

			/*double* finalPoint = gradientBasedNumeric(AinvKexp, ytrain, paramsKexp,xtrain, n_points, Dimension);*/

			finalPoint = moveBoundaries(xtrain, n_points, Dimension, finalPoint);
		}


		double** points = newDouble(1, Dimension);
		for (i = 0; i < Dimension; i++)
		{
			points[0][i] = finalPoint[i];
		}

		free(finalPoint);

		freeMatrix(xtrainLBFGS, n_points);

		double** xtrainTemp = newDouble(n_points, Dimension);
		copyMatrix(xtrainTemp, xtrain, n_points, Dimension);
		freeMatrix(xtrain, n_points);
		xtrain = newDouble(n_points + 1, Dimension);
		copyMatrix(xtrain, xtrainTemp, n_points, Dimension);
		for (j = 0; j < Dimension; j++)
		{
			xtrain[n_points][j] = points[0][j];
		}
		freeMatrix(xtrainTemp, n_points);

		double* ytrainTemp = newDoubleV(n_points);
		copyVector(ytrainTemp, ytrain, n_points);
		free(ytrain);
		ytrain = newDoubleV(n_points + 1);
		copyVector(ytrain, ytrainTemp, n_points);
		ytrain[n_points] = funcEval(xtrain[n_points]);
		free(ytrainTemp);


		freeMatrix(points, 1);

		freeMatrix(AinvKexp, n_points);



		FILE *file = fopen(tempName, "a");
		for (j = 0; j < Dimension; j++)
		{
			fprintf(file, "%.6f\t", xtrain[n_points][j]);
		}
		fprintf(file, "%.11f\t", ytrain[n_points]);

		if (ytrain[n_points]<bestValue)
		{
			bestValue = ytrain[n_points];
			if (systemACQ == 0)
			{
				bestIterEI = 0;
			}
			if (systemACQ == 1)
			{
				bestIterUCB = 0;
			}
		}
		else
		{
			if (systemACQ == 0)
			{
				bestIterEI++;
				if (bestIterEI == 300)
				{
					bestIterUCB = 0;
					systemACQ = 1;
				}
			}
			if (systemACQ == 1)
			{
				bestIterUCB++;
				if (bestIterUCB == 5)
				{
					systemACQ = 0;
					bestIterEI = 0;
				}
			}
		}

		fprintf(file, "%.11f\t", bestValue);
		fprintf(file, "%.11f\n", marginal_likelihoodKexp);

		fclose(file);


		printf("bestIterEI %d\t", bestIterEI);
		printf("bestIterUCB %d\t", bestIterUCB);

		printf("%.11f\t", ytrain[n_points]);
		printf("%.11f\t", kappa);

		printf("%d\t", iterations);
		printf("%.11f\t\n", bestValue);

		n_points++;
	}


}