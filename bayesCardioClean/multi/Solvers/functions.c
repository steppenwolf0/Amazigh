#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../Amazigh/supportC.h"
#define PI 3.141592653589793238462643 


//double griewankFunctionProblemEvaluation(double *parameters, int Dimension)
//{
//	int    i;
//	double yi, sum, prod, result;
//
//	sum = 0;
//	prod = 1.0;
//	for (i = 0; i < Dimension; i++)
//	{
//		double val = (100 - (-100))*parameters[i] + (-100);
//		yi = val;
//		sum += yi*yi;
//
//		yi = (val) / sqrt((double)(i + 1));
//		prod *= cos(yi);
//	}
//
//	result = sum / 4000.0 - prod + 1.0;
//
//	return result;
//}
//
//double sphereFunctionProblemEvaluation(double* xtrain, int Dimension)
//{
//	int j;
//	double tempY = 0;
//	for (j = 0; j < Dimension; j++)
//	{
//		double val = xtrain[j];
//		val = (20 - (-20))*val + (-20);
//		tempY += val*val;
//	}
//	return tempY;
//}
//double michalewiczFunctionProblemEvaluation(double *xtrain, int Dimension)
//{
//	int    i;
//	double result;
//
//	result = 0.0;
//	for (i = 0; i < Dimension; i++)
//	{
//		double val = xtrain[i];
//		val = (PI - (0))*val + (0);
//		result += -sin(val)*pow(sin(((i + 1)*val * val) / PI), 20.0);
//	}
//
//
//	return result;
//}

const double point[100] =  {
	52.34, -85.34, -50.21, 56.2, 89.55, 10.07, -63.64, -77.9, 66.73, 64.45,
		-98.51, -51.52, -41.99, -17.25, 48.74, -44.86, -87, 48.13, 30.18, -53.53, 30.49, 92.69, 10.46,
		52.01, 45.3, -87.43, 0.46, -71.61, 25.69, -49.87, -74.86, -8.04, -9.43, 3.78, 82.61, -57.76,
		11.6, -71.07, 99.42, -11.39, 73.25, 84.24, 22.3, 49.72, 4.86, -53.32, 98.37, -96.25, 96.65,
		37.24, 72.7, 95.84, 26.33, 62.58, 16.86, -86.87, -90.32, 28.12, 20.44, -6.52, 98.74, -60.42,
		4.48, 64.53, 7.89, 21.54, -5.37, 33.11, -10.98, -33.08, 34.87, 0.87, -0.93, 22.94, 55.97, 2,
		43.27, -94.82, -52.26, -65.53, 30.72, -54.08, -78.03, 42.36, -87.16, -85.26, -7.77, 51.21,
		38.54, 27.73, 3.31, 22.61, 87.26, 70.41, -57.51, -19.36, 33.78, -92.23, -35.92, -53.91
};

double griewankFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	
		double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}

	double yi, sum, prod, result;

	sum = 0;
	prod = 1.0;
	for (i = 0; i < Dimension; i++)
	{
		yi = parameters[i] - point[i];
		sum += yi*yi;

		yi = (parameters[i] - point[i]) / sqrt((double)(i + 1));
		prod *= cos(yi);
	}

	result = sum / 4000.0 - prod + 1.0;

	free(parameters);
	return result;
}

double sphereFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	
	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i]-point[i];
	}
	double result;

	result = 0.0;
	for (i = 0; i < Dimension; i++)
		result += parameters[i] * parameters[i];

	free(parameters);
	return result;

}

double rosenbrockFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	
	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i] - point[i];
	}
	double result;

	result = 0.0;
	for (i = 0; i < Dimension - 1; i++)
		result += 100 * (parameters[i + 1] - parameters[i] * parameters[i])*(parameters[i + 1] - parameters[i] * parameters[i]) + (1.0 - parameters[i])*(1.0 - parameters[i]);

	free(parameters);
	return result;
}

double ackleyFunctionProblemEvaluation(double *xtrain, int Dimension) {
	

	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i] - point[i];
	}
	double sum1, sum2;
	double result;



	sum1 = 0.0L;
	sum2 = 0.0L;
	for (i = 0; i < Dimension; i++) {
		sum1 += (parameters[i] * parameters[i]);
		sum2 += cos(2.0L * PI * parameters[i]);
	}
	double e = 2.7182818284590452353602874713526625L;
	result = -20.0L * exp(-0.2L * sqrt(sum1 / Dimension)) - exp(sum2 / Dimension) + 20.0L + e;

	free(parameters);
	return result;
}

double michalewiczFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i] - point[i];
	}
	double result;

	result = 0.0;
	for (i = 0; i < Dimension; i++)
		result += -sin(parameters[i])*pow(sin(((i + 1)*parameters[i] * parameters[i]) / PI), 20.0);

	free(parameters);
	return result;
}
