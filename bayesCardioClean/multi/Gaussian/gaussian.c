#ifdef _WIN32
	#include <time.h>
	#include "../Solvers/util/ToolsWIN.h"
#else
	#include "../Solvers/util/Tools.h"
	#if HAVE_SYS_TIME_H
		#include <sys/time.h>
	#endif
#endif

typedef int bool;
#define true 1
#define false 0

#include <math.h>
#include "../Amazigh/mathVector.h"
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/supportC.h"

#include "gaussian.h"
#define Pi 3.141592653589793238462643 
#define M_2_SQRTPI 1.12837916709551257390   // 2/sqrt(pi)
#define M_SQRT2    1.41421356237309504880   // sqrt(2)


double* randomNormalDistribution(int size)
{
	initializeRandomNumberGenerator();
	double* vector = newDoubleV(size);

	int i;
	for (i = 0; i < size; i++)
	{
		vector[i] = random1DNormalUnit();
	}

	return vector;
}

double exponentialKernelExp(double** xKernel, double** yKernel, int indexI, int indexJ, int Dimension, double* k_param)
{
	double* xI = newDoubleV(Dimension);
	double* xJ = newDoubleV(Dimension);
	int i;

	for (i = 0; i < Dimension; i++)
	{
		xI[i] = xKernel[indexI][i];
		xJ[i] = yKernel[indexJ][i];
	}

	double* temp = subVector(xI, xJ, Dimension);

	double** id = idMatrix(Dimension);
	for (i = 0; i < Dimension; i++)
	{
		id[i][i] = exp(-(k_param[i]* k_param[i]));
	}

	double* tempVector = multVector(id, temp, Dimension, Dimension);

	double norm = dotVector(temp, tempVector, Dimension);

	double r = sqrt(norm);

	free(temp);
	free(xI);
	free(xJ);
	free(tempVector);
	for (i = 0; i < Dimension; i++)
	{
		free(id[i]);
	}
	free(id);
	

	double exponential_kernel = exp(-(k_param[Dimension]* k_param[Dimension])) * exp(-0.5 * norm);

	return exponential_kernel;
}

double matern52KernelExp(double** xKernel, double** yKernel, int indexI, int indexJ, int Dimension, double* k_param)
{
	double* xI = newDoubleV(Dimension);
	double* xJ = newDoubleV(Dimension);
	int i;

	for (i = 0; i < Dimension; i++)
	{
		xI[i] = xKernel[indexI][i];
		xJ[i] = yKernel[indexJ][i];
	}

	double* temp = subVector(xI, xJ, Dimension);

	double** id = idMatrix(Dimension);
	for (i = 0; i < Dimension; i++)
	{
		id[i][i] = exp(-(k_param[i] * k_param[i]));
	}

	double* tempVector = multVector(id, temp, Dimension, Dimension);

	double norm = dotVector(temp, tempVector, Dimension);

	double r = sqrt(norm);

	free(temp);
	free(xI);
	free(xJ);
	free(tempVector);
	for (i = 0; i < Dimension; i++)
	{
		free(id[i]);
	}
	free(id);

	double matern52_kernel = exp(-(k_param[Dimension] * k_param[Dimension]))*(1 + sqrt(5 * norm) + 5 / 3 * norm)*exp(-sqrt(5 * norm));

	return matern52_kernel;
}


double exponentialKernel(double** xKernel, double** yKernel, int indexI, int indexJ, int Dimension, double* k_param)
{
	double* xI = newDoubleV(Dimension);
	double* xJ = newDoubleV(Dimension);
	int i;

	for (i = 0; i < Dimension; i++)
	{
		xI[i] = xKernel[indexI][i];
		xJ[i] = yKernel[indexJ][i];
	}

	double* temp = subVector(xI, xJ, Dimension);

	double** id = idMatrix(Dimension);
	for (i = 0; i < Dimension; i++)
	{
		id[i][i] = (1/(k_param[i] * k_param[i]));
	}

	double* tempVector = multVector(id, temp, Dimension, Dimension);

	double norm = dotVector(temp, tempVector, Dimension);

	double r = sqrt(norm);

	free(temp);
	free(xI);
	free(xJ);
	free(tempVector);
	for (i = 0; i < Dimension; i++)
	{
		free(id[i]);
	}
	free(id);


	double exponential_kernel = ((k_param[Dimension] * k_param[Dimension])) * exp(-0.5 * norm);

	return exponential_kernel;
}

double matern52Kernel(double** xKernel, double** yKernel, int indexI, int indexJ, int Dimension, double* k_param)
{

	double* xI = newDoubleV(Dimension);
	double* xJ = newDoubleV(Dimension);
	int i;

	for (i = 0; i < Dimension; i++)
	{
		xI[i] = xKernel[indexI][i];
		xJ[i] = yKernel[indexJ][i];
	}

	double* temp = subVector(xI, xJ, Dimension);

	double** id = idMatrix(Dimension);
	for (i = 0; i < Dimension; i++)
	{
		id[i][i] = (1 / (k_param[i] * k_param[i]));
	}

	double* tempVector = multVector(id, temp, Dimension, Dimension);

	double norm = dotVector(temp, tempVector, Dimension);

	double r = sqrt(norm);

	free(temp);
	free(xI);
	free(xJ);
	free(tempVector);
	for (i = 0; i < Dimension; i++)
	{
		free(id[i]);
	}
	free(id);

	double matern52_kernel = ((k_param[Dimension] * k_param[Dimension]))*(1 + sqrt(5 * norm) + 5 / 3 * norm)*exp(-sqrt(5 * norm));

	return matern52_kernel;
}


double** KernelMV(double** xKernel, double** x2Kernel, int size_x, int size_x2, int Dimension, double* k_param, 
	int type)
{


	double** cov = newDouble(size_x, size_x2);
	int i, j;

	for (i = 0; i< size_x; i++)
	{
		for (j = 0; j< size_x2; j++)
		{
		
			if (type == 0)
			{
				cov[i][j] = matern52KernelExp(xKernel, x2Kernel, i, j, Dimension, k_param);
			}
			if (type == 1)
			{
				cov[i][j] = exponentialKernelExp(xKernel, x2Kernel, i, j, Dimension, k_param);
			}

			if (type == 2)
			{
				cov[i][j] = matern52Kernel(xKernel, x2Kernel, i, j, Dimension, k_param);
			}
			if (type == 3)
			{
				cov[i][j] = exponentialKernel(xKernel, x2Kernel, i, j, Dimension, k_param);
			}

		}
	}


	return cov;
}

void shuffleIntV(int *vectorS, int n)
{
	initializeRandomNumberGenerator();
	if (n > 1)
	{
		int i;
		for (i = 0; i < n - 1; i++)
		{
			int j = i + rand() / (RAND_MAX / (n - i) + 1);

			int t = vectorS[j];
			vectorS[j] = vectorS[i];
			if (t < 0)
				t = 1;
			vectorS[i] = t;
		}
	}
}

double cdf(double x)
{
	 double RT2PI = sqrt(4.0*acos(0.0));

	static const double SPLIT = 7.07106781186547;

	static const double N0 = 220.206867912376;
	static const double N1 = 221.213596169931;
	static const double N2 = 112.079291497871;
	static const double N3 = 33.912866078383;
	static const double N4 = 6.37396220353165;
	static const double N5 = 0.700383064443688;
	static const double N6 = 3.52624965998911e-02;
	static const double M0 = 440.413735824752;
	static const double M1 = 793.826512519948;
	static const double M2 = 637.333633378831;
	static const double M3 = 296.564248779674;
	static const double M4 = 86.7807322029461;
	static const double M5 = 16.064177579207;
	static const double M6 = 1.75566716318264;
	static const double M7 = 8.83883476483184e-02;

	const double z = fabs(x);
	double c = 0.0;

	if (z <= 37.0)
	{
		const double e = exp(-z*z / 2.0);
		if (z<SPLIT)
		{
			const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
			const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
			c = e*n / d;
		}
		else
		{
			const double f = z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 13.0 / 20.0))));
			c = e / (RT2PI*f);
		}
	}
	return x <= 0.0 ? c : 1 - c;
}
//phi=cdf
double phi(double x)
{
	double L, K, w;
	/* constants */
	double  a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
	double  a4 = -1.821255978, a5 = 1.330274429;

	L = fabs(x);
	K = 1.0 / (1.0 + 0.2316419 * L);
	w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K, 3) + a4 * pow(K, 4) + a5 * pow(K, 5));

	if (x < 0) {
		w = 1.0 - w;
	}
	return w;
}

double pdf(double x) 
{
	double k = 0.398942280401433; // 1/sqrt(2*pi)
	double temp = x;
	return k  * exp(-0.5 * temp * temp);
}

double marginalLikelihood(int n_points, double** xtrain, int Dimension, double* params, double* ytrain, int type)
{
	double** id = idMatrix(n_points);
	double** noise = multConstant(id, 1e-12, n_points, n_points);

	double** K = KernelMV(xtrain, xtrain, n_points, n_points, Dimension, params, type);
	double** A = addMatrix(K, noise, n_points, n_points);
	double** Ainv = inverseMatrix(A, n_points);
	double* Ainvy = multVector(Ainv, ytrain, n_points, n_points);

	double yt_Ainvy = dotVector(ytrain, Ainvy, n_points);

	double logDeterminant = logDet(A, n_points);

	double marginal_likelihood = (-0.5*yt_Ainvy - 0.5*(logDeterminant) - 0.5 *log(2 * PI));
	free(Ainvy);
	freeMatrix(id, n_points);
	freeMatrix(noise, n_points);
	freeMatrix(K, n_points);
	freeMatrix(Ainv, n_points);
	freeMatrix(A, n_points);
	return marginal_likelihood;
}

double** getAinv(int n_points, double** xtrain, int Dimension, double* params, double* ytrain, int type)
{
	double** id = idMatrix(n_points);
	double** noise = multConstant(id, 1e-12, n_points, n_points);

	double** K = KernelMV(xtrain, xtrain, n_points, n_points, Dimension, params, type);
	double** A = addMatrix(K, noise, n_points, n_points);
	double** Ainv = inverseMatrix(A, n_points);

	freeMatrix(id, n_points);
	freeMatrix(noise, n_points);
	freeMatrix(A, n_points);
	freeMatrix(K, n_points);
	return Ainv;
}

double getMean(double* wi, int n_points, int Dimension, double* params,
	double** xtrain, double** points, int type)
{
	double** testpoint_kernel = KernelMV(xtrain, points, n_points, 1, Dimension, params, type);

	double** transtestpoint_kernel = transposeMatrix(testpoint_kernel, n_points, 1);

	double* mean_predTemp = multVector(transtestpoint_kernel, wi, 1, n_points);

	double mean_pred = mean_predTemp[0];

	freeMatrix(testpoint_kernel, n_points);
	free(mean_predTemp);
	freeMatrix(transtestpoint_kernel, 1);
	return mean_pred;
}

double getMeanWi(double** Ainv, double* ytrain, int n_points, int Dimension, double* params,
	double** xtrain, double** points, int type)
{
	double* wi = multVector(Ainv, ytrain, n_points, n_points);

	double** testpoint_kernel = KernelMV(xtrain, points, n_points, 1, Dimension, params, type);

	double** transtestpoint_kernel = transposeMatrix(testpoint_kernel, n_points, 1);

	double* mean_predTemp = multVector(transtestpoint_kernel, wi, 1, n_points);

	double mean_pred = mean_predTemp[0];

	free(wi);
	freeMatrix(testpoint_kernel, n_points);
	free(mean_predTemp);
	freeMatrix(transtestpoint_kernel, 1);
	return mean_pred;
}

double getVar(double** Ainv, int n_points, int Dimension, double* params, 
	double** xtrain, double** points, int type)
{
	double** testpoint_kernel = KernelMV(xtrain, points, n_points, 1, Dimension, params, type);

	double** transtestpoint_kernel = transposeMatrix(testpoint_kernel, n_points, 1);

	double** Ainvk = multMatrix(Ainv, testpoint_kernel, n_points, n_points, n_points, 1);

	double** ktAinvk = multMatrix(transtestpoint_kernel, Ainvk, 1, n_points, n_points, 1);

	//this is always one for one test point
	//double** k = KernelMV(points, points, 1, 1, Dimension, params, type);

	//double** ksubktAinvk = subMatrix(k, ktAinvk, 1, 1);

	//double var_pred = fabs(ksubktAinvk[0][0]);

	double var_pred = fabs(1 - ktAinvk[0][0])+ 1e-12;

	/*freeMatrix(ksubktAinvk, 1);
	freeMatrix(k, 1);*/
	freeMatrix(ktAinvk, 1);
	freeMatrix(testpoint_kernel, n_points);
	freeMatrix(transtestpoint_kernel, 1);
	freeMatrix(Ainvk, n_points);
	return var_pred;
}

double getErrorFunctionValue(double x)
{
	bool swap = (x < 0);

	if (swap)
	{
		x = -x;
	}

	const double p = 0.3275911;

	const double a1 = 0.254829592;
	const double a2 = -0.284496736;
	const double a3 = 1.421413741;
	const double a4 = -1.453152027;
	const double a5 = 1.061405429;

	double t1 = 1.0 / (1.0 + p * x);
	double t2 = t1 * t1;
	double t3 = t1 * t2;
	double t4 = t2 * t2;
	double t5 = t2 * t3;

	double value =
		1 - (a1 * t1 + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) *
		exp(-(x * x));

	return swap ? -value : value;
}

double getErrorFunctionDerivative(double x)
{
	return M_2_SQRTPI * exp(-(x * x));
}

double getCumulativeNormalValue(double x, double mean, double stdDev)
{
	double value = (x - mean) / (M_SQRT2 * stdDev);

	return 0.5 * (1 + getErrorFunctionValue(value));
}

double getNormalValue(double x, double mean, double stdDev)
{
	double sqrtTwoPi = 2.50662827463;

	double temp = x - mean;

	double dummyVar= exp(-(temp * temp) / (2 * stdDev * stdDev)) / (sqrtTwoPi * stdDev);
	return dummyVar;
}

// Determine average.

double averagePop(double* y, int size)
{
	double average = 0;
	double averageFactor = 1.0 / size;

	for (int i = size - 1; i >= 0; --i)
	{
		average += y[i] * averageFactor;
	}

	return average;
}

double standardDeviationPop(double* y, int size, double average, double averageFactor)
{
	// Determine standard deviation.
	double standardDeviation = 0;

	for (int i = size - 1; i >= 0; --i)
	{
		const double diff = y[i] - average;

		standardDeviation += diff * diff * averageFactor;
	}

	standardDeviation = sqrt(standardDeviation);

	return standardDeviation;
}


void moveY(double average , double standardDeviation, double* y, int size)
{
	// Scale data.
	double offset = average;
	double scaleFactor = 2 * standardDeviation;

	for (int i = size - 1; i >= 0; --i)
	{
		y[i] = (y[i] - offset) / scaleFactor;
	}
}

//getCumulativeNormalValue();
double acqFunctionEI2Kernels(double** AinvKexp, double** AinvKM52, double* wiKexp, double* wiKm52, int n_x,
	int number_of_parameters, double* paramKexp, double*paramKM52,  double** xtrain, 
	double** testpoint, double bestSolutionA)
{
	double	mean_predKexp = getMean(wiKexp, n_x, number_of_parameters, paramKexp, xtrain, testpoint, 1);

	double	var_predKexp = getVar(AinvKexp, n_x, number_of_parameters, paramKexp, xtrain, testpoint, 1);

	double	mean_predKM52 = getMean(wiKm52, n_x, number_of_parameters, paramKM52, xtrain, testpoint, 0);

	double	var_predKM52 = getVar(AinvKM52, n_x, number_of_parameters, paramKM52, xtrain, testpoint, 0);

	double mean_pred = (mean_predKexp + mean_predKM52) / 2;
	double var_pred = (var_predKexp + var_predKM52) / 2;

	double result;
	if (var_pred < abs(1e-12))
	{
		result = 0;
	}
	else
	{
		double Z = (bestSolutionA - mean_pred) / sqrt(var_pred);

		result = (bestSolutionA - mean_pred)*cdf(Z) + sqrt(var_pred)*pdf(Z);

		result = -result;

		/*if (cdf(Z) < abs(1e-8) && pdf(Z)< abs(1e-8))
		{
		double meanPop = averagePop(ytrain, n_x);
		double averageFactor = 1.0 / n_x;
		double standardPop = standardDeviationPop(ytrain, n_x, meanPop, averageFactor);

		double map = getCumulativeNormalValue(mean_pred, meanPop, standardPop);
		result = map * pow(var_pred, -power);
		ind->objective_value = result;
		}*/
	}

	/*if (result == 0)
	{
	double meanPop = averagePop(ytrain, n_x);
	double averageFactor = 1.0 / n_x;
	double standardPop = standardDeviationPop(ytrain, n_x, meanPop, averageFactor);

	double map = getCumulativeNormalValue(mean_pred, meanPop, standardPop);
	result = map * pow(var_pred, -power);
	ind->objective_value = result;
	}*/
	return result;
}

//getCumulativeNormalValue();
double acqFunctionUCB2Kernels(double** AinvKexp, double** AinvKM52, double* wiKexp, double* wiKm52, int n_x,
	int number_of_parameters, double* paramKexp, double*paramKM52, double** xtrain,
	double** testpoint, double kappa)
{
	double	mean_predKexp = getMean(wiKexp, n_x, number_of_parameters, paramKexp, xtrain, testpoint, 1);

	double	var_predKexp = getVar(AinvKexp, n_x, number_of_parameters, paramKexp, xtrain, testpoint, 1);

	double	mean_predKM52 = getMean(wiKm52, n_x, number_of_parameters, paramKM52, xtrain, testpoint, 0);

	double	var_predKM52 = getVar(AinvKM52, n_x, number_of_parameters, paramKM52, xtrain, testpoint, 0);


	double mean_pred = (mean_predKexp + mean_predKM52) / 2;
	double var_pred = (var_predKexp + var_predKM52) / 2;

	double result = mean_pred - kappa*sqrt(var_pred);

	return result;
}

double acqFunctionEI(double** Ainv,  double* wi, int n_x,
	int dimension, double* params, double** xtrain,
	double** testpoint, double bestSolutionA, int type)
{
	double mean_pred, var_pred;

		mean_pred = getMean(wi, n_x, dimension, params, xtrain, testpoint, type);

		var_pred = getVar(Ainv, n_x, dimension, params, xtrain, testpoint, type);


	double Z = (bestSolutionA - mean_pred);
	//double Z = (mean_pred - bestSolutionA);
	double result;
	if (var_pred ==0)
	{
		result = min(0,Z);
	}
	else
	{
		double tempZ=Z / sqrt(var_pred);

		result = Z*cdf(tempZ) + sqrt(var_pred)*pdf(tempZ);

		
	}
	result = -result;

	return result;
}

double acqFunctionUCB(double** Ainv, double* wi, int n_x,
	int dimension, double* params, double** xtrain,
	double** testpoint, double kappa, int type)
{
	double mean_pred, var_pred;
	
		mean_pred = getMean(wi, n_x, dimension, params, xtrain, testpoint, type);

		var_pred = getVar(Ainv, n_x, dimension, params, xtrain, testpoint, type);

		/*printf("mean \t %0.8f \n", mean_pred);
		printf("var \t %0.8f \n", var_pred);*/

	double result = mean_pred - kappa*sqrt(var_pred);

	for (int i = 0; i < dimension; i++)
	{
		if (testpoint[0][i] > 1)
		{
			result = 1e20;
		}
		if (testpoint[0][i] < 0)
		{
			result = 1e20;
		}
	}
	//printf("result \t %0.8f \n", result);

	return result;
}

double* getKseV(int n_points, int Dimension, double* params,
	double** xtrain, double** points, int type)
{
	double* KseV = newDoubleV(n_points);
	double** testpoint_kernel = KernelMV(xtrain, points, n_points, 1, Dimension, params, type);

	for (int i = 0; i < n_points; i++)
	{
		KseV[i] = testpoint_kernel[i][0];
	}

	freeMatrix(testpoint_kernel, n_points);
	return KseV;

}



double gradexponentialKernel(double** xKernel, double** yKernel, int indexI, int indexJ, int Dimension, double* k_param)
{
	double* xI = newDoubleV(Dimension);
	double* xJ = newDoubleV(Dimension);
	int i;

	for (i = 0; i < Dimension; i++)
	{
		xI[i] = xKernel[indexI][i];
		xJ[i] = yKernel[indexJ][i];
	}

	double* temp = subVector(xI, xJ, Dimension);

	double** id = idMatrix(Dimension);
	for (i = 0; i < Dimension; i++)
	{
		id[i][i] = 1 / (exp(-(k_param[i] * k_param[i])));
	}

	double* tempVector = multVector(id, temp, Dimension, Dimension);

	double norm = dotVector(temp, tempVector, Dimension);

	double r = sqrt(norm);

	free(temp);
	free(xI);
	free(xJ);
	free(tempVector);
	for (i = 0; i < Dimension; i++)
	{
		free(id[i]);
	}
	free(id);


	double exponential_kernel = exp(-(k_param[Dimension] * k_param[Dimension]))  * exp(-0.5 * norm);

	return exponential_kernel;
}


double** gradientCov(double** xKernel, double** x2Kernel, int size_x, int size_x2, int Dimension, 
	double* k_param)
{
	
	double** cov = newDouble(size_x, size_x2);
	int i, j;

	for (i = 0; i< size_x; i++)
	{
		for (j = 0; j< size_x2; j++)
		{
				cov[i][j] = gradexponentialKernel(xKernel, x2Kernel, i, j, Dimension, k_param);
		}
	}


	return cov;
}


double* derivateKseV(int n_points, int Dimension, double* params,
	double** xtrain, double** points, int type, int variableIndex)
{
	double** testpoint_kernel = gradientCov(xtrain, points, n_points, 1, Dimension, params);

	double* Kse= newDoubleV(n_points);

	for (int i = 0; i < n_points; i++)
	{
		Kse[i] = testpoint_kernel[i][0];
	}

	double* distances = newDoubleV(n_points);

	for (int i = 0; i < n_points; i++)
	{
		distances[i] = -(points[0][variableIndex]- xtrain[i][variableIndex])/
			exp(-(params[variableIndex] * params[variableIndex]));
	}

	double* dKsev=newDoubleV(n_points);
	
	for (int i = 0; i < n_points; i++)
	{
		dKsev[i] = distances[i] * Kse[i];
	}

	free(distances);
	free(Kse);
	freeMatrix(testpoint_kernel, n_points);
	
	return dKsev;

}

 double getDerUCB(double** Ainv, double* ytrain, int n_points, int Dimension, double* params,
	double** xtrain, double** points, int type, int variableIndex, double kappa)
{
	double* wi = multVector(Ainv, ytrain, n_points, n_points);

	double* dKse = derivateKseV(n_points, Dimension, params,
		xtrain, points, type, variableIndex);

	double var_pred = getVar(Ainv, n_points, Dimension, 
		params, xtrain, points, type);

	double* Kse = getKseV(n_points, Dimension, params,
		xtrain, points, type);


	double* AinvKse = multVector(Ainv, Kse, n_points, n_points);
	
	double coef1 = dotVector(dKse, AinvKse, n_points);

	double* AinvdKse= multVector(Ainv, dKse, n_points, n_points);

	double coef2 = dotVector(Kse, AinvdKse, n_points);
	
	double dmean = dotVector(wi, dKse, n_points);

	double dvarsqr = -(coef1 + coef2);

	//printf("dmean  \t %0.4f \n", dmean);
	//printf("dvarsqr  \t %0.4f \n", dvarsqr);
	//printf("var_pred  \t %0.4f \n", var_pred);
	double result = 0;
	if (fabs(var_pred) > 0)
	{
		double dvar =  pow(var_pred,-(1/2))/2*dvarsqr;
		result = dmean - kappa*dvar;
	}
	else
	{
		result = 0;
	}
	//printf("result  \t %0.4f \n", result);
	//printf("%d \t %0.4f \n", variableIndex, dmean);
	free(AinvdKse);
	free(AinvKse);
	free(Kse);
	free(dKse);
	free(wi);

	return result;
}
