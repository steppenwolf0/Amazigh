#include "problemDef.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/mathVector.h"
#include "../Gaussian/gaussian.h"
#ifdef OSisWindows
#include "./RV-GOMEA/RV-GOMEAWIN.h"
#else
#include "./RV-GOMEA/RV-GOMEA.h"
#endif
int maxEval ;

Individual* RVGOMEAUCB(double** Ainv, double** xtrain,double* param,double* ytrain,
	int n_param,int n_x, int Dimension, double kappa, int kernelType, int evals)
{
	maxEval = evals;
	type = kernelType;
	problem_Solving = 0;
	initializeproblemUCB( Ainv,  xtrain,  param,  ytrain, n_param, n_x, Dimension, kappa, kernelType);
	
	//-b -w -s -r -f 1 -1e2 1e2 0.35 10 25 0.9 1 0 1e-10 100 0.0  10000
	const char** argv[19];
	//argv[17]= { "","-s", "-r", "-f", "1", "-1e2", "1e2", "0.35", "10", "25", "0.9", "1", "0", "1e-10", "100", "0.0", "10" };
	argv[0] = " ";
	argv[1] = "-b";
	argv[2] = "-b";
	//argv[3] = "-s"; generational statistics
	argv[3] = "-b";
	argv[4] = "-b";
	argv[5] = "-f";
	argv[6] = "1";
	argv[7] = "0";
	argv[8] = "1";
	argv[9] = "0.35";
	argv[10] = "200"; //play with this one
	argv[11] = "5";
	argv[12] = "0.9";
	argv[13] = "1";
	argv[14] = "0";
	argv[15] = "-1e8"; //without -r this does not matter
	argv[16] = "100";
	argv[17] = "0.0";
	argv[18] = "30";

	//noError = noError && sscanf(argv[*index + 0], "%lf", &lower_user_range);
	//noError = noError && sscanf(argv[*index + 1], "%lf", &upper_user_range);
	//noError = noError && sscanf(argv[*index + 2], "%lf", &tau);
	//noError = noError && sscanf(argv[*index + 3], "%d", &base_population_size);
	//noError = noError && sscanf(argv[*index + 4], "%d", &maximum_number_of_populations);
	//noError = noError && sscanf(argv[*index + 5], "%lf", &distribution_multiplier_decrease);
	//noError = noError && sscanf(argv[*index + 6], "%lf", &st_dev_ratio_threshold);
	//noError = noError && sscanf(argv[*index + 7], "%lf", &maximum_number_of_evaluations);
	//noError = noError && sscanf(argv[*index + 8], "%lf", &vtr);
	//noError = noError && sscanf(argv[*index + 9], "%d", &maximum_no_improvement_stretch);
	//noError = noError && sscanf(argv[*index + 10], "%lf", &fitness_variance_tolerance);
	//noError = noError && sscanf(argv[*index + 11], "%lf", &maximum_number_of_seconds);
	
	int argc = 19;
	interpretCommandLine(argc, argv);
	
	Individual* ouput = run();

	closeProblem();
	
	return ouput;
}

Individual* RVGOMEAEI(double** Ainv, double** xtrain, double* param, double* ytrain, int n_param, int n_x, int Dimension, double bestSolution, int kernelType)
{
	type = kernelType;
	problem_Solving = 2;
	initializeproblemEI(Ainv, xtrain, param, ytrain, n_param, n_x, Dimension, bestSolution);

	
	//-b -w -s -r -f 1 -1e2 1e2 0.35 10 25 0.9 1 0 1e-10 100 0.0  10000
	const char** argv[19];
	//argv[17]= { "","-s", "-r", "-f", "1", "-1e2", "1e2", "0.35", "10", "25", "0.9", "1", "0", "1e-10", "100", "0.0", "10" };
	argv[0] = " ";
	argv[1] = "-b";
	argv[2] = "-b";
	//argv[3] = "-s"; generational statistics
	argv[3] = "-b";
	argv[4] = "-r";
	argv[5] = "-f";
	argv[6] = "1";
	argv[7] = "-1e2";
	argv[8] = "1e2";
	argv[9] = "0.35";
	argv[10] = "10";
	argv[11] = "25";
	argv[12] = "0.9";
	argv[13] = "1";
	argv[14] = "0";
	argv[15] = "-1e4";
	argv[16] = "100";
	argv[17] = "0.0";
	argv[18] = "10";

	

	/*printf("\n");
	printf("# Tau                     = %e\n", tau);
	printf("# Population size/normal  = %d\n", base_population_size);
	printf("# FOS element size        = %d\n", FOS_element_size);
	printf("# Max num of populations  = %d\n", maximum_number_of_populations);
	printf("# Dis. mult. decreaser    = %e\n", distribution_multiplier_decrease);
	printf("# St. dev. rat. threshold = %e\n", st_dev_ratio_threshold);
	printf("# Maximum numb. of eval.  = %lf\n", maximum_number_of_evaluations);
	printf("# Value to reach (vtr)    = %e\n", vtr);
	printf("# Max. no improv. stretch = %d\n", maximum_no_improvement_stretch);
	printf("# Fitness var. tolerance  = %e\n", fitness_variance_tolerance);
	printf("# Random seed             = %ld\n", random_seed);
	printf("#\n");*/

	int argc = 19;
	interpretCommandLine(argc, argv);

	Individual* ouput = run();

	closeProblem();

	return ouput;
}

