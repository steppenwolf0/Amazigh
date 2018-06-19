#ifndef PROBLEM_H
#define PROBLEM_H

typedef struct Individual{
	double constrained_value;
	double objective_value;
	double* parameters;
}Individual;

Individual ***populations;
Individual ***selections;
int number_of_parameters;

void copyIndividual(Individual *ind1, Individual* ind2, int number_of_parameters);
Individual* initializeIndividual(int number_of_parameters);
void FunctionProblemEvaluation(Individual *ind, int number_of_parameters  );
void FunctionPartialProblemEvaluation(Individual *ind, int number_of_touched_parameters, int *touched_parameters_indices, 
		double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );
double functionProblemLowerRangeBound( int dimension );
double functionProblemUpperRangeBound( int dimension );
void initializePopulations(int population_index, int number_of_parameters, int* population_sizes, int* selection_sizes);
void freeIndividual(Individual* ind);

//Problem Specific********************************************************************************
void initializeproblemUCB(double** _Ainv, double** _xtrain, double* _param, double* _ytrain, int _n_param, int _n_x,
	int Dimension, double _Kappa, int _type);
void initializeproblemEI(double** _Ainv, double** _xtrain, double* _param, double* _ytrain, int _n_param, int _n_x,
	int _Dimension, double _bestSolution);
void initializeproblem2A(double** _AinvKexp, double** _AinvKM52, double** _xtrain, double* _paramKexp,
	double* _paramKM52, double* _ytrain, int _n_param, int _n_x, int _Dimension, double _bestSolution, double power);
void initializeproblem0A(double** _AinvKexp, double** _AinvKM52, double** _xtrain, double* _paramKexp, double* _paramKM52,
	double* _ytrain, int _n_param, int _n_x, int _Dimension, double _Kappa);
void initializeIndividuals(int maximum_number_of_populations);
void closeProblem();
void closeProblemA();

int problem_Solving;
int type;
//Problem Specific********************************************************************************
#endif
