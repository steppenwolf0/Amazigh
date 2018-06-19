/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "../Amazigh/mathVector.h"
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/openFile.h"
#include "../Amazigh/showFile.h"
#include "../Amazigh/printFile.h"
#include "../Amazigh/supportC.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "../util/problem.h"
#include <math.h>
#include <stdlib.h>
double** transferMatrix;
int dim_1;
int dim_2;
double* b;
double** nodesMatrix;

void initializeproblem(){
	struct matDouble transfer=openMatrixDouble("Transfer.matrix");
	transferMatrix=transfer.matrix;
	dim_1=transfer.dim_1;
	dim_2=transfer.dim_2;
	printf("dim_1 %d\n", dim_1);
	printf("dim_2 %d\n", dim_2);
	
	struct matDouble ut=openMatrixDouble("ut.matrix");
	b=getVector(ut.matrix,85, dim_1);
	printf("dim_1 %d\n", ut.dim_1);
	printf("dim_2 %d\n", ut.dim_2);

	
	struct matDouble Nodes = openMatrixDouble("nodes.matrix");
	nodesMatrix = Nodes.matrix;
	printf("dim_1 %d\n", Nodes.dim_1);
	printf("dim_2 %d\n", Nodes.dim_2);

	int i, j;
	for (i = 0; i < Nodes.dim_1; i++)
	{
		for (j= 0; j < Nodes.dim_2; j++)
		{
			nodesMatrix[i][j] = nodesMatrix[i][j] / 1000.0;
		}

	}

	//printVector(b,dim_1,"b.txt");
	number_of_parameters = 5;
}

void FunctionProblemEvaluation( Individual *ind, int number_of_parameters  )
{
    int    i;
    double result;

	double* x = newDoubleV(dim_2);
	double* y = newDoubleV(dim_2);
	double* z = newDoubleV(dim_2);



	for (i = 0; i < dim_2; i++)
	{
		x[i] = (1.0/fabs(nodesMatrix[i][0] - ind->parameters[0]))*(1.0 / fabs(nodesMatrix[i][0] - ind->parameters[0]));
		y[i] = (1.0/fabs(nodesMatrix[i][1] - ind->parameters[1]))*(1.0 / fabs(nodesMatrix[i][1] - ind->parameters[1]));
		z[i] = (1.0/fabs(nodesMatrix[i][2] - ind->parameters[2]))*(1.0 / fabs(nodesMatrix[i][2] - ind->parameters[2]));
	}


		
	double* laplace = newDoubleV(dim_2);
	for (i = 0; i < dim_2; i++)
	{
		laplace[i] = (ind->parameters[4])*sqrt(x[i]) * sqrt(y[i]) * sqrt(z[i])+ ind->parameters[3];
	}

	
	double* temp=multVector(transferMatrix, laplace,dim_1,dim_2);
	double* subTemp=subVector(b,temp,dim_1);
	
    result=dotVector(subTemp,subTemp,dim_1);
	
	
  /*  result = 0.0;
    for( i = 0; i < number_of_parameters; i++ )
        result += ind->parameters[i]*ind->parameters[i];*/
        
	for (i = 0; i < dim_2; i++)
		ind->result[i]= laplace[i];
	
	ind->objective_value  = result;
	ind->constrained_value = 0;

	free(x);
	free(y);
	free(z);



	free(subTemp);


	free(temp);

	free(laplace);
}

void FunctionPartialProblemEvaluation( Individual* ind, int number_of_touched_parameters, int *touched_parameters_indices, 
double *touched_parameters, double *parameters_before, double objective_value_before, 
double constraint_value_before )
{
    int    i,j;
    double result;

  /*  result = objective_value_before;
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
        result += touched_parameters[i]*touched_parameters[i];
        result -= parameters_before[i]*parameters_before[i];
    }*/
    
    for (i=0;i<number_of_touched_parameters;i++)
    {
    	double constant=touched_parameters[i];
    	double constant2=parameters_before[i];
		
		if (touched_parameters_indices[i] == 3)
		{
			double* laplace = newDoubleV(dim_2);
			for (j = 0; j < dim_2; j++)
			{
				laplace[i] = ind->result[i] - constant2+constant;
			}

			double* temp = multVector(transferMatrix, laplace, dim_1, dim_2);
			double* subTemp = subVector(b, temp, dim_1);

			result = dotVector(subTemp, subTemp, dim_1);

			for (j = 0; j < dim_2; j++)
				ind->result[j] = laplace[j];
			
			free(laplace);
			free(temp);
			free(subTemp);

			ind->objective_value = result;
			ind->constrained_value = 0;
		}
		if (touched_parameters_indices[i] == 2)
		{
			double* x = newDoubleV(dim_2);
			double* y = newDoubleV(dim_2);
			double* z = newDoubleV(dim_2);

			for (j = 0; j < dim_2; j++)
			{
				x[j] = (1.0 / fabs(nodesMatrix[j][0] - ind->parameters[0]))*(1.0 / fabs(nodesMatrix[j][0] - ind->parameters[0]));
				y[j] = (1.0 / fabs(nodesMatrix[j][1] - ind->parameters[1]))*(1.0 / fabs(nodesMatrix[j][1] - ind->parameters[1]));
				z[j] = (1.0 / fabs(nodesMatrix[j][2] - constant))*(1.0 / fabs(nodesMatrix[j][2] - constant));
			}



			double* laplace = newDoubleV(dim_2);
			for (j = 0; j < dim_2; j++)
			{
				laplace[j] = -(1.0 / 4.0)*sqrt(x[j]) * sqrt(y[j]) * sqrt(z[j]) + ind->parameters[3];
			}


			double* temp = multVector(transferMatrix, laplace, dim_1, dim_2);
			double* subTemp = subVector(b, temp, dim_1);

			result = dotVector(subTemp, subTemp, dim_1);

			for (j = 0; j < dim_2; j++)
				ind->result[j] = laplace[j];

			free(x);
			free(y);
			free(z);
			free(laplace);
			free(temp);
			free(subTemp);

			ind->objective_value = result;
			ind->constrained_value = 0;
		}
		if (touched_parameters_indices[i] == 1)
		{
			double* x = newDoubleV(dim_2);
			double* y = newDoubleV(dim_2);
			double* z = newDoubleV(dim_2);

			for (j = 0; j < dim_2; j++)
			{
				x[j] = (1.0 / fabs(nodesMatrix[j][0] - ind->parameters[0]))*(1.0 / fabs(nodesMatrix[j][0] - ind->parameters[0]));
				y[j] = (1.0 / fabs(nodesMatrix[j][1] - constant))*(1.0 / fabs(nodesMatrix[j][1] - constant));
				z[j] = (1.0 / fabs(nodesMatrix[j][2] - ind->parameters[2]))*(1.0 / fabs(nodesMatrix[j][2] - ind->parameters[2]));
			}



			double* laplace = newDoubleV(dim_2);
			for (j = 0; j < dim_2; j++)
			{
				laplace[j] = -(1.0 / 4.0)*sqrt(x[j]) * sqrt(y[j]) * sqrt(z[j]) + ind->parameters[3];
			}


			double* temp = multVector(transferMatrix, laplace, dim_1, dim_2);
			double* subTemp = subVector(b, temp, dim_1);

			result = dotVector(subTemp, subTemp, dim_1);

			for (j = 0; j < dim_2; j++)
				ind->result[j] = laplace[j];

			free(x);
			free(y);
			free(z);
			free(laplace);
			free(temp);
			free(subTemp);

			ind->objective_value = result;
			ind->constrained_value = 0;
		}
		if (touched_parameters_indices[i] == 0)
		{
			double* x = newDoubleV(dim_2);
			double* y = newDoubleV(dim_2);
			double* z = newDoubleV(dim_2);

			for (j = 0; j < dim_2; j++)
			{
				x[j] = (1.0 / fabs(nodesMatrix[j][0] - constant))*(1.0 / fabs(nodesMatrix[j][0] - constant));
				y[j] = (1.0 / fabs(nodesMatrix[j][1] - ind->parameters[1]))*(1.0 / fabs(nodesMatrix[j][1] - ind->parameters[1]));
				z[j] = (1.0 / fabs(nodesMatrix[j][2] - ind->parameters[2]))*(1.0 / fabs(nodesMatrix[j][2] - ind->parameters[2]));
			}



			double* laplace = newDoubleV(dim_2);
			for (j = 0; j < dim_2; j++)
			{
				laplace[j] = -(1.0 / 4.0)*sqrt(x[j]) * sqrt(y[j]) * sqrt(z[j]) + ind->parameters[3];
			}


			double* temp = multVector(transferMatrix, laplace, dim_1, dim_2);
			double* subTemp = subVector(b, temp, dim_1);

			result = dotVector(subTemp, subTemp, dim_1);

			for (j = 0; j < dim_2; j++)
				ind->result[j] = laplace[j];

			free(x);
			free(y);
			free(z);
			free(laplace);
			free(temp);
			free(subTemp);

			ind->objective_value = result;
			ind->constrained_value = 0;
		}
	}

	
   
}

double sphereFunctionProblemLowerRangeBound( int dimension )
{
    return( -1e+7 );
}

double sphereFunctionProblemUpperRangeBound( int dimension )
{
    return( 1e+7 );
}

void copyIndividual(Individual *ind2, Individual *ind1, int number_of_parameters)
{
	int j;
	for (j = 0; j < number_of_parameters; j++)
		ind2->parameters[j] = ind1->parameters[j];

	for (j = 0; j < dim_2; j++)
		ind2->result[j] = ind1->result[j];

	ind2->objective_value = ind1->objective_value;
	ind2->constrained_value = ind1->constrained_value;
}

void initializeIndividuals(int maximum_number_of_populations)
{
	populations = (Individual ***)Malloc(maximum_number_of_populations*sizeof(Individual **));
	selections = (Individual ***)Malloc(maximum_number_of_populations*sizeof(Individual **));

	
}

void initializePopulations(int population_index, int number_of_parameters, int* population_sizes, int* selection_sizes)
{
	int j;
	
	populations[population_index] = (Individual **)Malloc(population_sizes[population_index] * sizeof(Individual*));
	for (j = 0; j < population_sizes[population_index]; j++)
		populations[population_index][j] = initializeIndividual(number_of_parameters);

	selections[population_index] = (Individual **)Malloc(selection_sizes[population_index] * sizeof(Individual*));
	for (j = 0; j < selection_sizes[population_index]; j++)
		//selections[population_index][j].parameters = (double *)Malloc(number_of_parameters*sizeof(double));
		selections[population_index][j] = initializeIndividual(number_of_parameters);
}

//Individual* declareIndividual(double* parameters, double objective_value, double contrained_value, int number_of_parameters)
//{
//	Individual* ind;
//	ind->parameters=(double *)Malloc(number_of_parameters*sizeof(double));
//	int i;
//	for (i = 0; i < number_of_parameters; i++)
//	{
//		ind->parameters[i] = parameters[i];
//	}
//	ind->constrained_value = contrained_value;
//	ind->objective_value = objective_value;
//
//
//	return ind;
//}
//
//Individual* initializeIndividual(int number_of_parameters)
//{
//	Individual *ind=(Individual*)malloc(sizeof(Individual));
//	ind->parameters = (double *)Malloc(number_of_parameters*sizeof(double));
//	int i;
//	for (i = 0; i < number_of_parameters; i++)
//	{
//		ind->parameters[i] = 0.0;
//	}
//	ind->constrained_value = 0.0;
//	ind->objective_value = 0.0;
//
//
//	return ind;
//}

Individual* declareIndividual(double* parameters, double objective_value, double contrained_value, int number_of_parameters,
	double* result)
{
	Individual* ind = (Individual *)Malloc(sizeof(Individual));
	ind->parameters = (double *)Malloc(number_of_parameters*sizeof(double));
	ind->result = (double *)Malloc(dim_2*sizeof(double));
	int i;
	for (i = 0; i < number_of_parameters; i++)
	{
		ind->parameters[i] = parameters[i];
	}
	ind->constrained_value = contrained_value;
	ind->objective_value = objective_value;
	for (i = 0; i < dim_2; i++)
	{
		ind->result[i] = result[i];
	}

	return ind;
}

Individual* initializeIndividual(int number_of_parameters)
{
	Individual* ind = (Individual *)Malloc(sizeof(Individual));
	ind->parameters = (double *)Malloc(number_of_parameters*sizeof(double));
	ind->result = (double *)Malloc(dim_2*sizeof(double));
	int i;
	for (i = 0; i < number_of_parameters; i++)
	{
		ind->parameters[i] = 0.0;
	}
	ind->constrained_value = 0.0;
	ind->objective_value = 0.0;

	for (i = 0; i < dim_2; i++)
	{
		ind->result[i] = 0.0;
	}

	return ind;
}


void freeIndividual(Individual* ind)
{
	if (ind)
	{
		free(ind->parameters);
		free(ind->result);
		free(ind);
	}
	
	
}
