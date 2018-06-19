#include <stdlib.h>
double** newDouble(int dim_1, int dim_2)
{
	double** matrixdouble =  malloc( dim_1*( sizeof( double* ) ) );
    for( int i = 0; i < dim_1; i++ )
        matrixdouble[i] =  malloc( dim_2*( sizeof( double ) ) );
        
    for ( int i=0;i<dim_1;i++)
            {
                for ( int j = 0; j < dim_2; j++)
                {
                	matrixdouble[i][j]=0.0;
                }
            }
    return matrixdouble;                	
}


double* newDoubleV(int dim)
{
	double* vector = (double*) malloc( dim*( sizeof( double ) ) );
        
        for ( int i=0;i<dim;i++)
            {
            	vector[i]=0.0;
            }
        return vector;
}

double* newDoubleVInit(int dim,double constant)
{
	double* vector = (double*)malloc(dim*(sizeof(double)));

	for (int i = 0; i<dim; i++)
	{
		vector[i] = constant;
	}
	return vector;
}


int* newIntV(int dim)
{
	int* vector = (int*)malloc(dim*(sizeof(int)));

	for (int i = 0; i<dim; i++)
	{
		vector[i] = 0;
	}
	return vector;
}

void freeMatrix(double** matrix, int dim)
{
	int i;
	for (i = 0; i < dim; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void copyVector(double* vector1, double* vector2, int size)
{
	int j;
	for (j = 0; j < size; j++)
	{
		vector1[j] = vector2[j];

	}
}

void copyMatrix(double** matrix1, double** matrix2, int size1, int size2)
{
	int i, j;
	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			matrix1[i][j] = matrix2[i][j];

		}
	}
}

int** newInt(int dim_1, int dim_2)
{
	int** matrixint = (int**) malloc(dim_1*sizeof(int*));
	for (int i = 0; i < dim_1; i++)
		matrixint[i] = (int*)malloc(dim_2*sizeof(int));

	for (int i = 0; i<dim_1; i++)
	{
		for (int j = 0; j < dim_2; j++)
		{
			matrixint[i][j] = 0;
		}
	}
	return matrixint;
}

void freeMatrixInt(int** matrix, int dim)
{
	int i;
	for (i = 0; i < dim; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void copyMatrixInt(int** matrix1, int** matrix2, int size1, int size2)
{
	for (int i = 0; i < size1; i++)
	{
		for (int j = 0; j < size2; j++)
		{
			matrix1[i][j] = matrix2[i][j];

		}
	}
}

//USAGE
//int size = 0;
//int** List;
//size = listInt2DExists(&List, 0, 5, size);
//printf("size %d \n", size);
//showMatrixInt(List, size, 2);


int listInt2D(int*** List, int index0, int index1, int size)
{
	if (size != 0)
	{
		int** tempList = newInt(size, 2);
		copyMatrixInt(tempList, *List, size, 2);

		freeMatrixInt(*List, size);
		int** tempPointer = newInt(size + 1, 2);
		tempPointer[size][0] = index0;
		tempPointer[size][1] = index1;
		copyMatrixInt(tempPointer, tempList, size, 2);
		freeMatrixInt(tempList, size);

		*List = tempPointer;
	}
	else
	{
		int** tempList = newInt(1, 2);
		tempList[0][0] = index0;
		tempList[0][1] = index1;
		*List = tempList;



	}


	return (size + 1);
}

int listInt2DExists(int*** List, int index0, int index1, int size)
{
	int flag = 0;
	int temp = size;
	if (size != 0)
	{
		for (int i = 0; i < size; i++)
		{
			int l0 = List[0][i][0];
			int l1 = List[0][i][1];

			if (l0 == index0 && l1 == index1)
				flag = 1;
		}
		if (flag == 0)
		{
			int** tempList = newInt(size, 2);
			copyMatrixInt(tempList, *List, size, 2);

			freeMatrixInt(*List, size);
			int** tempPointer = newInt(size + 1, 2);
			tempPointer[size][0] = index0;
			tempPointer[size][1] = index1;
			copyMatrixInt(tempPointer, tempList, size, 2);
			freeMatrixInt(tempList, size);

			*List = tempPointer;
			temp = size + 1;
		}
		
	}
	else
	{
		int** tempList = newInt(1, 2);
		tempList[0][0] = index0;
		tempList[0][1] = index1;
		*List = tempList;

		temp = 1;

	}


	return temp;
}


int** newIntValue(int dim_1, int dim_2, int value)
{
	int** matrixint = (int**)malloc(dim_1*sizeof(int*));
	for (int i = 0; i < dim_1; i++)
		matrixint[i] = (int*)malloc(dim_2*sizeof(int));

	if (matrixint == NULL)
	{
		printf("Memory allocation failed");
		return;
	}

	for (int i = 0; i<dim_1; i++)
	{
		for (int j = 0; j < dim_2; j++)
		{
			matrixint[i][j] = value;
		}
	}
	return matrixint;
}


void copyVectorConst(double* vector1, double* vector2, int size, double constant)
{
	int j;
	for (j = 0; j < size; j++)
	{
		vector1[j] = vector2[j]*constant;

	}
}
