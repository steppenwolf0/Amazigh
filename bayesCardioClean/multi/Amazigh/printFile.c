#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
void printMatrix(double** matrix, int n, int m, const char* name)
        {
        	FILE *f = fopen(name, "w");
			if (f == NULL)
			{
			    printf("Error opening file!\n");
			    exit(1);
			}
        	
        	int i,j;
            fprintf(f,"%d\n",n);
            fprintf(f,"%d\n",m);
            for ( i = 0; i < n; i++)
            {
                
                for ( j = 0; j < m; j++)
                {
                    fprintf(f,"%.18f ",matrix[i][j]);
                }
                fprintf(f,"\n");
                
            }
            fclose(f);
        }
void printVector(double* vector, int n, const char* name)
        {
        	FILE *f = fopen(name, "w");
			if (f == NULL)
			{
			    printf("Error opening file!\n");
			    exit(1);
			}
        	
        	int i;
            fprintf(f,"%d\n",n);
            for ( i = 0; i < n; i++)
            {

                fprintf(f,"%.18f",vector[i]);
                fprintf(f,"\n");
                
            }
            fclose(f);
        }

void printMatrixInt(int** matrix, int n, int m, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i, j;
	fprintf(f, "%d\n", n);
	fprintf(f, "%d\n", m);
	for (i = 0; i < n; i++)
	{

		for (j = 0; j < m; j++)
		{
			fprintf(f, "%d ", matrix[i][j]);
		}
		fprintf(f, "\n");

	}
	fclose(f);
}

void printVectorInt(int* vector, int n, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i;
	fprintf(f, "%d\n", n);
	for (i = 0; i < n; i++)
	{

		fprintf(f, "%d", vector[i]);
		fprintf(f, "\n");

	}
	fclose(f);
}