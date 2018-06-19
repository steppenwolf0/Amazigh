#include <stdio.h>
#include <string.h>

void showMatrix(double** matrix, int n, int m)
        {
        	int i,j;
            
            for ( i = 0; i < n; i++)
            {
                
                for ( j = 0; j < m; j++)
                {
                    printf("%.6f\t" ,matrix[i][j]);
                }
                printf("\n");
                
            }
            
        }
        
void showMatrixInt(int** matrix, int n, int m)
        {
        	            
            for ( int i = 0; i < n; i++)
            {
                
                for ( int j = 0; j < m; j++)
                {
                    printf("%d ",matrix[i][j]);
                }
                printf("\n");
                
            }
            
        }

void showDouble(double value)
        {
            printf("%.4f ",value);
            printf("\n");
        }

void showVector(double* vector, int n)
{
	int j;

	for (j = 0; j < n; j++)
		{
			printf("%.6f\t ", vector[j]);
		}
	printf("\n");

	

}

void showVectorInt(int* vector, int n)
{
	int j;

	for (j = 0; j < n; j++)
	{
		printf("%d\t ", vector[j]);
	}
	printf("\n");



}