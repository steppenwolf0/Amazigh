#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "supportC.h"
double **idMatrix(int n)
        {
        	int i,j;
            double **I=newDouble(n,n);
			
            for ( i=0;i<n;i++)
            {
                for ( j = 0; j < n; j++)
                {
                    if (i==j)
                    {
                        I[i][j] = 1;
                    }
                    else
                    {
                        I[i][j] = 0;
                    }
                    
                }
                
            }

            return I;
        }

double** multMatrix(double** Matrix1, double** Matrix2, int dim_1M1,int dim_2M1,int dim_1M2,int dim_2M2)
        {
        	
            //double[,] matrixdouble = new double[Matrix1.GetLength(0), Matrix2.GetLength(1)];
            int i,j,k;
            double** matrixdouble = newDouble(dim_1M1,dim_2M2);
			

            for ( i = 0; i < dim_1M1; i++)
            {
                for ( j = 0; j < dim_2M2; j++)
                {
                    double suma = 0;
                    for ( k = 0; k < dim_2M1; k++)
                    {
                        suma += Matrix1[i][ k] * Matrix2[k ][j];
                    }
                    matrixdouble[i ][j] = suma;
                }
            }
            return matrixdouble;
        }
void swapRows(double** A, int row1, int row2,int dim_2)
        {
            if (row1 == row2) 
			return;  // Only swap if indicies are different
			int i;
            for ( i = 0; i < dim_2; i++)
            {
                double r1val = A[row1][ i];
                double r2val = A[row2 ][i];

                A[row1][ i] = r2val;
                A[row2][ i] = r1val;
            }
        }

void multiply(double** A, int row, double f,int dim_2)
        {
        	int i;
            for ( i = 0; i < dim_2; i++)
            {
                A[row][ i] *= f;
            }
        }

void zeroCoefficient(double** A, int pvt, int row,int dim_2)
        {
            // Multiply row[ r1 ] by the opposite of the coefficient at [ r2, r1 ]
            // and add it to r2
            int j;
            double multiplier = -A[row][ pvt];
            for ( j = 0; j < dim_2; j++)
            {
                A[row][j] += multiplier * A[pvt][ j];
            }
        }

void zeroCoefficient2(double** A, int pvt, int row, double** A_Inv,int dim_2)
        {
            // Multiply row[ r1 ] by the opposite of the coefficient at [ r2, r1 ]
            // and add it to r2
            int j;
            double multiplier = -A[row][ pvt];
            for (j = 0; j < dim_2; j++)
            {
                A_Inv[row][ j] += multiplier * A_Inv[pvt][ j];
            }
        }

double** inverseMatrix(double** matrix, int dim)
        {
            // matrix.GetLength(0).Show();
            int i;
            double** A = newDouble(dim,dim);
            
			double** A_inv = newDouble(dim,dim);
			
            int j;
            for ( i = 0; i < dim; i++)
            {
                for ( j = 0; j < dim; j++)
                {
                    A[i][ j] = matrix[i][ j];
                    if (i == j)
                    {
                        A_inv[i][ j] = 1;
                    }
                    else
                    {
                    	A_inv[i][ j] = 0;
					}
                }


            }

            // Use gaussian elimination to get matrix in row echelon form.
            int pvt;
            for ( pvt = 0; pvt < dim; pvt++)
            {
                // Within col, find coefficient with highest abs. from rows i to n - 1.
                // Swap with that row. Dividing by highest multiple reduces rounding errors
                int bestRow = pvt;
                int row;
                for ( row = pvt + 1; row < dim; row++)
                {
                    if (fabs(A[row][ pvt]) > fabs(A[bestRow][ pvt]))
                    {
                        bestRow = row;
                    }
                }

                swapRows(A_inv, pvt, bestRow,dim);
                swapRows(A, pvt, bestRow,dim);

                // If the coefficient is zero, then the variable probably isn't used
                if (A[pvt][ pvt] != 0)
                {
                    // Divide by the reciprocal of the pivot to get it to 1
                    double f = 1.0 / A[pvt][ pvt];
                    multiply(A, pvt, f,dim);
                    multiply(A_inv, pvt, f,dim);

                    // Zero out the rest of the coefficients in the same column, working from top to bottom
                    int k;
                    for ( k = pvt + 1; k < dim; k++)
                    {
                        zeroCoefficient2(A, pvt, k, A_inv,dim);
                        zeroCoefficient(A, pvt, k,dim);

                    }
                }
            }

            // Back substitution from right to left to get matrix in reduced row echelon form
            int k;
            for ( i = dim - 1; i >= 0; i--)
            {
                for ( k = i - 1; k >= 0; k--)
                {
                    zeroCoefficient2(A, i, k, A_inv,dim);
                    zeroCoefficient(A, i, k,dim);

                }
            }
			freeMatrix(A, dim);
            return A_inv;

        }

double** multConstant(double** Matrix1, double constant, int dim_1,int dim_2)
        {
            double** matrixdouble = newDouble(dim_1,dim_2);
			int i,j;
            for ( i = 0; i < dim_1; i++)
            {
                for ( j = 0; j < dim_2; j++)
                {
                    matrixdouble[i][ j] = Matrix1[i][ j] * constant;
                }
            }
            return matrixdouble;
        }
    
double** addMatrix(double** Matrix1, double** Matrix2, int dim_1,int dim_2 )
        {
            double**outMatrix = newDouble(dim_1,dim_2);
            int i,j;
            for ( i = 0; i < dim_1; i++)
            {
                for ( j = 0; j < dim_2; j++)
                {
                    outMatrix[i][ j] = Matrix1[i][ j] + Matrix2[i][ j];
                }
            }

            return outMatrix;
        }

double** subMatrix(double** Matrix1, double** Matrix2, int dim_1,int dim_2 )
        {
            double**outMatrix = newDouble(dim_1,dim_2);
            int i,j;
            for ( i = 0; i < dim_1; i++)
            {
                for ( j = 0; j < dim_2; j++)
                {
                    outMatrix[i][ j] = Matrix1[i][ j] - Matrix2[i][ j];
                }
            }

            return outMatrix;
        }

double** transposeMatrix(double** Matrix, int dim_1,int dim_2)
        {
            double** matrixOut = newDouble(dim_2,dim_1);

            for (int i=0;i<dim_2;i++)
            {
                for (int j=0;j<dim_1;j++)
                {
                    matrixOut[i][ j] = Matrix[j][ i];
                }
            }

            return matrixOut;
        }

double detMatrix( double** matrix, int n)
        {
            
            int i, j, j1, j2;
            double det = 0;
            double** m;

            if (n < 1)
            { /* Error */

            }
            else if (n == 1)
            { /* Shouldn't get used */
                det = matrix[0][ 0];
            }
            else if (n == 2)
            {
                det = matrix[0][ 0] * matrix[1][1] - matrix[1][ 0] * matrix[0][ 1];
            }
            else
            {
                det = 0;
                for (j1 = 0; j1 < n; j1++)
                {
                    m = newDouble(n - 1, n - 1);
                    for (i = 1; i < n; i++)
                    {
                        j2 = 0;
                        for (j = 0; j < n; j++)
                        {
                            if (j == j1)
                                continue;
                            m[i - 1][ j2] = matrix[i][ j];
                            j2++;
                        }
                    }
                    det += pow(-1.0, 1.0 + j1 + 1.0) * matrix[0][ j1] * detMatrix(m,n - 1);
                    freeMatrix(m,n-1);
                }
            }

			//NOTE FREE DET
            return (det);
        }
        
double* multVector(double** Matrix1, double* Vector, int dim_1,int dim_2)
        {
            double* vectorDouble = newDoubleV(dim_1);

            for (int i = 0; i < dim_1; i++)
            {
                double suma = 0;
                for (int j = 0; j < dim_2; j++)
                {
                    suma += Matrix1[i][ j] * Vector[j];                 
                    
                }
                vectorDouble[i] = suma;
            }
            return vectorDouble;
        }
        
double** cutMatrix(double** Matrix, int x0, int y0, int xn, int yn, int dim_1,int dim_2 )
        {
            double** MatrixOut = newDouble(yn - y0, xn - x0);

            for (int i = 0; i < dim_1; i++)
            {

                for (int j = 0; j < dim_2; j++)
                {
                    if (i >= y0 && i < yn)
                    {
                        if (j >= x0 && j < xn)
                        {
                            MatrixOut[i - y0][ j - x0] = Matrix[i][ j];
                        }
                    }
                }

            }
            return MatrixOut;
        }

double* getVector(double** Matrix,int index, int dim)
{
	double* vector=newDoubleV(dim);
	for (int i=0;i<dim;i++)
	{
		vector[i]=Matrix[i][index];
	}
	return vector;
}

double frobeniusNorm(double** Matrix,int dim_1,int dim_2)
{
	double outVal = 0;
	int i, j;
	for (i = 0; i < dim_1; i++)
	{
		for (j = 0; j < dim_2; j++)
		{
			outVal += Matrix[i][j] * Matrix[i][j];
		}
	}

	return sqrt(outVal);
}

double matrixNorm(double** Matrix, int dim_1, int dim_2)
{
	double best = -1e20;
	
	int i, j;
	for (i = 0; i < dim_1; i++)
	{
		double outVal = 0;
		for (j = 0; j < dim_2; j++)
		{
			outVal += fabs(Matrix[j][i]);
		}
		if (outVal > best)
			best = outVal;
	}

	return best;
}

double traceMatrix(double** Matrix, int size)
{
	int i = 0;
	double result = 0;
	for (i = 0; i < size; i++)
	{
		result += Matrix[i][i];
	}
	return result;
}

double** Cholesky(double** matrix, int size)
{
	double** L = newDouble(size, size);
	L[0][0] = sqrt(matrix[0][0]);

	int i, j, k;

	/*for ( i = 0; i<size; i++)
	{
	for ( j = 0; j<size; j++)
	{
	L[i][j] = 0;
	}
	}*/

	double sumk = 0;
	double sumj = 0;

	for (i = 0; i<size; i++)
	{
		sumk = 0;
		for (k = 0; k< i; k++)
		{
			sumk += pow(L[i][k], 2.0);
		}

		L[i][i] = sqrt(matrix[i][i] - sumk);

		for (j = i + 1; j<size; j++)
		{
			sumj = 0;
			for (k = 0; k< i; k++)
			{
				sumj += L[i][k] * L[j][k];
			}
			L[j][i] = (matrix[j][i] - sumj) / L[i][i];
		}
	}

	return L;
}

double logDet(double **a, int n)
{
	//       log(det(A)) = 2*sum(log(vecdiag(G))) 
//https://blogs.sas.com/content/iml/2012/10/31/compute-the-log-determinant-of-a-matrix.html
	double** chol = Cholesky(a, n);
	double result = 0;
	for (int i = 0; i < n; i++)
	{
		result += log(chol[i][i]);
	}
	freeMatrix(chol, n);
	result = result * 2;
	return(result);
}

double det4(double** m) {
	return
		m[0][3] * m[1][2] * m[2][1] * m[3][0] - m[0][2] * m[1][3] * m[2][1] * m[3][0] -
		m[0][3] * m[1][1] * m[2][2] * m[3][0] + m[0][1] * m[1][3] * m[2][2] * m[3][0] +
		m[0][2] * m[1][1] * m[2][3] * m[3][0] - m[0][1] * m[1][2] * m[2][3] * m[3][0] -
		m[0][3] * m[1][2] * m[2][0] * m[3][1] + m[0][2] * m[1][3] * m[2][0] * m[3][1] +
		m[0][3] * m[1][0] * m[2][2] * m[3][1] - m[0][0] * m[1][3] * m[2][2] * m[3][1] -
		m[0][2] * m[1][0] * m[2][3] * m[3][1] + m[0][0] * m[1][2] * m[2][3] * m[3][1] +
		m[0][3] * m[1][1] * m[2][0] * m[3][2] - m[0][1] * m[1][3] * m[2][0] * m[3][2] -
		m[0][3] * m[1][0] * m[2][1] * m[3][2] + m[0][0] * m[1][3] * m[2][1] * m[3][2] +
		m[0][1] * m[1][0] * m[2][3] * m[3][2] - m[0][0] * m[1][1] * m[2][3] * m[3][2] -
		m[0][2] * m[1][1] * m[2][0] * m[3][3] + m[0][1] * m[1][2] * m[2][0] * m[3][3] +
		m[0][2] * m[1][0] * m[2][1] * m[3][3] - m[0][0] * m[1][2] * m[2][1] * m[3][3] -
		m[0][1] * m[1][0] * m[2][2] * m[3][3] + m[0][0] * m[1][1] * m[2][2] * m[3][3];
}