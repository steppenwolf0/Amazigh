#include <math.h>
#include "supportC.h"
double* idVector(int n, double constant)
        {
            double* I = newDoubleV( n);

            for (int i=0;i<n;i++)
            {
                I[i] = constant;
            }

            return I;
        }



double dotVector(double* vector1, double* vector2, int dim)
        {
            double suma = 0;
            double vectordotproductout;
            for (int i = 0; i < dim; i++)
            {
                suma += vector1[i] * vector2[i];

            }
            vectordotproductout = suma;
            return vectordotproductout;
        }
        
double* multConstantV(double* vector, double constant, int dim)
        {
            double* vectorOutput =  newDoubleV(dim);
            for (int i = 0; i < dim; i++)
            {
                vectorOutput[i] = constant * vector[i];
            }
            return vectorOutput;
        }

double* cross3(double* a, double* b)
        {
            double* output = newDoubleV(3);
            output[0] = a[1] * b[2] - a[2] * b[1];
            output[1] = a[2] * b[0] - a[0] * b[2];
            output[2] = a[0] * b[1] - a[1] * b[0];
            return output;
        }

double norm(double* vector, int dim)
        {
            double suma = 0;
            for (int i = 0; i < dim; i++)
            {
                suma += vector[i] * vector[i];
            }
            return sqrt(suma);
        }

double* addVector(double* Vector1, double* Vector2, int dim)
        {
            double* outVector = newDoubleV(dim);
            for (int i = 0; i < dim; i++)
            {

                outVector[i] = Vector1[i] + Vector2[i];

            }

            return outVector;
        }
        
double* subVector(double* Vector1, double* Vector2, int dim)
        {
            double* outVector = newDoubleV(dim);
            for (int i = 0; i < dim; i++)
            {

                outVector[i] = Vector1[i] - Vector2[i];
                
            }

            return outVector;
        }
        
double** transVector(double* vector, int dim)
        {
            double** vT = newDouble(1, dim);

            for (int i=0;i<dim;i++)
            {
                vT[0][ i] = vector[i];
            }
            return vT;
        }

double maxVector(double* Vector, int dim)
{
	double maxValue = Vector[0];
	for (int i = 0; i < dim; i++)
	{

		if (Vector[i] > maxValue)
			maxValue = Vector[i];
		

	}

	return maxValue;
}
        

