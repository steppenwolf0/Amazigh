#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mathMatrix.h"
#include "mathVector.h"
#include "printFile.h"
#include "supportC.h"

double* PBCGSTABC(double** A, double* b, int size)
{
	int iterations = 3000;
	int i;
	double** M = Cholesky(A, size);

	double* x = newDoubleV(size);
	//r=b-A*x;
	double* Ax = multVector(A, x, size, size);
	double* r = subVector(b, Ax, size);
	//rbis=r;
	double* rbis = newDoubleV(size);
	for (i = 0; i < size; i++)
	{
		rbis[i] = r[i];
	}


	double ro = 1;
	double roant = 1;
	double w = 1;
	double alpha = 1;

	double* v = newDoubleV(size);
	double* p = newDoubleV(size);

	double res = dotVector(r, r, size);

	int its = 0;
	while (its < iterations && res > 1e-20)
	{
		//?i = (r?0, ri?1)
		ro = dotVector(rbis, r, size);
		//? = (?i/?i?1)(?/?i?1)
		double beta = (ro / roant) * (alpha / w);
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//temp1=?i?1vi?1
		//temp2=pi?1 ? temp1
		//temp3= ?(temp2)
		double* temp1 = multConstantV(v, w, size);
		double* temp2 = subVector(p, temp1, size);
		double* temp3 = multConstantV(temp2, beta, size);
		//pi = ri?1 +temp3
		free(p);
		p = addVector(r, temp3, size);
		//y = K^(?1)*pi
		double* y = multVector(M, p, size, size);
		//vi = Ay
		free(v);
		v = multVector(A, y, size, size);
		//? = ?i/(r?0, vi)

		alpha = ro / dotVector(rbis, v, size);
		//s = ri?1 ? ?vi
		//temp1=?vi
		//s = ri?1 ? temp1
		double* temp4 = multConstantV(v, alpha, size);
		double* s = subVector(r, temp4, size);
		//z = K^(?1)*s
		double* z = multVector(M, s, size, size);
		//t = Az
		double* t = multVector(A, z, size, size);
		//z2=K^(?1)_1*t
		double* z2 = multVector(M, t, size, size);
		//?i = (z2, K^(?1)_1s)/(z2, z2)
		w = dotVector(z2, z, size) / dotVector(z2, z2, size);
		//xi = xi?1 + ?y + ?iz
		double* temp5 = multConstantV(z, w, size);
		double* temp6 = multConstantV(y, alpha, size);
		double* temp7 = addVector(temp5, temp6, size);
		double* xtemp = addVector(x, temp7, size);
		for (i = 0; i < size; i++)
		{
			x[i] = xtemp[i];
		}

		//ri = s ? ?it
		double* temp8 = multConstantV(t, w, size);
		free(r);
		r = subVector(s, temp8, size);
		roant = ro;
		res = dotVector(s, s, size);
		its++;

		free(s);
		free(z);
		free(t);
		free(z2);
		free(y);
		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
		free(temp5);
		free(temp6);
		free(temp7);
		free(temp8);
		free(xtemp);
	}


	free(v);
	free(p);
	free(rbis);
	free(r);
	free(Ax);

	for (i = 0; i < size; i++)
	{
		free(M[i]);
	}
	free(M);


	return x;

}

double* multVectorSparse(double* Avalues, int** Apositions, int Asize, double* x, int sizeVector)
{
	double* out = newDoubleV(sizeVector);
	for (int i = 0; i < Asize; i++)
	{
		out[Apositions[i][0]] += x[Apositions[i][1]] * Avalues[i];
	}

	return out;
}

double* multTransposeVector(double* M, double* x, int size)
{
	double* out = newDoubleV(size);
	for (int i = 0; i < size; i++)
	{
		out[i] = M[i] * x[i];
	}
	return out;
}

double* PBCGSTABCSparse(double* Avalues, int** Apositions, int Asize, double* b, int size)
{
	int iterations = 3000;

	double* M = newDoubleV(size);


	for (int i = 0; i < Asize; i++)
	{
		if (Apositions[i][0] == Apositions[i][1])
		{
			M[Apositions[i][0]] = 1.0 / Avalues[i];
		}
	}

	//*M[i,i]=A[i,i];
	//Minv[i,i]=1/A[i,i];
	//M=Minv;

	double* x = newDoubleV(size);
	//r=b-A*x;
	double* Ax = multVectorSparse(Avalues, Apositions, Asize, x, size);
	double* r = subVector(b, Ax, size);
	//rbis=r;
	double* rbis = newDoubleV(size);
	for (int i = 0; i < size; i++)
	{
		rbis[i] = r[i];
	}


	double ro = 1;
	double roant = 1;
	double w = 1;
	double alpha = 1;

	double* v = newDoubleV(size);
	double* p = newDoubleV(size);

	double res = dotVector(r, r, size);

	int its = 0;
	while (its < iterations && res > 1e-20)
	{
		
		//?i = (r?0, ri?1)
		ro = dotVector(rbis, r, size);
		//? = (?i/?i?1)(?/?i?1)
		double beta = (ro / roant) * (alpha / w);
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//temp1=?i?1vi?1
		//temp2=pi?1 ? temp1
		//temp3= ?(temp2)
		double* temp1 = multConstantV(v, w, size);
		double* temp2 = subVector(p, temp1, size);
		double* temp3 = multConstantV(temp2, beta, size);
		//pi = ri?1 +temp3
		free(p);
		p = addVector(r, temp3, size);
		//y = K^(?1)*pi
		double* y = multTransposeVector(M, p, size);
		//vi = Ay
		free(v);
		v = multVectorSparse(Avalues, Apositions, Asize, y, size);
		//? = ?i/(r?0, vi)

		alpha = ro / dotVector(rbis, v, size);
		//s = ri?1 ? ?vi
		//temp1=?vi
		//s = ri?1 ? temp1
		double* temp4 = multConstantV(v, alpha, size);
		double* s = subVector(r, temp4, size);
		//z = K^(?1)*s
		double* z = multTransposeVector(M, s, size);
		//t = Az
		double* t = multVectorSparse(Avalues, Apositions, Asize, z, size);
		//z2=K^(?1)_1*t
		double* z2 = multTransposeVector(M, t, size);
		//?i = (z2, K^(?1)_1s)/(z2, z2)
		w = dotVector(z2, z, size) / dotVector(z2, z2, size);
		//xi = xi?1 + ?y + ?iz
		double* temp5 = multConstantV(z, w, size);
		double* temp6 = multConstantV(y, alpha, size);
		double* temp7 = addVector(temp5, temp6, size);
		double* xtemp = addVector(x, temp7, size);
		for (int i = 0; i < size; i++)
		{
			x[i] = xtemp[i];
		}

		//ri = s ? ?it
		double* temp8 = multConstantV(t, w, size);
		free(r);
		r = subVector(s, temp8, size);
		roant = ro;
		res = dotVector(s, s, size);
		its++;

		free(s);
		free(z);
		free(t);
		free(z2);
		free(y);
		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
		free(temp5);
		free(temp6);
		free(temp7);
		free(temp8);
		free(xtemp);
	}
	printf("iterations %d \n", its);

	free(v);
	free(p);
	free(rbis);
	free(r);
	free(Ax);

	
	free(M);


	return x;

}


void eulerODE(double* RATES, double* STATES, int sizeSTATES, double dt)
{

	for (int j = 0; j < sizeSTATES; j++)
	{
		STATES[j] = RATES[j] * dt + STATES[j];
	}

}


//public static double[] SolveSparsePCG(int[, ] MIJ, double[] MV, int[, ] IJ,
//	double[] V, int length, double[] B, double[] x0, double ptol)
//{
//
//
//
//	double[] b = B.Minus();
//
//	
//
//
//	int its = 0;
//	int MaxItn = 3000;
//	double[] f = x0; //initial guess
//	double[] g = SparseMatrixMult(IJ, V, f, length).Add(b);
//	
//	double[] z = SparseMatrixMult(MIJ, MV, g, length);
//	double[] p = z.Minus();
//	double s = g.Dot(z);
//	double[] h;
//	while (its < MaxItn && s > ptol)
//	{
//		s = g.Dot(z);
//		h = SparseMatrixMult(IJ, V, p, length);
//		double t = s / (p.Dot(h));
//		f = f.Add(p.Mult(t));
//		g = g.Add(h.Mult(t));
//		
//		z = SparseMatrixMult(MIJ, MV, g, length);
//		double s1 = g.Dot(z);
//		double beta = s1 / s;
//		p = z.Minus().Add(p.Mult(beta));
//		its++;
//
//	}
//
//	
//	Console.WriteLine(its.ToString());
//	return f;
//}

double* SolveSparsePCG( int** IJ, 
	double* AV, int length, double* B, double* x0,int size, double ptol)
{

	double* MV = newDoubleV(length);


	for (int i = 0; i < length; i++)
	{
		if (IJ[i][0] == IJ[i][1])
		{
			MV[i] = 1.0 / AV[i];
		}
	}

	double* b = multConstantV(B, -1, size);

	int its = 0;
	int MaxItn = 3000;
	double* f = newDoubleV(size);
	copyVector(f, x0, size);
	double* gTemp = multVectorSparse(AV, IJ, length, f, size);
	double* g = addVector(gTemp, b, size);

	double* z = multVectorSparse(MV, IJ, length, g, size);
	double* p = multConstantV(z, -1, size);
	double s = dotVector(g, z, size);
	
	

	while (its < MaxItn && s > ptol)
	{
		s = dotVector(g, z, size);
		
		double* h = multVectorSparse(AV, IJ,  length, p, size);
		double t = s / (dotVector(p, h, size));
		
		double* pt = multConstantV(p, t, size);
		double* tempf = addVector(f, pt,size);
		copyVector(f, tempf, size);

		double* ht = multConstantV(h, t, size);
		double* tempg = addVector(g, ht, size);
		copyVector(g, tempg, size);

		double* tempz = multVectorSparse( MV, IJ, length, g, size);
		copyVector(z, tempz, size);

		double s1 = dotVector(g, z, size);
		double beta = s1 / s;

		double* pbeta = multConstantV(p, beta, size);
		double* zMinus = multConstantV(z, -1, size);
		double* tempp = addVector(zMinus, pbeta, size);
		copyVector(p, tempp, size);
		
		free(ht);
		free(tempg);
		free(tempz);
		free(pbeta);
		free(zMinus);
		free(tempp);
		free(h);
		free(pt);
		free(tempf);
		its++;
	}


	printf("its %d \n", its);
	//end
	
	free(gTemp);
	free(g);
	free(z);
	free(p);
	free(b);
	free(MV);
	return f;
}