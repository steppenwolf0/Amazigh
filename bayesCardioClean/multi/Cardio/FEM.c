#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Amazigh/supportC.h"
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/mathVector.h"
#include <math.h>
#include "FEMTests.h"

//*********************************************************Regular FEM**************************************

double massFunc(double x, double y, double z)
{
	double value=1;
	return value;
}

double** Mass3DforFEM(double** Nodes, int** Triangles, int** Tetrahedra,
	int** NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize, 
	int constrainedsurface)
{
	//Supposing that the surfaces does not intersect
	double** K = newDouble(nodeSize, nodeSize);

	for (int i = 0; i < tetrahedraSize; i++)
	{
		//Generate the element matrix for each element
		//[1 x1 y1 z1]
		//[1 x2 y2 z2]
		//[1 x3 y3 z3]
		//[1 x4 y4 z4]
		double** M = newDouble(4, 4);

		for (int l = 0; l < 4; l++)
		{

			for (int m = 0; m < 4; m++)
			{
				if (l == 0)
				{
					M[m][l] = 1;
				}
				else
				{
					M[m][ l] = Nodes[Tetrahedra[i][ m]][ l - 1];
				}

			}

		}

		double** G = newDouble(4, 4);

		//ERROR***********
		for (int l = 0; l < 4; l++)
		{
			for (int m = 0; m < 4; m++)
			{
				if (m == l)
				{
					G[m][l] = 3.0 / 9.0;

				}
				else
				{
					G[m][l] = 2.0 / 9.0;
				}
			}
		}
		double** Gtemp = multConstant(G, 0.25,4,4);
		copyMatrix(G, Gtemp, 4, 4);
		freeMatrix(Gtemp, 4);


		//ERROR***********
		double detM = det4(M);
		double V=fabs(detM) / 6.0;

		// Add variables and values to hashtable
		//*********************************this is used to have the same as the example erase later

		double I;
		

			double sumc1 = M[1][1] + M[2][1] + M[3][1] + M[0][1];
			double sumc2 = M[1][2] + M[2][2] + M[3][2] + M[0][2];
			double sumc3 = M[1][3] + M[2][3] + M[3][3] + M[0][3];

			double x = sumc1 / 4;
			double y = sumc2 / 4;
			double z = sumc3 / 4;

			double result = massFunc(x, y, z);
	
			I = V * result;

		


		//The triangle contributes to at most 6 entries in the (upper triangle
		//of the) stiffness matrix.  We compute these six entries in the
		//following double loop.
		for (int s = 0; s < 4; s++)
		{
			int lls = Tetrahedra[i][s];
			if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
			{
				for (int r = 0; r <= s; r++)
				{
					//If both vertices are free, then there is a contribution
					//to the stiffness matrix
					int llr = Tetrahedra[i][r];
					if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
					{
						if (llr <= lls)
						{
							K[llr][lls] = K[llr][lls] + G[r][s] * I;
						}
						else
						{
							K[lls][llr] = K[lls][llr] + G[r][s] * I;
						}
					}
				}
			}
		}
		//end of principal for    
		freeMatrix(G, 4);
		freeMatrix(M, 4);
	}

	for (int l = 0; l < nodeSize; l++)
	{
		for (int m = 0; m < nodeSize; m++)
		{
			if (K[l][m] != 0)
			{
				K[m][l] = K[l][m];

			}
			else
			{
				if (m == l)
				{
					K[m][l] = 1;
				}
			}
		}

	}

	return K;
}

double** Stiffness3DforFEM(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface)
{
	//Supposing that the surfaces does not intersect
	double** K = newDouble(nodeSize, nodeSize);
	
	for (int i = 0; i < tetrahedraSize; i++)
	{
		//Generate the element matrix for each element
		//[1 x1 y1 z1]
		//[1 x2 y2 z2]
		//[1 x3 y3 z3]
		//[1 x4 y4 z4]
		double** M = newDouble(4, 4);

		for (int l = 0; l < 4; l++)
		{

			for (int m = 0; m < 4; m++)
			{
				if (l == 0)
				{
					M[m][ l] = 1;
				}
				else
				{
					M[m][l] = Nodes[Tetrahedra[i][m]][l - 1];
				}

			}

		}

		//Get the coeficientes for each basis function
		//[a1 a2 a3 a4]
		//[b1 b2 b3 b4]
		//[c1 c2 c3 c4]
		//[d1 d2 d3 d4]
		double** C = inverseMatrix(M, 4);
			
		//The typical integral we must compute is
		//integral over Ti (k(x,y)(grad Ø1).(grad Ø2).
		//The gradients are constant vectors, since the basis functions
		//are linear over this triangle.   So we just compute all the
		//dot products and store them in a matrix G:
		

		double** C34 = newDouble(3, 4);

		for (int l = 1; l < 4; l++)
		{
			for (int m = 0; m < 4; m++)
			{
				C34[l - 1][m] = C[l][m];

			}
		}


		//multiply gradient(Øi)*gradient(Øj)
		double** C34T = transposeMatrix(C34, 3, 4);
		double** G = multMatrix(C34T, C34, 4, 3, 3, 4);
		//G = C34.T().Mult(C34);
		double detM = det4(M);
		
		double V = fabs(detM) / 6.0;
		//double V = Math.Abs(M.Determinant()) / 6.0;
		
		// Add variables and values to hashtable
		//*********************************this is used to have the same as the example erase later

			double sumc1 = M[1][1] + M[2][1] + M[3][1] + M[0][1];
			double sumc2 = M[1][2] + M[2][2] + M[3][2] + M[0][2];
			double sumc3 = M[1][3] + M[2][3] + M[3][3] + M[0][3];

			double x = sumc1 / 4;
			double y = sumc2 / 4;
			double z = sumc3 / 4;

			// Parse and write the result
			double result = stiffFunc(x, y, z);

			double	I = V * result;

		

		//The triangle contributes to at most 6 entries in the (upper triangle
		//of the) stiffness matrix.  We compute these six entries in the
		//following double loop.
		for (int s = 0; s < 4; s++)
		{
			int lls = Tetrahedra[i][s];
			if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
			{
				for (int r = 0; r <= s; r++)
				{
					//If both vertices are free, then there is a contribution
					//to the stiffness matrix
					int llr = Tetrahedra[i][r];
					if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
					{
						if (llr <= lls)
						{
							K[llr][lls] = K[llr][lls] + G[r][s] * I;
						}
						else
						{
							K[lls][llr] = K[lls][llr] + G[r][s] * I;
						}
					}
				}
			}
		}
		//end of principal for    
		freeMatrix(M, 4);
		freeMatrix(G, 4);
		freeMatrix(C, 4);
		freeMatrix(C34, 3);
		freeMatrix(C34T, 4);
	}

	for (int l = 0; l < nodeSize; l++)
	{
		for (int m = 0; m < nodeSize; m++)
		{
			if (K[l][m] != 0)
			{
				K[m][l] = K[l][m];

			}
			else
			{
				if (m == l)
				{
					K[m][l] = 1;
				}
			}
		}



	}


	return K;


}

double** GetNodesNormals(double** Nodes, int** Triangles, int nodeSize, int triangleSize)
{
	double** NodesNormals = newDouble(nodeSize, 3);

	double** TriangleNormals = newDouble(triangleSize, 3);


		for (int i = 0; i < triangleSize; i++)
		{

			double* temp1 = newDoubleV(3);
			double* temp2 = newDoubleV(3);
			double* X = newDoubleV(3);
			double* Y = newDoubleV(3);
			double* Z = newDoubleV(3);

			for (int j = 0; j < 3; j++)
			{

				X[j] = Nodes[Triangles[i][ j]][0];
				Y[j] = Nodes[Triangles[i][j]][1];
				Z[j] = Nodes[Triangles[i][j]][2];
				temp1[j] = (Nodes[Triangles[i][0]][j] - Nodes[Triangles[i][1]][j]);
				temp2[j] = (Nodes[Triangles[i][0]][j] - Nodes[Triangles[i][2]][j]);

			}
			
			double* nv = cross3(temp1, temp2);
			double n = norm(nv, 3);


			for (int j = 0; j < 3; j++)
			{
				nv[j] = nv[j] / n;
			}



			TriangleNormals[i][0] = nv[0];
			TriangleNormals[i][1] = nv[1];
			TriangleNormals[i][2] = nv[2];
			
			free(temp1);
			free(temp2);
			free(X);
			free(Y);
			free(Z);
			free(nv);
		}

		int* NodesAppearance = newIntV(nodeSize);

		for (int i = 0; i < triangleSize; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				NodesAppearance[Triangles[i][j]]++;
			}
		}

		for (int i = 0; i < triangleSize; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				NodesNormals[Triangles[i][j]][0] += TriangleNormals[i][0];
				NodesNormals[Triangles[i][j]][1] += TriangleNormals[i][1];
				NodesNormals[Triangles[i][j]][2] += TriangleNormals[i][2];
			}
		}

		for (int i = 0; i < nodeSize; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				NodesNormals[i][ j] = NodesNormals[i][j] / NodesAppearance[i];
			}
		}

		freeMatrix(TriangleNormals, triangleSize);
		free(NodesAppearance);
	return NodesNormals;
}


int neumanTriangles(int* triangleSurface,int triangleSize, int neumannSurface)
{
	int nTriangles = 0;
	for (int i = 0; i < triangleSize; i++)
	{
		if (triangleSurface[i] == neumannSurface)
		{
			nTriangles++;
		}

	}
	return nTriangles;
}

void getNeumanData3DNEW (double** Nodes, int** Triangles, int nodeSize, int triangleSize,
	 int* triangleSurface, int neumannSurface, double** NeumanData, int** Neuman, 
	int ntriangles, double t)
{

	
	int temp = 0;
	for (int i = 0; i < triangleSize; i++)
	{
		if (triangleSurface[i] == neumannSurface)
		{
			Neuman[temp][ 0] = Triangles[i][0];
			Neuman[temp][1] = Triangles[i][1];
			Neuman[temp][2] = Triangles[i][2];
			temp++;
		}

	}

	
	
	if (ntriangles>0)
	{
		

		double* g0 = newDoubleV(3);
		double* g1 = newDoubleV(3);
		double* g2 = newDoubleV(3);

		double* n0 = newDoubleV(3);
		double* n1 = newDoubleV(3);
		double* n2 = newDoubleV(3);

		double* x = newDoubleV(3);
		double* y = newDoubleV(3);
		double* z = newDoubleV(3);

		double* k = newDoubleV(3);


		double** NodeNormals = GetNodesNormals(Nodes, Triangles, nodeSize,triangleSize);

		for (int i = 0; i < ntriangles; i++)
		{


			for (int j = 0; j < 3; j++)
			{
				n0[j] = NodeNormals[Neuman[i][ 0]][j];
				n1[j] = NodeNormals[Neuman[i][1]][j];
				n2[j] = NodeNormals[Neuman[i][2]][j];
			}




			for (int j = 0; j < 3; j++)
			{

				x[j] = Nodes[Neuman[i][j]][0];
				y[j] = Nodes[Neuman[i][j]][1];
				z[j] = Nodes[Neuman[i][j]][2];



			}

			double result = 0;
			result = ux(x[0], y[0], z[0], t);
			g0[0] = result;
			result = uy(x[0], y[0], z[0], t);
			g0[1] = result;
			result = uz(x[0], y[0], z[0], t);
			g0[2] = result;

			result = ux(x[1], y[1], z[1], t);
			g1[0] = result;
			result = uy(x[1], y[1], z[1], t);
			g1[1] = result;
			result = uz(x[1], y[1], z[1], t);
			g1[2] = result;

			result = ux(x[2], y[2], z[2], t);
			g2[0] = result;
			result = uy(x[2], y[2], z[2], t);
			g2[1] = result;
			result = uz(x[2], y[2], z[2], t);
			g2[2] = result;
			
			//k0
			result = stiffFunc(x[0], y[0], z[0]);
			k[0] = result;
			//k1
			result = stiffFunc(x[1], y[1], z[1]);
			k[1] = result;
			//k2
			result = stiffFunc(x[2], y[2], z[2]);
			k[2] = result;

			

			NeumanData[i][0] = dotVector(g0, n0, 3) * k[0];
			NeumanData[i][1] = dotVector(g1, n1, 3) * k[1];
			NeumanData[i][2] = dotVector(g2, n2, 3) * k[2];


		}

		free(g0);
		free(g1);
		free(g2);
		free(n0);
		free(n1);
		free(n2);
		free(x);
		free(y);
		free(z);
		free(k);
		freeMatrix(NodeNormals, nodeSize);

	}


	

}

double* getDirichletData(double** Nodes, int nodeSize, int constrainedSurface, int* nodeSurface, double t)
{
	double* data = newDoubleV(nodeSize);
	for (int i = 0; i < nodeSize; i++)
	{
		if (nodeSurface[i] == constrainedSurface)
		{
			double x = Nodes[i][0];
			double y = Nodes[i][1];
			double z = Nodes[i][2];

			data[i] = uFunc(x, y, z, t);
		}
	}
	return data;
}

double** Load3D(double** Nodes, int** Tetrahedra, int* NodeSurface, int nodeSize, 
	int tetrahedraSize, int constrainedsurface, double* nodeDirichlet, double** NeumanData, 
	int** Neuman, int ntriangles, double t)
{
	double* F = newDoubleV(nodeSize);

	//Add the contribution from each element

	for (int i = 0; i < tetrahedraSize; i++)
	{
		//Generate the element matrix for each element
		//[1 x1 y1 z1]
		//[1 x2 y2 z2]
		//[1 x3 y3 z3]
		//[1 x4 y4 z4]
		double** M = newDouble(4, 4);

		for (int l = 0; l < 4; l++)
		{

			for (int m = 0; m < 4; m++)
			{
				if (l == 0)
				{
					M[m][ l] = 1;
				}
				else
				{
					M[m][l] = Nodes[Tetrahedra[i][m]][l - 1];
				}

			}

		}

		//Get the coeficientes for each basis function
		//[a1 a2 a3 a4]
		//[b1 b2 b3 b4]
		//[c1 c2 c3 c4]
		//[d1 d2 d3 d4]
		double** C = inverseMatrix(M, 4);

		//The typical integral we must compute is
		//integral over Ti (k(x,y)(grad Ø1).(grad Ø2).
		//The gradients are constant vectors, since the basis functions
		//are linear over this triangle.   So we just compute all the
		//dot products and store them in a matrix G:


		double** C34 = newDouble(3, 4);

		for (int l = 1; l < 4; l++)
		{
			for (int m = 0; m < 4; m++)
			{
				C34[l - 1][m] = C[l][m];

			}
		}


		//multiply gradient(Øi)*gradient(Øj)
		double** C34T = transposeMatrix(C34, 3, 4);
		double** G = multMatrix(C34T, C34, 4, 3, 3, 4);
		//G = C34.T().Mult(C34);
		double detM = det4(M);

		double V = fabs(detM) / 6.0;
		//Apply the quadrature rule:

		
		double I = 0;
		double x = 0;
		double y = 0;
		double z = 0;
		double result = 0;

		double sumc1 = M[1][1] + M[2][1] + M[3][1] + M[0][1];
		double sumc2 = M[1][2] + M[2][2] + M[3][2] + M[0][2];
		double sumc3 = M[1][3] + M[2][3] + M[3][3] + M[0][3];

		x = sumc1 / 4.0;
		y = sumc2 / 4.0;
		z = sumc3 / 4.0;

		
		result = fFunc(x, y, z, t);
		
			//I = (V * result) / 4.0;
			I = (V * result) / 4.0;

			for (int r = 0; r < 4; r++)
			{
				int llr = Tetrahedra[i][r];
				if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
				{
					F[llr] = F[llr] + I;
				}
			}
		
		//*******************************************************************************************
		//Add the contribution from the Dirichlet data to the right hand side
		//(if the Dirichlet condition is not homogeneous).


			//transpose the matrix 



			int any = 0;
			for (int m = 0; m < 4; m++)
			{
				if (NodeSurface[Tetrahedra[i][ m]] == constrainedsurface)
				{
					any++;
				}
			}
			if (any >0)
			{

				// Get the nodal values of the piecewise linear function G
				//agreeing with the Dirichlet data at the constrained nodes
				//and having value zero at each free node.

				//w=zeros(3,1);
				double* w = newDoubleV(4);//ç

				for (int j = 0; j < 4; j++)
				{
					if (NodeSurface[Tetrahedra[i][j]] == constrainedsurface)
					{
						w[j] = nodeDirichlet[Tetrahedra[i][j]];
					}
					else
					{
						w[j] = 0;
					}
				}

				//% The inverse of M contains the coefficients of the basis
				//% functions on this triangle:

				double* C34W = multVector(C34, w, 3, 4);
				double* GV = multVector(C34T, C34W, 4, 3);
				//double[] GV = (C34.T().Mult(C34.Mult(w)));



				double* I4 = newDoubleV(4);
				for (int m = 0; m < 4; m++)
				{					
					result = stiffFunc(x, y, z);
					I4[m] = (result * GV[m] * V);
				}

				//% Now add the contributions to the right hand side
				for (int r = 0; r < 4; r++)
				{
					int llr = Tetrahedra[i][r];
					if (NodeSurface[Tetrahedra[i][ r]] != constrainedsurface)
					{
						F[llr] = F[llr] - I4[r];
					}
				}
				//if any==true
				free(w);
				free(C34W);
				free(GV);
				free(I4);
			}

			freeMatrix(M, 4);
			freeMatrix(C, 4);
			freeMatrix(C34, 3);
			freeMatrix(C34T, 4);
			freeMatrix(G, 4);
		//FIN DE DIRICHLET
	}


	if (ntriangles>0)
	{
		for (int i = 0; i < ntriangles; i++)
		{

			double* temp1 = newDoubleV(3);
			double* temp2 = newDoubleV(3);
			for (int j = 0; j < 3; j++)
			{

				temp1[j] = (Nodes[Neuman[i][0]][j] - Nodes[Neuman[i][1]][j]);
				temp2[j] = (Nodes[Neuman[i][0]][j] - Nodes[Neuman[i][2]][j]);

			}
			double* nv = cross3(temp1, temp2);
			
			double n = norm(nv, 3);

			for (int j = 0; j < 3; j++)
			{
				nv[j] = nv[j] / n;
			}

			double A = n / 2.0;
			//double len = FEM3D.Norm(x1, x2, y1, y2);
			//Evaluate the integral (use the midpoint rule):
			//double I3 = (hval*len) / 3.0;
			//double I3 = (hval * A) / 3.0;

			//double I = (NeumanData[i, 0] + NeumanData[i, 1] + NeumanData[i, 2])*A / 9;

			double I3a = (NeumanData[i][0] * A) / 3.0;
			double I3b = (NeumanData[i][1] * A) / 3.0;
			double I3c = (NeumanData[i][2] * A) / 3.0;

			F[Neuman[i][0]] = F[Neuman[i][0]] + I3a;
			F[Neuman[i][1]] = F[Neuman[i][1]] + I3b;
			F[Neuman[i][2]] = F[Neuman[i][2]] + I3c;


			free(temp1);
			free(temp2);
			free(nv);
		}
	}
	
		for (int i = 0; i < nodeSize; i++)
		{
			if (NodeSurface[i] == constrainedsurface)
			{
				F[i] = nodeDirichlet[i];
			}

		}
	



	return F;
}

//*********************************************************Regular FEM**************************************
//*********************************************************Sparse FEM***************************************


//int getPositionsOriginal(double** Nodes, int** Triangles, int** Tetrahedra,
//	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
//	int constrainedsurface, int*** outPositions)
//{
//	int size = 0;
//	int** Positions;
//	for (int i = 0; i < tetrahedraSize; i++)
//	{
//		//Generate the element matrix for each element
//		//[1 x1 y1 z1]
//		//[1 x2 y2 z2]
//		//[1 x3 y3 z3]
//		//[1 x4 y4 z4]
//		double** M = newDouble(4, 4);
//
//		for (int l = 0; l < 4; l++)
//		{
//
//			for (int m = 0; m < 4; m++)
//			{
//				if (l == 0)
//				{
//					M[m][l] = 1;
//				}
//				else
//				{
//					M[m][l] = Nodes[Tetrahedra[i][m]][l - 1];
//				}
//
//			}
//
//		}
//
//		//Get the coeficientes for each basis function
//		//[a1 a2 a3 a4]
//		//[b1 b2 b3 b4]
//		//[c1 c2 c3 c4]
//		//[d1 d2 d3 d4]
//		double** C = inverseMatrix(M, 4);
//
//		//The typical integral we must compute is
//		//integral over Ti (k(x,y)(grad Ø1).(grad Ø2).
//		//The gradients are constant vectors, since the basis functions
//		//are linear over this triangle.   So we just compute all the
//		//dot products and store them in a matrix G:
//
//
//		double** C34 = newDouble(3, 4);
//
//		for (int l = 1; l < 4; l++)
//		{
//			for (int m = 0; m < 4; m++)
//			{
//				C34[l - 1][m] = C[l][m];
//
//			}
//		}
//
//
//		//multiply gradient(Øi)*gradient(Øj)
//		double** C34T = transposeMatrix(C34, 3, 4);
//		double** G = multMatrix(C34T, C34, 4, 3, 3, 4);
//		//G = C34.T().Mult(C34);
//		double detM = det4(M);
//
//		double V = fabs(detM) / 6.0;
//		//double V = Math.Abs(M.Determinant()) / 6.0;
//
//		// Add variables and values to hashtable
//		//*********************************this is used to have the same as the example erase later
//
//		double sumc1 = M[1][1] + M[2][1] + M[3][1] + M[0][1];
//		double sumc2 = M[1][2] + M[2][2] + M[3][2] + M[0][2];
//		double sumc3 = M[1][3] + M[2][3] + M[3][3] + M[0][3];
//
//		double x = sumc1 / 4;
//		double y = sumc2 / 4;
//		double z = sumc3 / 4;
//
//		// Parse and write the result
//		double result = stiffFunc(x, y, z);
//
//		double	I = V * result;
//
//
//
//		//The triangle contributes to at most 6 entries in the (upper triangle
//		//of the) stiffness matrix.  We compute these six entries in the
//		//following double loop.
//		for (int s = 0; s < 4; s++)
//		{
//			int lls = Tetrahedra[i][s];
//			if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
//			{
//				for (int r = 0; r <= s; r++)
//				{
//					//If both vertices are free, then there is a contribution
//					//to the stiffness matrix
//					int llr = Tetrahedra[i][r];
//					if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
//					{
//						if (llr <= lls)
//						{
//							
//							size = listInt2DExists(&Positions, llr, lls, size);
//							size = listInt2DExists(&Positions, lls, llr, size);
//							//K[llr][lls] = K[llr][lls] + G[r][s] * I;
//						}
//						else
//						{
//							size = listInt2DExists(&Positions, lls, llr, size);
//							size = listInt2DExists(&Positions, llr, lls, size);
//							//K[lls][llr] = K[lls][llr] + G[r][s] * I;
//						}
//					}
//				}
//			}
//		}
//		//end of principal for    
//		freeMatrix(M, 4);
//		freeMatrix(G, 4);
//		freeMatrix(C, 4);
//		freeMatrix(C34, 3);
//		freeMatrix(C34T, 4);
//	}
//
//	for (int l = 0; l < nodeSize; l++)
//	{
//		/*for (int m = 0; m < nodeSize; m++)
//		{
//			if (K[l][m] != 0)
//			{
//				K[m][l] = K[l][m];
//
//			}
//			else
//			{
//				if (m == l)
//				{
//					K[m][l] = 1;
//				}
//			}
//		}*/
//
//		size = listInt2DExists(&Positions, l, l, size);
//
//	}
//
//	*outPositions = Positions;
//	return size;
//
//
//}
//
//int getPositionsV2(double** Nodes, int** Triangles, int** Tetrahedra,
//	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
//	int constrainedsurface, int*** outPositions)
//{
//	int size = 0;
//	int** Positions;
//	for (int i = 0; i < tetrahedraSize; i++)
//	{
//		
//		//printf("%d tetra \n", i);
//		//The triangle contributes to at most 6 entries in the (upper triangle
//		//of the) stiffness matrix.  We compute these six entries in the
//		//following double loop.
//		for (int s = 0; s < 4; s++)
//		{
//			int lls = Tetrahedra[i][s];
//			if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
//			{
//				for (int r = 0; r <= s; r++)
//				{
//					//If both vertices are free, then there is a contribution
//					//to the stiffness matrix
//					int llr = Tetrahedra[i][r];
//					if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
//					{
//						/*if (llr <= lls)
//						{*/
//
//							size = listInt2DExists(&Positions, llr, lls, size);
//							size = listInt2DExists(&Positions, lls, llr, size);
//							//K[llr][lls] = K[llr][lls] + G[r][s] * I;
//						//}
//						//else
//						//{
//						//	size = listInt2DExists(&Positions, lls, llr, size);
//						//	size = listInt2DExists(&Positions, llr, lls, size);
//						//	//K[lls][llr] = K[lls][llr] + G[r][s] * I;
//						//}
//					}
//				}
//			}
//		}
//		//end of principal for    
//	}
//
//	for (int l = 0; l < nodeSize; l++)
//	{
//		/*for (int m = 0; m < nodeSize; m++)
//		{
//		if (K[l][m] != 0)
//		{
//		K[m][l] = K[l][m];
//
//		}
//		else
//		{
//		if (m == l)
//		{
//		K[m][l] = 1;
//		}
//		}
//		}*/
//
//		size = listInt2DExists(&Positions, l, l, size);
//
//	}
//
//	*outPositions = Positions;
//	return size;
//
//
//}

double* getStiffnessValues(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int** Positions, int sizePositions)
{
	//Supposing that the surfaces does not intersect
	double* stiffValues = newDoubleV(sizePositions);

	for (int i = 0; i < tetrahedraSize; i++)
	{
		//Generate the element matrix for each element
		//[1 x1 y1 z1]
		//[1 x2 y2 z2]
		//[1 x3 y3 z3]
		//[1 x4 y4 z4]
		double** M = newDouble(4, 4);

		for (int l = 0; l < 4; l++)
		{

			for (int m = 0; m < 4; m++)
			{
				if (l == 0)
				{
					M[m][l] = 1;
				}
				else
				{
					M[m][l] = Nodes[Tetrahedra[i][m]][l - 1];
				}

			}

		}

		//Get the coeficientes for each basis function
		//[a1 a2 a3 a4]
		//[b1 b2 b3 b4]
		//[c1 c2 c3 c4]
		//[d1 d2 d3 d4]
		double** C = inverseMatrix(M, 4);

		//The typical integral we must compute is
		//integral over Ti (k(x,y)(grad Ø1).(grad Ø2).
		//The gradients are constant vectors, since the basis functions
		//are linear over this triangle.   So we just compute all the
		//dot products and store them in a matrix G:


		double** C34 = newDouble(3, 4);

		for (int l = 1; l < 4; l++)
		{
			for (int m = 0; m < 4; m++)
			{
				C34[l - 1][m] = C[l][m];

			}
		}


		//multiply gradient(Øi)*gradient(Øj)
		double** C34T = transposeMatrix(C34, 3, 4);
		double** G = multMatrix(C34T, C34, 4, 3, 3, 4);
		//G = C34.T().Mult(C34);
		double detM = det4(M);

		double V = fabs(detM) / 6.0;
		//double V = Math.Abs(M.Determinant()) / 6.0;

		// Add variables and values to hashtable
		//*********************************this is used to have the same as the example erase later

		double sumc1 = M[1][1] + M[2][1] + M[3][1] + M[0][1];
		double sumc2 = M[1][2] + M[2][2] + M[3][2] + M[0][2];
		double sumc3 = M[1][3] + M[2][3] + M[3][3] + M[0][3];

		double x = sumc1 / 4;
		double y = sumc2 / 4;
		double z = sumc3 / 4;

		// Parse and write the result
		double result = stiffFunc(x, y, z);

		double	I = V * result;



		//The triangle contributes to at most 6 entries in the (upper triangle
		//of the) stiffness matrix.  We compute these six entries in the
		//following double loop.
		for (int s = 0; s < 4; s++)
		{
			int lls = Tetrahedra[i][s];
			if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
			{
				for (int r = 0; r <= s; r++)
				{
					//If both vertices are free, then there is a contribution
					//to the stiffness matrix
					int llr = Tetrahedra[i][r];
					if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
					{
						if (llr <= lls)
						{
							for (int o = 0; o < sizePositions; o++)
							{
								if (Positions[o][0] == llr && Positions[o][1] == lls)
								{
									stiffValues[o] += G[r][s] * I;
								}
								
							}
						}
						else
						{
							for (int o = 0; o < sizePositions; o++)
							{
								if (Positions[o][0] == lls && Positions[o][1] == llr)
								{
									stiffValues[o] += G[r][s] * I;
								}
								
									
							}
						}
					}
				}
			}
		}
		//end of principal for    
		freeMatrix(M, 4);
		freeMatrix(G, 4);
		freeMatrix(C, 4);
		freeMatrix(C34, 3);
		freeMatrix(C34T, 4);
	}


	for (int o = 0; o < sizePositions; o++)
	{
		if (Positions[o][0] == Positions[o][1] && stiffValues[o] == 0)
			stiffValues[o] = 1;
	}

	for (int i = 0; i < sizePositions; i++)
	{
		if (stiffValues[i] == 0)
			for (int j = 0; j < sizePositions; j++)
			{
				if (Positions[i][0] == Positions[j][1] && Positions[i][1] == Positions[j][0])
					stiffValues[i] = stiffValues[j];
			}
	}

	//for (int l = 0; l < nodeSize; l++)
	//{
	//	for (int m = 0; m < nodeSize; m++)
	//	{
	//		if (K[l][m] != 0)
	//		{
	//			K[m][l] = K[l][m];

	//		}
	//		else
	//		{
	//			if (m == l)
	//			{
	//				K[m][l] = 1;
	//			}
	//		}
	//	}



	//}


	return stiffValues;


}

int getPositions(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int*** outPositions)
{
	int size = 0;

	int** positionsTemp = newIntValue(tetrahedraSize*32, 2,-1);
	int positionsCount = 0;

	
	for (int i = 0; i < tetrahedraSize; i++)
	{

		//printf("%d tetra \n", i);
		//The triangle contributes to at most 6 entries in the (upper triangle
		//of the) stiffness matrix.  We compute these six entries in the
		//following double loop.
		for (int s = 0; s < 4; s++)
		{
			int lls = Tetrahedra[i][s];
			if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
			{
				for (int r = 0; r <= s; r++)
				{
					//If both vertices are free, then there is a contribution
					//to the stiffness matrix
					int llr = Tetrahedra[i][r];
					if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
					{
						
						positionsTemp[positionsCount][0] = llr;
						positionsTemp[positionsCount][1] = lls;
							positionsCount++;
						positionsTemp[positionsCount][0] = lls;
						positionsTemp[positionsCount][1] = llr;
							positionsCount++;
					}
				}
			}
		}
		//end of principal for    
	}

	for (int l = 0; l < nodeSize; l++)
	{
		positionsTemp[positionsCount][0] = l;
		positionsTemp[positionsCount][1] = l;
		positionsCount++;
	}

	printf("connections count %d \n", positionsCount);

	int** positionsTemp2 = newIntValue(positionsCount, 2, -1);
	int positionsCount2 = 0;


	for (int i = 0; i < positionsCount; i++)
	{
		int pos0 = positionsTemp[i][0];
		int pos1 = positionsTemp[i][1];
		int exists = 0;
		if (i != 0)
		{
			for (int j = 0; j < i; j++)
			{
				if (pos0 == positionsTemp[j][0] && pos1 == positionsTemp[j][1])
				{
					exists = 1;
					break;
				}
			}
		}
		
		if (exists == 0)
		{
			positionsTemp2[positionsCount2][0] = positionsTemp[i][0];
			positionsTemp2[positionsCount2][1] = positionsTemp[i][1];
			positionsCount2++;
		}
	}
	printf("non repeated connections %d \n", positionsCount2);

	int** Positions=newIntValue(positionsCount2, 2, -1);

	for (int l = 0; l < positionsCount2; l++)
	{
		Positions[l][0] = positionsTemp2[l][0];
		Positions[l][1] = positionsTemp2[l][1];
		//size = listInt2DExists(&Positions, positionsTemp2[l][0], positionsTemp2[l][1], size);
	}

	freeMatrixInt(positionsTemp, tetrahedraSize * 20);
	freeMatrixInt(positionsTemp2, positionsCount);

	*outPositions = Positions;
	return positionsCount2;


}


double* getMassValues(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int** Positions, int sizePositions)
	{
		//Supposing that the surfaces does not intersect
		double* massValues = newDoubleV(sizePositions);

		for (int i = 0; i < tetrahedraSize; i++)
		{
			//Generate the element matrix for each element
			//[1 x1 y1 z1]
			//[1 x2 y2 z2]
			//[1 x3 y3 z3]
			//[1 x4 y4 z4]
			double** M = newDouble(4, 4);

			for (int l = 0; l < 4; l++)
			{

				for (int m = 0; m < 4; m++)
				{
					if (l == 0)
					{
						M[m][l] = 1;
					}
					else
					{
						M[m][l] = Nodes[Tetrahedra[i][m]][l - 1];
					}

				}

			}

			double** G = newDouble(4, 4);

			//ERROR***********
			for (int l = 0; l < 4; l++)
			{
				for (int m = 0; m < 4; m++)
				{
					if (m == l)
					{
						G[m][l] = 3.0 / 9.0;

					}
					else
					{
						G[m][l] = 2.0 / 9.0;
					}
				}
			}
			double** Gtemp = multConstant(G, 0.25, 4, 4);
			copyMatrix(G, Gtemp, 4, 4);
			freeMatrix(Gtemp, 4);


			//ERROR***********
			double detM = det4(M);
			double V = fabs(detM) / 6.0;

			// Add variables and values to hashtable
			//*********************************this is used to have the same as the example erase later

			double I;


			double sumc1 = M[1][1] + M[2][1] + M[3][1] + M[0][1];
			double sumc2 = M[1][2] + M[2][2] + M[3][2] + M[0][2];
			double sumc3 = M[1][3] + M[2][3] + M[3][3] + M[0][3];

			double x = sumc1 / 4;
			double y = sumc2 / 4;
			double z = sumc3 / 4;

			double result = massFunc(x, y, z);

			I = V * result;



			//The triangle contributes to at most 6 entries in the (upper triangle
			//of the) stiffness matrix.  We compute these six entries in the
			//following double loop.
			for (int s = 0; s < 4; s++)
			{
				int lls = Tetrahedra[i][s];
				if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
				{
					for (int r = 0; r <= s; r++)
					{
						//If both vertices are free, then there is a contribution
						//to the stiffness matrix
						int llr = Tetrahedra[i][r];
						if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
						{
							if (llr <= lls)
							{
								for (int o = 0; o < sizePositions; o++)
								{
									if (Positions[o][0] == llr && Positions[o][1] == lls)
									{
										massValues[o] += G[r][s] * I;
									}

								}
							}
							else
							{
								for (int o = 0; o < sizePositions; o++)
								{
									if (Positions[o][0] == lls && Positions[o][1] == llr)
									{
										massValues[o] += G[r][s] * I;
									}


								}
							}
						}
					}
				}
			}
			//end of principal for    
			freeMatrix(M, 4);
			freeMatrix(G, 4);
			
		}


		for (int o = 0; o < sizePositions; o++)
		{
			if (Positions[o][0] == Positions[o][1] && massValues[o] == 0)
				massValues[o] = 1;
		}

		for (int i = 0; i < sizePositions; i++)
		{
			if (massValues[i] == 0)
				for (int j = 0; j < sizePositions; j++)
				{
					if (Positions[i][0] == Positions[j][1] && Positions[i][1] == Positions[j][0])
						massValues[i] = massValues[j];
				}
		}

		

		return massValues;


	}

	double* getStiffnessValuesTensor(double** Nodes, int** Triangles, int** Tetrahedra,
		int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
		int constrainedsurface, int** Positions, int sizePositions, 
		double alongfibers_cond, double perpenfibers_cond, double normalfibers_cond, double** fibers)
	{
		//Supposing that the surfaces does not intersect
		double* stiffValues = newDoubleV(sizePositions);

		for (int i = 0; i < tetrahedraSize; i++)
		{
			//Generate the element matrix for each element
			//[1 x1 y1 z1]
			//[1 x2 y2 z2]
			//[1 x3 y3 z3]
			//[1 x4 y4 z4]
			double** M = newDouble(4, 4);

			for (int l = 0; l < 4; l++)
			{

				for (int m = 0; m < 4; m++)
				{
					if (l == 0)
					{
						M[m][l] = 1;
					}
					else
					{
						M[m][l] = Nodes[Tetrahedra[i][m]][l - 1];
					}

				}

			}

			//Get the coeficientes for each basis function
			//[a1 a2 a3 a4]
			//[b1 b2 b3 b4]
			//[c1 c2 c3 c4]
			//[d1 d2 d3 d4]
			double** C = inverseMatrix(M, 4);

			//The typical integral we must compute is
			//integral over Ti (k(x,y)(grad Ø1).(grad Ø2).
			//The gradients are constant vectors, since the basis functions
			//are linear over this triangle.   So we just compute all the
			//dot products and store them in a matrix G:


			double** C34 = newDouble(3, 4);

			for (int l = 1; l < 4; l++)
			{
				for (int m = 0; m < 4; m++)
				{
					C34[l - 1][m] = C[l][m];

				}
			}


			// Parse and write the result
			double** A = newDouble(3, 3);



			for (int m = 0; m < 4; m++)
			{
				//in the fiber direction
				A[0][0] += 0.25 * fibers[Tetrahedra[i][m]][0];
				A[1][0] += 0.25 * fibers[Tetrahedra[i][m]][1];
				A[2][0] += 0.25 * fibers[Tetrahedra[i][m]][2];
				//perpendicular to the fibers
				A[0][1] += 0.25 * fibers[Tetrahedra[i][m]][3];
				A[1][1] += 0.25 * fibers[Tetrahedra[i][m]][4];
				A[2][1] += 0.25 * fibers[Tetrahedra[i][m]][5];
				//normal to the fibers
				A[0][2] += 0.25 * fibers[Tetrahedra[i][m]][6];
				A[1][2] += 0.25 * fibers[Tetrahedra[i][m]][7];
				A[2][2] += 0.25 * fibers[Tetrahedra[i][m]][8];
			}

			double** Mi = newDouble(3, 3);

			Mi[0][0] = alongfibers_cond;
			Mi[1][1] = perpenfibers_cond;
			Mi[2][2] = normalfibers_cond;

			// double[,] conductivity = A.Mult(Mi.Mult(A.T()));
			double** At = transposeMatrix(A,3,3);
			double** MiAt = multMatrix(Mi, At, 3, 3,3,3);
			double** conductivity = multMatrix(A, MiAt, 3, 3, 3, 3);
			double** conductivityC34 = multMatrix(conductivity, C34, 3, 3, 3, 4);


			//multiply gradient(Øi)*gradient(Øj)
			double** C34T = transposeMatrix(C34, 3, 4);
			double** G = multMatrix(C34T, conductivityC34, 4, 3, 3, 4);

			freeMatrix(A, 3);
			freeMatrix(Mi, 3);
			freeMatrix(At, 3);
			freeMatrix(MiAt, 3);
			freeMatrix(conductivity, 3);
			freeMatrix(conductivityC34, 3);
			//G = C34.T().Mult(C34);
			double detM = det4(M);

			double V = fabs(detM) / 6.0;
			//double V = Math.Abs(M.Determinant()) / 6.0;

			


			double	I = V ;



			//The triangle contributes to at most 6 entries in the (upper triangle
			//of the) stiffness matrix.  We compute these six entries in the
			//following double loop.
			for (int s = 0; s < 4; s++)
			{
				int lls = Tetrahedra[i][s];
				if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
				{
					for (int r = 0; r <= s; r++)
					{
						//If both vertices are free, then there is a contribution
						//to the stiffness matrix
						int llr = Tetrahedra[i][r];
						if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
						{
							if (llr <= lls)
							{
								for (int o = 0; o < sizePositions; o++)
								{
									if (Positions[o][0] == llr && Positions[o][1] == lls)
									{
										stiffValues[o] += G[r][s] * I;
									}

								}
							}
							else
							{
								for (int o = 0; o < sizePositions; o++)
								{
									if (Positions[o][0] == lls && Positions[o][1] == llr)
									{
										stiffValues[o] += G[r][s] * I;
									}


								}
							}
						}
					}
				}
			}
			//end of principal for    
			freeMatrix(M, 4);
			freeMatrix(G, 4);
			freeMatrix(C, 4);
			freeMatrix(C34, 3);
			freeMatrix(C34T, 4);
		}


		for (int o = 0; o < sizePositions; o++)
		{
			if (Positions[o][0] == Positions[o][1] && stiffValues[o] == 0)
				stiffValues[o] = 1;
		}

		for (int i = 0; i < sizePositions; i++)
		{
			if (stiffValues[i] == 0)
				for (int j = 0; j < sizePositions; j++)
				{
					if (Positions[i][0] == Positions[j][1] && Positions[i][1] == Positions[j][0])
						stiffValues[i] = stiffValues[j];
				}
		}

		//for (int l = 0; l < nodeSize; l++)
		//{
		//	for (int m = 0; m < nodeSize; m++)
		//	{
		//		if (K[l][m] != 0)
		//		{
		//			K[m][l] = K[l][m];

		//		}
		//		else
		//		{
		//			if (m == l)
		//			{
		//				K[m][l] = 1;
		//			}
		//		}
		//	}



		//}


		return stiffValues;


	}

	double* getStiffnessValuesTensorLesion(double** Nodes, int** Triangles, int** Tetrahedra,
		int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
		int constrainedsurface, int** Positions, int sizePositions,
		double alongfibers_cond, double perpenfibers_cond, double normalfibers_cond, double** fibers, 
		int lesionPosition, double reductionCoeff)
	{
		//Supposing that the surfaces does not intersect
		double* stiffValues = newDoubleV(sizePositions);

		for (int i = 0; i < tetrahedraSize; i++)
		{
			//Generate the element matrix for each element
			//[1 x1 y1 z1]
			//[1 x2 y2 z2]
			//[1 x3 y3 z3]
			//[1 x4 y4 z4]
			double** M = newDouble(4, 4);

			for (int l = 0; l < 4; l++)
			{

				for (int m = 0; m < 4; m++)
				{
					if (l == 0)
					{
						M[m][l] = 1;
					}
					else
					{
						M[m][l] = Nodes[Tetrahedra[i][m]][l - 1];
					}

				}

			}

			//Get the coeficientes for each basis function
			//[a1 a2 a3 a4]
			//[b1 b2 b3 b4]
			//[c1 c2 c3 c4]
			//[d1 d2 d3 d4]
			double** C = inverseMatrix(M, 4);

			//The typical integral we must compute is
			//integral over Ti (k(x,y)(grad Ø1).(grad Ø2).
			//The gradients are constant vectors, since the basis functions
			//are linear over this triangle.   So we just compute all the
			//dot products and store them in a matrix G:


			double** C34 = newDouble(3, 4);

			for (int l = 1; l < 4; l++)
			{
				for (int m = 0; m < 4; m++)
				{
					C34[l - 1][m] = C[l][m];

				}
			}


			// Parse and write the result
			double** A = newDouble(3, 3);



			for (int m = 0; m < 4; m++)
			{
				//in the fiber direction
				A[0][0] += 0.25 * fibers[Tetrahedra[i][m]][0];
				A[1][0] += 0.25 * fibers[Tetrahedra[i][m]][1];
				A[2][0] += 0.25 * fibers[Tetrahedra[i][m]][2];
				//perpendicular to the fibers
				A[0][1] += 0.25 * fibers[Tetrahedra[i][m]][3];
				A[1][1] += 0.25 * fibers[Tetrahedra[i][m]][4];
				A[2][1] += 0.25 * fibers[Tetrahedra[i][m]][5];
				//normal to the fibers
				A[0][2] += 0.25 * fibers[Tetrahedra[i][m]][6];
				A[1][2] += 0.25 * fibers[Tetrahedra[i][m]][7];
				A[2][2] += 0.25 * fibers[Tetrahedra[i][m]][8];
			}

			double** Mi = newDouble(3, 3);

			Mi[0][0] = alongfibers_cond;
			Mi[1][1] = perpenfibers_cond;
			Mi[2][2] = normalfibers_cond;

			// double[,] conductivity = A.Mult(Mi.Mult(A.T()));
			double** At = transposeMatrix(A, 3, 3);
			double** MiAt = multMatrix(Mi, At, 3, 3, 3, 3);
			double** conductivity = multMatrix(A, MiAt, 3, 3, 3, 3);
			double** conductivityC34 = multMatrix(conductivity, C34, 3, 3, 3, 4);


			//multiply gradient(Øi)*gradient(Øj)
			double** C34T = transposeMatrix(C34, 3, 4);
			double** G = multMatrix(C34T, conductivityC34, 4, 3, 3, 4);

			freeMatrix(A, 3);
			freeMatrix(Mi, 3);
			freeMatrix(At, 3);
			freeMatrix(MiAt, 3);
			freeMatrix(conductivity, 3);
			freeMatrix(conductivityC34, 3);
			//G = C34.T().Mult(C34);
			double detM = det4(M);

			double V = fabs(detM) / 6.0;
			//double V = Math.Abs(M.Determinant()) / 6.0;




			double	I = V;


			if (
				(Tetrahedra[i][ 0] == lesionPosition) ||
				(Tetrahedra[i][1] == lesionPosition) ||
				(Tetrahedra[i][2] == lesionPosition) ||
				(Tetrahedra[i][3] == lesionPosition)
				)
				I = I *reductionCoeff;


			//The triangle contributes to at most 6 entries in the (upper triangle
			//of the) stiffness matrix.  We compute these six entries in the
			//following double loop.
			for (int s = 0; s < 4; s++)
			{
				int lls = Tetrahedra[i][s];
				if (NodeSurface[Tetrahedra[i][s]] != constrainedsurface)
				{
					for (int r = 0; r <= s; r++)
					{
						//If both vertices are free, then there is a contribution
						//to the stiffness matrix
						int llr = Tetrahedra[i][r];
						if (NodeSurface[Tetrahedra[i][r]] != constrainedsurface)
						{
							if (llr <= lls)
							{
								for (int o = 0; o < sizePositions; o++)
								{
									if (Positions[o][0] == llr && Positions[o][1] == lls)
									{
										stiffValues[o] += G[r][s] * I;
									}

								}
							}
							else
							{
								for (int o = 0; o < sizePositions; o++)
								{
									if (Positions[o][0] == lls && Positions[o][1] == llr)
									{
										stiffValues[o] += G[r][s] * I;
									}


								}
							}
						}
					}
				}
			}
			//end of principal for    
			freeMatrix(M, 4);
			freeMatrix(G, 4);
			freeMatrix(C, 4);
			freeMatrix(C34, 3);
			freeMatrix(C34T, 4);
		}


		for (int o = 0; o < sizePositions; o++)
		{
			if (Positions[o][0] == Positions[o][1] && stiffValues[o] == 0)
				stiffValues[o] = 1;
		}

		for (int i = 0; i < sizePositions; i++)
		{
			if (stiffValues[i] == 0)
				for (int j = 0; j < sizePositions; j++)
				{
					if (Positions[i][0] == Positions[j][1] && Positions[i][1] == Positions[j][0])
						stiffValues[i] = stiffValues[j];
				}
		}

		//for (int l = 0; l < nodeSize; l++)
		//{
		//	for (int m = 0; m < nodeSize; m++)
		//	{
		//		if (K[l][m] != 0)
		//		{
		//			K[m][l] = K[l][m];

		//		}
		//		else
		//		{
		//			if (m == l)
		//			{
		//				K[m][l] = 1;
		//			}
		//		}
		//	}



		//}


		return stiffValues;


	}

