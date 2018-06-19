#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "geometry.h"
#include "FEM.h"
#include "../Amazigh/printFile.h"
#include "../Amazigh/mathVector.h"
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/solvers.h"
#include "../Amazigh/plot2D.h"
#include "../Amazigh/supportC.h"
#define PI 3.141592653589793238462643 
int globalExample = 999;

double stiffFunc(double x, double y, double z)
{
	double value;
	if (globalExample == 1)
	{
		value = 1;
	}
	if (globalExample == 2)
	{
		value = (x + 2);
	}
	if (globalExample == 10)
	{
		value = 1.0;
	}
	if (globalExample == 20)
	{
		value = 1.0;
	}
	return value;
}

double ux(double x, double y, double z, double t)
{
	double value;
	if (globalExample == 1)
	{
		value = 2 * (x - 2)*(y - 2)*(y - 2)*(z - 2)*(z - 2);
	}
	if (globalExample == 2)
	{
		value = 2 * (x - 2)*(y - 2)*(y - 2)*(z - 2)*(z - 2);
	}
	if (globalExample == 10)
	{
		value= -(2 * x + 2 * y + 2 * z) / (8 * t * exp(((x + y + z) * (x + y + z)) / (4 * t)) * sqrt((PI * t)));
	}
	if (globalExample == 20)
	{
		value = -(2 * x + 2 * y + 2 * z) / (8 * t * exp(((x + y + z) * (x + y + z)) / (4 * t)) * sqrt((PI * t)));
	}
	return value;
}

double uy(double x, double y, double z, double t)
{
	double value;
	if (globalExample == 1)
	{
		value = 2 * (y - 2)*(x - 2)*(x - 2)*(z - 2)*(z - 2);
	}
	if (globalExample == 2)
	{
		value = 2 * (y - 2)*(x - 2)*(x - 2)*(z - 2)*(z - 2);
	}
	if (globalExample == 10)
	{
		value = -(2 * x + 2 * y + 2 * z) / (8 * t * exp(((x + y + z) * (x + y + z)) / (4 * t)) * sqrt((PI * t)));
	}
	if (globalExample == 20)
	{
		value = -(2 * x + 2 * y + 2 * z) / (8 * t * exp(((x + y + z) * (x + y + z)) / (4 * t)) * sqrt((PI * t)));
	}
	return value;
}

double uz(double x, double y, double z, double t)
{
	double value;
	if (globalExample == 1)
	{
		value = 2 * (z - 2)*(x - 2)*(x - 2)*(y - 2)*(y - 2);
	}
	if (globalExample == 2)
	{
		value = 2 * (z - 2)*(x - 2)*(x - 2)*(y - 2)*(y - 2);
	}
	if (globalExample == 10)
	{
		value = -(2 * x + 2 * y + 2 * z) / (8 * t * exp(((x + y + z) * (x + y + z)) / (4 * t)) * sqrt((PI * t)));
	}
	if (globalExample == 20)
	{
		value = -(2 * x + 2 * y + 2 * z) / (8 * t * exp(((x + y + z) * (x + y + z)) / (4 * t)) * sqrt((PI * t)));
	}
	return value;
}

double fFunc(double x, double y, double z, double t)
{
	double value;
	if (globalExample == 1)
	{
		value = -2 * (y - 2)*(y - 2)*(z - 2)*(z - 2) -
			2 * (x - 2)*(x - 2)*(z - 2)*(z - 2) - 
			2 * (x - 2)*(x - 2)*(y - 2)*(y - 2);
	}
	if (globalExample == 2)
	{
		value = -(2 * x) * 2 * (y - 2)*(y - 2)*(z - 2)*(z - 2) - 
			(x + 2) * 2 * (x - 2)*(x - 2)*(z - 2)*(z - 2) - 
			(x + 2) * 2 * (x - 2)*(x - 2)*(y - 2)*(y - 2);
	}
	if (globalExample == 10)
	{

		value = (exp(-pow((x + y + z), 2) / (4 * t))*pow((x + y + z), 2)) / (8 * pow(PI, (1 / 2))*pow(t, (5 / 2))) - exp(-pow((x + y + z), 2) / (4 * t)) / (4 * pow(PI, (1 / 2))*pow(t, (3 / 2)));
		value = -(3.0/4.0)*value;
	}
	if (globalExample == 20)
	{

		value = 0;
	}
	return value;


}

double uFunc(double x, double y, double z, double t)
{
	double value;
	if (globalExample == 1)
	{
		value = (x - 2)*(x - 2)*(y - 2)*(y - 2)*(z - 2)*(z - 2);
	}
	if (globalExample == 2)
	{
		value = (x - 2)*(x - 2)*(y - 2)*(y - 2)*(z - 2)*(z - 2);
	}
	if (globalExample == 10)
	{
		double u = x + y + z;
		value = (exp((-(u * u)) / (4 * t))) / (2 * sqrt(PI * t));
	}
	if (globalExample == 20)
	{
		double u = x + y + z;
		value = (exp((-(u * u)) / (4 * t))) / (2 * sqrt(PI * t));
	}

	
	return value;
}

void FEMTest1()
{
	globalExample = 1;
	geometryVolume geo = openVolume("../Geometries/test.msh");

	//Stiffness Matrix
	double** K = Stiffness3DforFEM(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 0);
	//Neuman Data
	int ntriangles = neumanTriangles(geo.triangleSurface, geo.triangleSize, 1);
	int** Neuman = newInt(ntriangles, 3);
	double** NeumanData = newDouble(ntriangles, 3);
	getNeumanData3DNEW(geo.Nodes, geo.Triangles, geo.nodeSize, geo.triangleSize,
		geo.triangleSurface, 1, NeumanData, Neuman, ntriangles, 0);
	//Dirichlet Data
	double* nodeDirichlet = getDirichletData(geo.Nodes, geo.nodeSize, 0, geo.nodeSurface,0);
	//Load Data
	double* Load = Load3D(geo.Nodes, geo.Tetrahedras, geo.nodeSurface, geo.nodeSize,
		geo.tetrahedraSize, 0, nodeDirichlet, NeumanData, Neuman, ntriangles, 0 );

	double* x = PBCGSTABC(K, Load, geo.nodeSize);
	//printVector(x, geo.nodeSize, "x.txt");
	//printVector(Load, geo.nodeSize, "Load.txt");

	double* exact = newDoubleV(geo.nodeSize);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		double x = geo.Nodes[i][0];
		double y = geo.Nodes[i][1];
		double z = geo.Nodes[i][2];

		exact[i] = uFunc(x, y, z,0);
	}

	double** Values = newDouble(geo.nodeSize, 2);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		Values[i][0] = exact[i];
		Values[i][1] = x[i];
	}
	char** Names[2];
	Names[0] = "Exact";
	Names[1] = "Calculated";

	plot2DComparisonNamesNew(Values, 2, geo.nodeSize, "Results", "Graph.html", Names, "Value");

	printVolumeValues("Results.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values, 2, "Data");
}

void FEMTest2()
{
	globalExample = 1;
	geometryVolume geo = openVolume("../Geometries/test.msh");
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);
	//Stiffness Matrix

	int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 0,
		&Positions);
	printf(" Positions %d \n", sizePositions);

	double* stiffValues=getStiffnessValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		 0,  Positions,  sizePositions);
	printf("stiffness Matrix \n");
	//printVector(stiffValues, sizePositions, "stiffValues.txt");

    //Neuman Data
	int ntriangles = neumanTriangles(geo.triangleSurface, geo.triangleSize, 1);
	int** Neuman = newInt(ntriangles, 3);
	double** NeumanData = newDouble(ntriangles, 3);
	getNeumanData3DNEW(geo.Nodes, geo.Triangles, geo.nodeSize, geo.triangleSize,
		geo.triangleSurface, 1, NeumanData, Neuman, ntriangles, 0);
	//Dirichlet Data
	double* nodeDirichlet = getDirichletData(geo.Nodes, geo.nodeSize, 0, geo.nodeSurface,0);
	//Load Data
	double* Load = Load3D(geo.Nodes, geo.Tetrahedras, geo.nodeSurface, geo.nodeSize,
		geo.tetrahedraSize, 0, nodeDirichlet, NeumanData, Neuman, ntriangles, 0);
	printf("Load Vector \n");
	double* xSparse = PBCGSTABCSparse(stiffValues, Positions, sizePositions, Load, geo.nodeSize);
	//printVector(xSparse, geo.nodeSize, "xSparse.txt");

	double* exact = newDoubleV(geo.nodeSize);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		double x = geo.Nodes[i][0];
		double y = geo.Nodes[i][1];
		double z = geo.Nodes[i][2];

		exact[i] = uFunc(x, y, z, 0);
	}

	double** Values = newDouble(geo.nodeSize, 2);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		Values[i][0] = exact[i];
		Values[i][1] = xSparse[i];
	}
	char** Names[2];
	Names[0] = "Exact";
	Names[1] = "Calculated";

	plot2DComparisonNamesNew(Values, 2, geo.nodeSize, "Results", "Graph.html", Names,"Value");

	printVolumeValues("Results.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values, 2, "Data");

	
}

void FEMTest3()
{
	globalExample = 2;
	geometryVolume geo = openVolume("../Geometries/test3.msh");
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);
	//Stiffness Matrix

	int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 0,
		&Positions);
	printf(" Positions %d \n", sizePositions);

	double* stiffValues = getStiffnessValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		0, Positions, sizePositions);
	printf("stiffness Matrix \n");
	//printVector(stiffValues, sizePositions, "stiffValues.txt");

	//Neuman Data
	int ntriangles = neumanTriangles(geo.triangleSurface, geo.triangleSize, 1);
	int** Neuman = newInt(ntriangles, 3);
	double** NeumanData = newDouble(ntriangles, 3);
	getNeumanData3DNEW(geo.Nodes, geo.Triangles, geo.nodeSize, geo.triangleSize,
		geo.triangleSurface, 1, NeumanData, Neuman, ntriangles, 0);
	//Dirichlet Data
	double* nodeDirichlet = getDirichletData(geo.Nodes, geo.nodeSize, 0, geo.nodeSurface,0);
	//Load Data
	double* Load = Load3D(geo.Nodes, geo.Tetrahedras, geo.nodeSurface, geo.nodeSize,
		geo.tetrahedraSize, 0, nodeDirichlet, NeumanData, Neuman, ntriangles, 0);
	printf("Load Vector \n");
	double* xSparse = PBCGSTABCSparse(stiffValues, Positions, sizePositions, Load, geo.nodeSize);
	//printVector(xSparse, geo.nodeSize, "xSparse.txt");

	double* exact = newDoubleV(geo.nodeSize);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		double x = geo.Nodes[i][0];
		double y = geo.Nodes[i][1];
		double z = geo.Nodes[i][2];

		exact[i] = uFunc(x,y,z, 0);
	}

	double** Values = newDouble(geo.nodeSize, 2);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		Values[i][0] = exact[i];
		Values[i][1] = xSparse[i];
	}
	char** Names[2];
	Names[0] = "Exact";
	Names[1] = "Calculated";

	plot2DComparisonNamesNew(Values, 2, geo.nodeSize, "Results", "Graph.html", Names, "Value");

	printVolumeValues("Results.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values, 2, "Data");


}

double** Heateq3D(double** Nodes,  int nodeSize, double t)
{

	double* nodeDirichlet = newDoubleV(nodeSize);

	for (int i = 0; i < nodeSize; i++)
	{

		double x;
		double y;
		double z;
		double result;
		x = Nodes[i][ 0];
		y = Nodes[i][1];
		z = Nodes[i][2];
		double u = x + y + z;
		result = (exp((-(u * u)) / (4 * t))) / (2 * sqrt(PI * t));

		nodeDirichlet[i] = result;

	}


	return nodeDirichlet;

}

void FEMTestTime1()
{
	globalExample = 10;
	geometryVolume geo = openVolume("../Geometries/test.msh");

	//Mass Matrix
	double** M = Mass3DforFEM(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 5);
	//Stiffness Matrix
	double** K = Stiffness3DforFEM(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 5);
	
	
	double dt = 0.01;
	double heatCoef = 1 / 3;
	double coef = dt *heatCoef;

	double** coefStiff = multConstant(K, coef, geo.nodeSize, geo.nodeSize);
	double** A = addMatrix(M, coefStiff, geo.nodeSize, geo.nodeSize);

	double limit = 2.5;

	double t = 1;
	double* u0 = Heateq3D(geo.Nodes, geo.nodeSize, t);
	double* v = newDoubleV(geo.nodeSize);
	double* vant = multVector(M, u0, geo.nodeSize, geo.nodeSize);
	int steps = (int)((limit - t) / dt);

	int k = 0;
	double** Values = newDouble(geo.nodeSize, steps + 1);
	//NodeSurface = NodeSurface.Initialize(1);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		geo.nodeSurface[i] = 1;
	}

	printf("Label 0 ");

	
	printf("steps %d \n", steps);
	do
	{

		//Neuman Data
		int ntriangles = neumanTriangles(geo.triangleSurface, geo.triangleSize, 1);
		int** Neuman = newInt(ntriangles, 3);
		double** NeumanData = newDouble(ntriangles, 3);
		getNeumanData3DNEW(geo.Nodes, geo.Triangles, geo.nodeSize, geo.triangleSize,
			geo.triangleSurface, 1, NeumanData, Neuman, ntriangles, t);
		//Dirichlet Data
		double* nodeDirichlet = getDirichletData(geo.Nodes, geo.nodeSize, 0, geo.nodeSurface,0);
		//Load Data
		double* Load = Load3D(geo.Nodes, geo.Tetrahedras, geo.nodeSurface, geo.nodeSize,
			geo.tetrahedraSize, 5, nodeDirichlet, NeumanData, Neuman, ntriangles, t);

		double* coefF = multConstantV(Load, coef, geo.nodeSize);
		double* vantTemp = addVector(vant, coefF,geo.nodeSize);
		double* v = PBCGSTABC(A, vantTemp, geo.nodeSize);

		

		for (int i = 0; i < geo.nodeSize; i++)
		{
			Values[i][ k] = v[i];

		}

		double* Mv = multVector(M, v, geo.nodeSize, geo.nodeSize);

		copyVector(vant, Mv, geo.nodeSize);

		t += dt;
		printf("t %.11f \n", t);

		k++;

		
		free(nodeDirichlet);
		free(Load);
		free(coefF);
		free(vantTemp);
		free(v);
		free(Mv);
		freeMatrixInt(Neuman, ntriangles);
		freeMatrix(NeumanData, ntriangles);

	} while (t < limit);

	
	
	printf("Label 1 ");

	double* exact = Heateq3D(geo.Nodes,geo.nodeSize,t-dt);

	double** Values2 = newDouble(geo.nodeSize, 2);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		Values2[i][0] = exact[i];
		Values2[i][1] = Values[i][k-1];
	}
	char** Names[2];
	Names[0] = "Exact";
	Names[1] = "Calculated";

	plot2DComparisonNamesNew(Values2, 2, geo.nodeSize, "Results", "Graph.html", Names, "Value");

	printVolumeValues("Results.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values2, 2, "Data");
}

void FEMTestTime2()
{
	//that is why it is a 1/3
	//diff(u, t) - 1 / 3(diff(ux, x) + diff(uy, y) + diff(uz, z)) = 0;

	globalExample = 10;
	geometryVolume geo = openVolume("../../Geometries/test.msh");

	//Mass Matrix
	double** M = Mass3DforFEM(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 0);
	//Stiffness Matrix
	double** K = Stiffness3DforFEM(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 0);


	double dt = 0.01;
	double heatCoef = 1.0 / 3.0;
	double coef = dt *heatCoef;

	double** A = newDouble(geo.nodeSize, geo.nodeSize);

	for (int i = 0; i < geo.nodeSize; i++)
	{
		for (int j = 0; j < geo.nodeSize; j++)
		{
			
				A[i][j] = M[i][ j] + K[i][j] * coef;
			
		}
	}

	double limit = 2.5;

	double t = 1;
	double* u0 = Heateq3D(geo.Nodes, geo.nodeSize, t);
	double* v = newDoubleV(geo.nodeSize);
	double* vant = multVector(M, u0, geo.nodeSize, geo.nodeSize);
	int steps = (int)((limit - t) / dt);

	int k = 0;
	double** Values = newDouble(geo.nodeSize, steps + 1);
	
	//printVectorInt(geo.nodeSurface, geo.nodeSize, "Nodesurface.txt");

	printf("steps %d \n", steps);
	do
	{

		//Neuman Data
		int ntriangles = neumanTriangles(geo.triangleSurface, geo.triangleSize, 1);
		int** Neuman = newInt(ntriangles, 3);
		double** NeumanData = newDouble(ntriangles, 3);
		getNeumanData3DNEW(geo.Nodes, geo.Triangles, geo.nodeSize, geo.triangleSize,
			geo.triangleSurface, 1, NeumanData, Neuman, ntriangles, t);
		//Dirichlet Data
		double* nodeDirichlet = getDirichletData(geo.Nodes, geo.nodeSize, 0, geo.nodeSurface,t);
		//Load Data
		double* Load = Load3D(geo.Nodes, geo.Tetrahedras, geo.nodeSurface, geo.nodeSize,
			geo.tetrahedraSize, 0, nodeDirichlet, NeumanData, Neuman, ntriangles, t);

		double* coefF = multConstantV(Load, coef, geo.nodeSize);
		double* vantTemp = addVector(vant, coefF, geo.nodeSize);

		
		double* v = PBCGSTABC(A, vantTemp, geo.nodeSize);


		for (int i = 0; i < geo.nodeSize; i++)
		{
			if (geo.nodeSurface[i] == 0)
			{
				Values[i][k] = nodeDirichlet[i];
			}
			else
			{
				Values[i][k] = v[i];
			}
				
			

		}

		double* Mv = multVector(M, v, geo.nodeSize, geo.nodeSize);

		copyVector(vant, Mv, geo.nodeSize);

		t += dt;
		printf("t %.11f \n", t);

		k++;


		free(nodeDirichlet);
		free(Load);
		free(coefF);
		free(vantTemp);
		free(v);
		free(Mv);
		freeMatrixInt(Neuman, ntriangles);
		freeMatrix(NeumanData, ntriangles);

	} while (t < limit);



	printf("Label 1 ");

	double* exact = Heateq3D(geo.Nodes, geo.nodeSize, t-dt );

	

	double** Values2 = newDouble(geo.nodeSize, 2);
	double* result = newDoubleV(geo.nodeSize);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		Values2[i][0] = exact[i];
		Values2[i][1] = Values[i][k - 1];
		if (geo.nodeSurface[i]==1)
		{
			result[i] = Values[i][k - 1];
		}
		else
		{
			exact[i] = 0;
		}
		
	}

	double* subResults = subVector(result, exact, geo.nodeSize);
	double normResults = norm(subResults, geo.nodeSize);
	printf("Norm %.10f \n", normResults);

	char** Names[2];
	Names[0] = "Exact";
	Names[1] = "Calculated";

	plot2DComparisonNamesNew(Values2, 2, geo.nodeSize, "Results", "Graph.html", Names, "Value");

	printVolumeValues("Results.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values2, 2, "Data");
}

void FEMTestTime3()
{
	//that is why it is a 1/3
	//diff(u, t) - 1 / 3(diff(ux, x) + diff(uy, y) + diff(uz, z)) = 0;

	globalExample = 10;
	geometryVolume geo = openVolume("../../Geometries/test3.msh");
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);
	//Stiffness Matrix

	int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 0,
		&Positions);
	printf(" Positions %d \n", sizePositions);

	double* stiffValues = getStiffnessValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		0, Positions, sizePositions);
	printf("stiffness Matrix \n");

	double* massValues = getMassValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		0, Positions, sizePositions);
	printf("Mass Matrix \n");


	double dt = 0.01;
	double heatCoef = 1.0 / 1.0;
	double coef = dt *heatCoef;

	double* Avalues = newDoubleV(sizePositions);

	for (int i = 0; i < sizePositions; i++)
	{
		Avalues[i] = massValues[i] + stiffValues[i] * coef;
	
	}

	double limit = 2.5;

	double t = 1;
	double* u0 = Heateq3D(geo.Nodes, geo.nodeSize, t);
	double* v = newDoubleV(geo.nodeSize);
	
	double* vant = multVectorSparse(massValues, Positions, sizePositions, u0, geo.nodeSize);
	int steps = (int)((limit - t) / dt);

	int k = 0;
	double** Values = newDouble(geo.nodeSize, steps + 1);

	//printVectorInt(geo.nodeSurface, geo.nodeSize, "Nodesurface.txt");

	printf("steps %d \n", steps);
	do
	{

		//Neuman Data
		int ntriangles = neumanTriangles(geo.triangleSurface, geo.triangleSize, 1);
		int** Neuman = newInt(ntriangles, 3);
		double** NeumanData = newDouble(ntriangles, 3);
		getNeumanData3DNEW(geo.Nodes, geo.Triangles, geo.nodeSize, geo.triangleSize,
			geo.triangleSurface, 1, NeumanData, Neuman, ntriangles, t);
		//Dirichlet Data
		double* nodeDirichlet = getDirichletData(geo.Nodes, geo.nodeSize, 0, geo.nodeSurface, t);
		//Load Data
		double* Load = Load3D(geo.Nodes, geo.Tetrahedras, geo.nodeSurface, geo.nodeSize,
			geo.tetrahedraSize, 0, nodeDirichlet, NeumanData, Neuman, ntriangles, t);

		double* coefF = multConstantV(Load, coef, geo.nodeSize);
		double* vantTemp = addVector(vant, coefF, geo.nodeSize);


		//double* v = PBCGSTABC(A, vantTemp, geo.nodeSize);
		double* v = PBCGSTABCSparse(Avalues, Positions, sizePositions, vantTemp, geo.nodeSize);

		for (int i = 0; i < geo.nodeSize; i++)
		{
			if (geo.nodeSurface[i] == 0)
			{
				Values[i][k] = nodeDirichlet[i];
			}
			else
			{
				Values[i][k] = v[i];
			}



		}

		double* Mv = multVectorSparse(massValues, Positions, sizePositions, v, geo.nodeSize);

		copyVector(vant, Mv, geo.nodeSize);

		t += dt;
		printf("t %.11f \n", t);

		k++;


		free(nodeDirichlet);
		free(Load);
		free(coefF);
		free(vantTemp);
		free(v);
		free(Mv);
		freeMatrixInt(Neuman, ntriangles);
		freeMatrix(NeumanData, ntriangles);

	} while (t < limit);



	printf("Label 1 ");

	double* exact = Heateq3D(geo.Nodes, geo.nodeSize, t - dt);



	double** Values2 = newDouble(geo.nodeSize, 2);
	double* result = newDoubleV(geo.nodeSize);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		Values2[i][0] = exact[i];
		Values2[i][1] = Values[i][k - 1];
		if (geo.nodeSurface[i] == 1)
		{
			result[i] = Values[i][k - 1];
		}
		else
		{
			exact[i] = 0;
		}

	}

	double* subResults = subVector(result, exact, geo.nodeSize);
	double normResults = norm(subResults, geo.nodeSize);
	printf("Norm %.10f \n", normResults);

	char** Names[2];
	Names[0] = "Exact";
	Names[1] = "Calculated";

	plot2DComparisonNamesNew(Values2, 2, geo.nodeSize, "Results", "Graph.html", Names, "Value");

	printVolumeValues("Results.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values2, 2, "Data");
}

void FEMTestTime4()
{
	//that is why it is a 1/3
	//diff(u, t) - 1 / 3(diff(ux, x) + diff(uy, y) + diff(uz, z)) = 0;

	globalExample = 20;
	geometryVolume geo = openVolume("../../Geometries/test3A.msh");
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);
	//Stiffness Matrix

	int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 0,
		&Positions);
	printf(" Positions %d \n", sizePositions);

	double* stiffValues = getStiffnessValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		0, Positions, sizePositions);
	printf("stiffness Matrix \n");

	double* massValues = getMassValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		0, Positions, sizePositions);
	printf("Mass Matrix \n");


	double dt = 0.01;
	double heatCoef = 1.0 / 3.0;
	double coef = dt * heatCoef;

	double* Avalues = newDoubleV(sizePositions);

	for (int i = 0; i < sizePositions; i++)
	{
		Avalues[i] = massValues[i] + stiffValues[i] * coef;

	}

	double limit = 2.5;

	double t = 1;
	double* u0 = Heateq3D(geo.Nodes, geo.nodeSize, t);
	double* v = newDoubleV(geo.nodeSize);

	double* vant = multVectorSparse(massValues, Positions, sizePositions, u0, geo.nodeSize);
	int steps = (int)((limit - t) / dt);

	int k = 0;
	double** Values = newDouble(geo.nodeSize, steps + 1);

	//printVectorInt(geo.nodeSurface, geo.nodeSize, "Nodesurface.txt");

	printf("steps %d \n", steps);
	do
	{

		//Neuman Data
		int ntriangles = neumanTriangles(geo.triangleSurface, geo.triangleSize, 1);
		int** Neuman = newInt(ntriangles, 3);
		double** NeumanData = newDouble(ntriangles, 3);
		getNeumanData3DNEW(geo.Nodes, geo.Triangles, geo.nodeSize, geo.triangleSize,
			geo.triangleSurface, 1, NeumanData, Neuman, ntriangles, t);
		//Dirichlet Data
		double* nodeDirichlet = getDirichletData(geo.Nodes, geo.nodeSize, 0, geo.nodeSurface, t);
		//Load Data
		double* Load = Load3D(geo.Nodes, geo.Tetrahedras, geo.nodeSurface, geo.nodeSize,
			geo.tetrahedraSize, 0, nodeDirichlet, NeumanData, Neuman, ntriangles, t);

		double* coefF = multConstantV(Load, coef, geo.nodeSize);
		double* vantTemp = addVector(vant, coefF, geo.nodeSize);


		//double* v = PBCGSTABC(A, vantTemp, geo.nodeSize);
		double* v = PBCGSTABCSparse(Avalues, Positions, sizePositions, vantTemp, geo.nodeSize);

		for (int i = 0; i < geo.nodeSize; i++)
		{
			if (geo.nodeSurface[i] == 0)
			{
				Values[i][k] = nodeDirichlet[i];
			}
			else
			{
				Values[i][k] = v[i];
			}



		}

		double* Mv = multVectorSparse(massValues, Positions, sizePositions, v, geo.nodeSize);

		copyVector(vant, Mv, geo.nodeSize);

		t += dt;
		printf("t %.11f \n", t);

		k++;


		free(nodeDirichlet);
		free(Load);
		free(coefF);
		free(vantTemp);
		free(v);
		free(Mv);
		freeMatrixInt(Neuman, ntriangles);
		freeMatrix(NeumanData, ntriangles);

	} while (t < limit);



	printf("Label 1 ");

	double* exact = Heateq3D(geo.Nodes, geo.nodeSize, t - dt);



	double** Values2 = newDouble(geo.nodeSize, 2);
	double* result = newDoubleV(geo.nodeSize);
	for (int i = 0; i < geo.nodeSize; i++)
	{
		Values2[i][0] = exact[i];
		Values2[i][1] = Values[i][k - 1];
		if (geo.nodeSurface[i] == 1)
		{
			result[i] = Values[i][k - 1];
		}
		else
		{
			exact[i] = 0;
		}

	}

	double* subResults = subVector(result, exact, geo.nodeSize);
	double normResults = norm(subResults, geo.nodeSize);
	printf("Norm %.10f \n", normResults);

	char** Names[2];
	Names[0] = "Exact";
	Names[1] = "Calculated";

	plot2DComparisonNamesNew(Values2, 2, geo.nodeSize, "Results", "Graph.html", Names, "Value");

	printVolumeValues("Results.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values2, 2, "Data");
}