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
#include "../Amazigh/openFile.h"
#include "cell.h"
#ifdef  _WIN32
#define OSVAR 0
#else
#define OSVAR 1 
#endif
void CreatePreconditioner(double* AV, int** APositions, int Asize, double* CondV)
{
	
	for (int i = 0; i < Asize; i++)
	{
		if (APositions[i][ 0] == APositions[i][ 1] && AV[i] != 0)
		{
			CondV[i] = 1.0 / AV[i];
		}
	}
}

void CreateCrankNicholsonMatrix(double* MassMatrixV, 
	double* StiffMatrixV, int** Positions, int sizePositions, double coef,
	double omega, double* AV, double* BV)
{
	//A[i, j] = Mass[i, j] + Stiff[i, j] * coef * omega;
	for (int i = 0; i < sizePositions; i++)
	{
		AV[i] = MassMatrixV[i] + StiffMatrixV[i] * coef * omega;
	}

	//B[i, j] = Mass[i, j] - Stiff[i, j] * coef * (1 - omega);
	
	for (int i = 0; i < sizePositions; i++)
	{
		BV[i] = MassMatrixV[i] -StiffMatrixV[i] * coef * (1 - omega);
	}
}


void monodomainHeart()
{
	geometryVolume geo = openVolume("../Geometries/heartVol.msh");
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);

	double** fibers = openMatrix("../Geometries/fibersHeart.matrix");

	double beta = 400;
	double dt = 0.25;
	double phi = 0.8;
	double coef = (dt * (phi / (phi + 1)));
	double omega = 1.0 / 2.0;

	double cli = 0.003 / beta;
	double cti = 0.001 / beta;
	double cni = 0.00031525 / beta;

	double cle = 0.002 / beta;
	double cte = 0.00165 / beta;
	double cne = 0.00013514 / beta;

	int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 5,
		&Positions);
	printf(" Positions %d \n", sizePositions);

	double* stiffValues = getStiffnessValuesTensor(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli, cti, cni,fibers);
	printf("stiffness Matrix \n");

	double* massValues = getMassValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions);
	printf("Mass Matrix \n");

	double* MVtemp = getStiffnessValuesTensor(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli, cti, cni, fibers);
	printf("MV Matrix \n");

	double* MV = multConstantV(MVtemp, -1, sizePositions);

	double* NV = getStiffnessValuesTensor(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli + cle, cti + cte, cni + cne, fibers);
	printf("NV Matrix \n");

	double* AV = newDoubleV(sizePositions);
	double* BV = newDoubleV(sizePositions);

	CreateCrankNicholsonMatrix(massValues,
		stiffValues, Positions, sizePositions, coef,
		omega, AV, BV);
	printf("Create AV BV \n");

	double* u0 = newDoubleV(geo.nodeSize);

	double* vant = multVectorSparse(BV, Positions, sizePositions, u0, geo.nodeSize);

	double limit = 1000;
	double t = 0;
	double steps = (limit - t) / dt;

	int stepsint = (int) steps + 1;
	int k = 0;

	double** Values = newDouble(geo.nodeSize, stepsint);
	double** Valuesue = newDouble(geo.nodeSize, stepsint);

	double** ALGEBRAICTotal = newDouble(geo.nodeSize, 8);
	double** CONSTANTSTotal = newDouble(geo.nodeSize, 22);
	double** RATESTotal = newDouble(geo.nodeSize, 3);
	double** STATESTotal = newDouble(geo.nodeSize, 3);
	
	for (int i = 0; i < geo.nodeSize; i++)
	{
		initConsts_FemtomKarmaMod(CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i]);
	}
	printf("initConsts_FemtomKarmaMod \n");

	do
	{
		double* FuncVnm1 = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			STATESTotal[i][0] = u0[i];
			computeRates_FemtomKarmaMod(t, CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i], ALGEBRAICTotal[i], 0);
			eulerODE(RATESTotal[i], STATESTotal[i], 3, dt);

			FuncVnm1[i] = STATESTotal[i][0];
		}
		double* FuncCurrent = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			double current = 0;

			if (i == 887 && t < 1)
			{
				current = (-2);
			}
			if (i == 790 && t < 1)
			{
				current = (-2);
			}
			if (i == 197 && t < 2)
			{
				current = (-2);
			}
			FuncCurrent[i] = -dt * current;
		}

		double* BFuncVnm1 = multVectorSparse(BV, Positions, sizePositions, FuncVnm1, geo.nodeSize);
		copyVector(FuncVnm1, BFuncVnm1,geo.nodeSize);
		
		double* BFuncCurrent = multVectorSparse(BV, Positions, sizePositions, FuncCurrent, geo.nodeSize);
		copyVector(FuncCurrent, BFuncCurrent, geo.nodeSize);

		double* FuncVnm1plusCurrent = addVector(FuncVnm1, FuncCurrent, geo.nodeSize);
		copyVector(vant, FuncVnm1plusCurrent, geo.nodeSize);

		double* v = PBCGSTABCSparse(AV, Positions, sizePositions, vant, geo.nodeSize);

		double* vm = multVectorSparse(MV, Positions, sizePositions, v, geo.nodeSize);

		double* ue= PBCGSTABCSparse(NV, Positions, sizePositions, vm, geo.nodeSize);

		
		copyVector(u0, v,geo.nodeSize);

		for (int i = 0; i < geo.nodeSize; i++)
		{
			Values[i][ k] = v[i];
			Valuesue[i][k] = ue[i];
		}
		
		printf(" t %0.4f \n", t);
		double maxValue = maxVector(v, geo.nodeSize);
		printf(" maxValue %0.8f \n", maxValue);
		
		
		
		t += dt;
		k++;

		//end
		free(FuncVnm1);
		free(FuncCurrent);
		free(BFuncVnm1);
		free(BFuncCurrent);
		free(FuncVnm1plusCurrent);
		free(v);
		free(vm);
		free(ue);
	}
	while (t < limit);

	printVolumeValues("Transmembrane.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values, stepsint, "Transmembrane Potential");

	printVolumeValues("Extracellular.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Valuesue, stepsint, "Extracellular Potential");
}

void monodomainBrain()
{
	char* geometryName;
	char* fiberName;
	char* positionName;
	if (OSVAR == 0)
	{
		geometryName = "../../Geometries/volBrain.msh";
		fiberName = "../../Geometries/fibersBrain.matrix";
		positionName = "../../Geometries/positionsBrain.txt";
	}
	else
	{
		geometryName = "../Geometries/volBrain.msh";
		fiberName = "../Geometries/fibersBrain.matrix";
		positionName = "../Geometries/positionsBrain.txt";

	}

	geometryVolume geo = openVolume(geometryName);
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);

	double** fibers = openMatrix(fiberName);

	 double beta_brain = 0.08;
	 double dt = 0.5;
	 double phi = 0.8;
	 double coef = (dt * (phi / (phi + 1)));
	 double omega = 1.0 / 2.0;

	//double cli = 0.003 / beta;
	 double cli_brain = 0.30 / beta_brain;
	 double cti_brain = 0.01 / beta_brain;
	 double cni_brain = 0.0031525 / beta_brain;

	//double cle = 0.002 / beta;
	 double cle_brain = 0.2 / beta_brain;
	 double cte_brain = 0.0165 / beta_brain;
	 double cne_brain = 0.0013514 / beta_brain;

	 double beta_muscle = 0.4;
	 double cli_muscle = 0.30 / beta_muscle;
	 double cti_muscle = 0.01 / beta_muscle;
	 double cni_muscle = 0.0031525 / beta_muscle;

	//double cle = 0.002 / beta;
	 double cle_muscle = 0.2 / beta_muscle;
	 double cte_muscle = 0.0165 / beta_muscle;
	 double cne_muscle = 0.0013514 / beta_muscle;

	int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 5,
		&Positions);
	printf(" Positions %d \n", sizePositions);

	double* stiffValues = getStiffnessValuesTensor(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli_brain, cti_brain, cni_brain, fibers);
	printf("stiffness Matrix \n");

	double* massValues = getMassValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions);
	printf("Mass Matrix \n");

	double* MVtemp = getStiffnessValuesTensor(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli_brain, cti_brain, cni_brain, fibers);
	printf("MV Matrix \n");

	double* MV = multConstantV(MVtemp, -1, sizePositions);

	double* NV = getStiffnessValuesTensor(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli_brain + cle_brain, cti_brain + cte_brain, 
		cni_brain + cne_brain, fibers);
	printf("NV Matrix \n");

	double* AV = newDoubleV(sizePositions);
	double* BV = newDoubleV(sizePositions);

	CreateCrankNicholsonMatrix(massValues,
		stiffValues, Positions, sizePositions, coef,
		omega, AV, BV);
	printf("Create AV BV \n");

	double* u0 = newDoubleV(geo.nodeSize);

	double* vant = multVectorSparse(BV, Positions, sizePositions, u0, geo.nodeSize);

	double limit = 150;
	double t = 0;
	double steps = (limit - t) / dt;

	int stepsint = (int)steps + 1;
	int k = 0;

	double** Values = newDouble(geo.nodeSize, stepsint);
	double** Valuesue = newDouble(geo.nodeSize, stepsint);

	double** ALGEBRAICTotal = newDouble(geo.nodeSize, 8);
	double** CONSTANTSTotal = newDouble(geo.nodeSize, 22);
	double** RATESTotal = newDouble(geo.nodeSize, 3);
	double** STATESTotal = newDouble(geo.nodeSize, 3);

	for (int i = 0; i < geo.nodeSize; i++)
	{
		initConsts_FemtomKarmaMod(CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i]);
	}
	printf("initConsts_FemtomKarmaMod \n");

	do
	{
		double* FuncVnm1 = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			STATESTotal[i][0] = u0[i];
			computeRates_FemtomKarmaMod(t, CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i], ALGEBRAICTotal[i], 0);
			eulerODE(RATESTotal[i], STATESTotal[i], 3, dt);

			FuncVnm1[i] = STATESTotal[i][0];
		}
		double* FuncCurrent = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			double current_brain = 0;

			if (i == 3049 && t < 1)
				//if (i == 1411 && t < 1)
			{
				current_brain = (-5);

			}
			FuncCurrent[i] = -dt * current_brain;
		}

		double* BFuncVnm1 = multVectorSparse(BV, Positions, sizePositions, FuncVnm1, geo.nodeSize);
		copyVector(FuncVnm1, BFuncVnm1, geo.nodeSize);

		double* BFuncCurrent = multVectorSparse(BV, Positions, sizePositions, FuncCurrent, geo.nodeSize);
		copyVector(FuncCurrent, BFuncCurrent, geo.nodeSize);

		double* FuncVnm1plusCurrent = addVector(FuncVnm1, FuncCurrent, geo.nodeSize);
		copyVector(vant, FuncVnm1plusCurrent, geo.nodeSize);

		double* v = PBCGSTABCSparse(AV, Positions, sizePositions, vant, geo.nodeSize);

		double* vm = multVectorSparse(MV, Positions, sizePositions, v, geo.nodeSize);

		double* ue = PBCGSTABCSparse(NV, Positions, sizePositions, vm, geo.nodeSize);


		copyVector(u0, v, geo.nodeSize);

		for (int i = 0; i < geo.nodeSize; i++)
		{
			Values[i][k] = v[i];
			Valuesue[i][k] = ue[i];
		}

		printf(" t %0.4f \n", t);
		double maxValue = maxVector(v, geo.nodeSize);
		printf(" maxValue %0.8f \n", maxValue);



		t += dt;
		k++;

		//end
		free(FuncVnm1);
		free(FuncCurrent);
		free(BFuncVnm1);
		free(BFuncCurrent);
		free(FuncVnm1plusCurrent);
		free(v);
		free(vm);
		free(ue);
	} while (t < limit);

	printVolumeValues("Transmembrane.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values, stepsint, "Transmembrane Potential Brain");

	/*printVolumeValues("Extracellular.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Valuesue, stepsint, "Extracellular Potential Brain");*/
}

void monodomainBrainPatientBackUp()
{
    printf("test \n");
	geometryVolume geo = openVolume("../Geometries/volBrain.msh");
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);

	double** fibers = openMatrix("../Geometries/fibersBrain.matrix");

	double beta_brain = 0.08;
	double dt = 0.5;
	double phi = 0.8;
	double coef = (dt * (phi / (phi + 1)));
	double omega = 1.0 / 2.0;

	//double cli = 0.003 / beta;
	double cli_brain = 0.30 / beta_brain;
	double cti_brain = 0.01 / beta_brain;
	double cni_brain = 0.0031525 / beta_brain;

	//double cle = 0.002 / beta;
	double cle_brain = 0.2 / beta_brain;
	double cte_brain = 0.0165 / beta_brain;
	double cne_brain = 0.0013514 / beta_brain;

	double beta_muscle = 0.04;
	double cli_muscle = 0.30 / beta_muscle;
	double cti_muscle = 0.01 / beta_muscle;
	double cni_muscle = 0.0031525 / beta_muscle;

	//double cle = 0.002 / beta;
	double cle_muscle = 0.2 / beta_muscle;
	double cte_muscle = 0.0165 / beta_muscle;
	double cne_muscle = 0.0013514 / beta_muscle;

	/*int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 5,
		&Positions);
	printf(" Positions %d \n", sizePositions)*/;

	int** Positions= openMatrixInt("../Positions.txt");
	int sizePositions = 78635;
	printf(" Positions %d \n", sizePositions);

	double* massValues = getMassValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions);
	printf("Mass Matrix \n");

	//printMatrixInt(Positions, sizePositions, 2, "Positions.txt");
	printVector(massValues, sizePositions, "massValues.txt");

	double* stiffValues = getStiffnessValuesTensorLesion(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli_brain, cti_brain, cni_brain, fibers,1599,0.126);
	printf("stiffness Matrix \n");

	

	

	

	double* MVtemp = getStiffnessValuesTensorLesion(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli_brain, cti_brain, cni_brain, fibers,1599, 0.126);
	printf("MV Matrix \n");

	double* MV = multConstantV(MVtemp, -1, sizePositions);

	double* NV = getStiffnessValuesTensorLesion(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli_brain + cle_brain, cti_brain + cte_brain,
		cni_brain + cne_brain, fibers, 1599, 0.126);
	printf("NV Matrix \n");

	double* AV = newDoubleV(sizePositions);
	double* BV = newDoubleV(sizePositions);

	CreateCrankNicholsonMatrix(massValues,
		stiffValues, Positions, sizePositions, coef,
		omega, AV, BV);
	printf("Create AV BV \n");

	double* u0 = newDoubleV(geo.nodeSize);

	double* vant = multVectorSparse(BV, Positions, sizePositions, u0, geo.nodeSize);

	double limit = 150;
	double t = 0;
	double steps = (limit - t) / dt;

	int stepsint = (int)steps + 1;
	int k = 0;

	double** Values = newDouble(geo.nodeSize, stepsint);
	double** Valuesue = newDouble(geo.nodeSize, stepsint);

	double** ALGEBRAICTotal = newDouble(geo.nodeSize, 8);
	double** CONSTANTSTotal = newDouble(geo.nodeSize, 22);
	double** RATESTotal = newDouble(geo.nodeSize, 3);
	double** STATESTotal = newDouble(geo.nodeSize, 3);

	for (int i = 0; i < geo.nodeSize; i++)
	{
		initConsts_FemtomKarmaMod(CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i]);
	}
	printf("initConsts_FemtomKarmaMod \n");

	do
	{
		double* FuncVnm1 = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			STATESTotal[i][0] = u0[i];
			computeRates_FemtomKarmaMod(t, CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i], ALGEBRAICTotal[i], 0);
			eulerODE(RATESTotal[i], STATESTotal[i], 3, dt);

			FuncVnm1[i] = STATESTotal[i][0];
		}
		double* FuncCurrent = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			double current_brain = 0;

			if (i == 3049 && t < 1)
				//if (i == 1411 && t < 1)
			{
				current_brain = (-5);

			}
			FuncCurrent[i] = -dt * current_brain;
		}

		double* BFuncVnm1 = multVectorSparse(BV, Positions, sizePositions, FuncVnm1, geo.nodeSize);
		copyVector(FuncVnm1, BFuncVnm1, geo.nodeSize);

		double* BFuncCurrent = multVectorSparse(BV, Positions, sizePositions, FuncCurrent, geo.nodeSize);
		copyVector(FuncCurrent, BFuncCurrent, geo.nodeSize);

		double* FuncVnm1plusCurrent = addVector(FuncVnm1, FuncCurrent, geo.nodeSize);
		copyVector(vant, FuncVnm1plusCurrent, geo.nodeSize);

		double* v = PBCGSTABCSparse(AV, Positions, sizePositions, vant, geo.nodeSize);

		double* vm = multVectorSparse(MV, Positions, sizePositions, v, geo.nodeSize);

		double* ue = PBCGSTABCSparse(NV, Positions, sizePositions, vm, geo.nodeSize);


		copyVector(u0, v, geo.nodeSize);

		for (int i = 0; i < geo.nodeSize; i++)
		{
			Values[i][k] = v[i];
			Valuesue[i][k] = ue[i];
		}

		printf(" t %0.4f \n", t);
		double maxValue = maxVector(v, geo.nodeSize);
		printf(" maxValue %0.8f \n", maxValue);



		t += dt;
		k++;

		//end
		free(FuncVnm1);
		free(FuncCurrent);
		free(BFuncVnm1);
		free(BFuncCurrent);
		free(FuncVnm1plusCurrent);
		free(v);
		free(vm);
		free(ue);
	} while (t < limit);

	printVolumeValues("Transmembrane.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values, stepsint, "Transmembrane Potential Brain Patient");

	/*printVolumeValues("Extracellular.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Valuesue, stepsint, "Extracellular Potential Brain");*/
}

void monodomainBrainPatient()
{
	char* geometryName;
	char* fiberName;
	char* positionName;
	if (OSVAR == 0)
	{
		geometryName = "../../Geometries/volBrain.msh";
		fiberName = "../../Geometries/fibersBrain.matrix";
		positionName = "../../Geometries/positionsBrain.txt";
	}
	else
	{
		geometryName = "../Geometries/volBrain.msh";
		fiberName = "../Geometries/fibersBrain.matrix";
		positionName = "../Geometries/positionsBrain.txt";
	
	}

	geometryVolume geo = openVolume(geometryName);
	printf(" geo.tetrahedraSize %d \n", geo.tetrahedraSize);

	double** fibers = openMatrix(fiberName);

	double beta_brain = 0.08;
	double dt = 0.5;
	double phi = 0.8;
	double coef = (dt * (phi / (phi + 1)));
	double omega = 1.0 / 2.0;

	//double cli = 0.003 / beta;
	double cli_brain = 0.30 / beta_brain;
	double cti_brain = 0.01 / beta_brain;
	double cni_brain = 0.0031525 / beta_brain;

	//double cle = 0.002 / beta;
	double cle_brain = 0.2 / beta_brain;
	double cte_brain = 0.0165 / beta_brain;
	double cne_brain = 0.0013514 / beta_brain;

	/*int** Positions;
	int sizePositions = getPositions(geo.Nodes, geo.Triangles, geo.Tetrahedras,
	geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize, 5,
	&Positions);
	printf(" Positions %d \n", sizePositions)*/;

	int** Positions = openMatrixInt(positionName);
	int sizePositions = 78635;
	printf(" Positions %d \n", sizePositions);

	double* massValues = getMassValues(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions);
	printf("Mass Matrix \n");

	//printMatrixInt(Positions, sizePositions, 2, "Positions.txt");

	double* stiffValues = getStiffnessValuesTensorLesion(geo.Nodes, geo.Triangles, geo.Tetrahedras,
		geo.nodeSurface, geo.nodeSize, geo.triangleSize, geo.tetrahedraSize,
		5, Positions, sizePositions, cli_brain, cti_brain, cni_brain, fibers, 1599, 0.126);
	printf("stiffness Matrix \n");

	double* AV = newDoubleV(sizePositions);
	double* BV = newDoubleV(sizePositions);

	CreateCrankNicholsonMatrix(massValues,
		stiffValues, Positions, sizePositions, coef,
		omega, AV, BV);
	printf("Create AV BV \n");

	double* u0 = newDoubleV(geo.nodeSize);

	double* vant = multVectorSparse(BV, Positions, sizePositions, u0, geo.nodeSize);

	double limit =150.0;
	double t = 0;
	double steps = (limit - t) / dt;

	int stepsint = (int)steps + 1;
	int k = 0;

	double** Values = newDouble(geo.nodeSize, stepsint);

	double** ALGEBRAICTotal = newDouble(geo.nodeSize, 8);
	double** CONSTANTSTotal = newDouble(geo.nodeSize, 22);
	double** RATESTotal = newDouble(geo.nodeSize, 3);
	double** STATESTotal = newDouble(geo.nodeSize, 3);

	for (int i = 0; i < geo.nodeSize; i++)
	{
		initConsts_FemtomKarmaMod(CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i]);
	}
	printf("initConsts_FemtomKarmaMod \n");

	do
	{
		double* FuncVnm1 = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			STATESTotal[i][0] = u0[i];
			computeRates_FemtomKarmaMod(t, CONSTANTSTotal[i], RATESTotal[i], STATESTotal[i], ALGEBRAICTotal[i], 0);
			eulerODE(RATESTotal[i], STATESTotal[i], 3, dt);

			FuncVnm1[i] = STATESTotal[i][0];
		}
		double* FuncCurrent = newDoubleV(geo.nodeSize);
		for (int i = 0; i < geo.nodeSize; i++)
		{
			double current_brain = 0;

			if (i == 3049 && t < 1)
				//if (i == 1411 && t < 1)
			{
				current_brain = (-5);

			}
			FuncCurrent[i] = -dt * current_brain;
		}


		double* BFuncVnm1 = multVectorSparse(BV, Positions, sizePositions, FuncVnm1, geo.nodeSize);
		copyVector(FuncVnm1, BFuncVnm1, geo.nodeSize);

		double* BFuncCurrent = multVectorSparse(BV, Positions, sizePositions, FuncCurrent, geo.nodeSize);
		copyVector(FuncCurrent, BFuncCurrent, geo.nodeSize);

		double* FuncVnm1plusCurrent = addVector(FuncVnm1, FuncCurrent, geo.nodeSize);
		copyVector(vant, FuncVnm1plusCurrent, geo.nodeSize);

		//double* v = PBCGSTABCSparse(AV, Positions, sizePositions, vant, geo.nodeSize);

		double* v = SolveSparsePCG(Positions, AV, sizePositions, vant,
			u0, geo.nodeSize, 1e-14);

		copyVector(u0, v, geo.nodeSize);

		for (int i = 0; i < geo.nodeSize; i++)
		{
			Values[i][k] = v[i];
		}

		printf(" t %0.4f \n", t);
		double maxValue = maxVector(v, geo.nodeSize);
		printf(" maxValue %0.8f \n", maxValue);



		t += dt;
		k++;

		//end
		free(FuncVnm1);
		free(FuncCurrent);
		free(BFuncVnm1);
		free(BFuncCurrent);
		free(FuncVnm1plusCurrent);
		free(v);

	} while (t < limit);

	printMatrix(Values, geo.nodeSize, stepsint, "Transmembrane.txt");

	printVolumeValues("Transmembrane.msh", geo.Nodes, geo.nodeSize, geo.nodeSurface,
		geo.Triangles, geo.triangleSize, geo.triangleSurface,
		geo.Tetrahedras, geo.tetrahedraSize, geo.tetrahedraSurface, Values, stepsint, "Transmembrane Potential Brain Patient");

}
