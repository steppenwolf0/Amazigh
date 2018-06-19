#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "geometry.h"
#include "../Amazigh/supportC.h"
#include "../Amazigh/printFile.h"

void freeGeometrySurface(geometrySurface surface)
{
	free(surface.nodeSurface);
	freeMatrix(surface.Nodes, surface.nodeSize);
	free(surface.triangleSurface);
	freeMatrixInt(surface.Triangles, surface.triangleSize);

}

geometrySurface openSurface(char* fileName)
{
	//Triangles are move to base 0
	FILE* file = fopen(fileName, "r"); /* should check the result */
	int n = 100;
	char line[100];
	for(int i = 0; i < 4;i++)
	{
		printf("%s", fgets(line, sizeof(line), file));
	}
	int nodeSize = atoi(fgets(line, sizeof(line), file));
	printf("%d\n", nodeSize);
	double** Nodes = newDouble(nodeSize, 3);
	char line2[1000];
	for (int i = 0; i<nodeSize; i++)
	{
		fgets(line2, sizeof(line2), file);
		char *p;
		p = strtok(line2, " ");
		for (int j = 0; j<4; j++)
		{
			if (p)
				if (j!=0)
					Nodes[i][j-1] = atof(p);
			p = strtok(NULL, " ");
		}
	}
	for (int i = 0; i < 2; i++)
	{
		printf("%s", fgets(line, sizeof(line), file));
	}
	int triangleSize = atoi(fgets(line, sizeof(line), file));
	printf("%d\n", triangleSize);
	int** Triangles = newInt(triangleSize, 3);
	int* triangleSurface = newIntV(triangleSize);

	for (int i = 0; i<triangleSize; i++)
	{
		fgets(line2, sizeof(line2), file);
		char *p;
		p = strtok(line2, " ");
		for (int j = 0; j<8; j++)
		{
			if (p)
			{
				if (j == 4)
					triangleSurface[i] = atoi(p);
				if (j >4)
					Triangles[i][j - 5] = atoi(p) - 1;
			}
				
			p = strtok(NULL, " ");
		}
	}
	
	int* nodeSurface = newIntV(nodeSize);

	for (int i = 0; i < triangleSize; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			nodeSurface[Triangles[i][j]] = triangleSurface[i];
		}
	}
		


	struct geometrySurface outGeometry;
	outGeometry.nodeSize = nodeSize;
	outGeometry.Nodes = Nodes;
	outGeometry.nodeSurface = nodeSurface;
	
	outGeometry.triangleSize = triangleSize;
	outGeometry.Triangles = Triangles;
	outGeometry.triangleSurface = triangleSurface;
	
	return outGeometry;
}

void printSurface(char* fileName, double** Nodes, int nodeSize, int* nodeSurface,
	int** Triangles, int triangleSize, int* triangleSurface)
{
	FILE *f = fopen(fileName, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f, "$MeshFormat\n");
	fprintf(f, "2.2 0 8\n");
	fprintf(f, "$EndMeshFormat\n");
	fprintf(f, "$Nodes\n");
	fprintf(f, "%d\n", nodeSize);
	for (int i = 0; i < nodeSize; i++)
	{
		fprintf(f, "%d ",i+1);
		for (int j = 0; j < 3;j++)
			fprintf(f, "%.18f ", Nodes[i][j]);
		fprintf(f, "\n");
	}
	fprintf(f, "$EndNodes\n");
	fprintf(f, "$Elements\n");
	fprintf(f, "%d\n", triangleSize);
	for (int i = 0; i < triangleSize; i++)
	{
		fprintf(f, "%d ", i + 1);
		fprintf(f, "2 2 0 ");
		fprintf(f, "%d ", triangleSurface[i]);
		for (int j = 0; j < 3; j++)
			fprintf(f, "%d ", Triangles[i][j]+1);
		fprintf(f, "\n");
	}
	fprintf(f, "$EndElements\n");
	fclose(f);
}

geometryVolume openVolume(char* fileName)
{
	//Triangles are move to base 0
	FILE* file = fopen(fileName, "r"); /* should check the result */
	int n = 100;
	char line[100];
	for (int i = 0; i < 4; i++)
	{
		printf("%s", fgets(line, sizeof(line), file));
	}
	int nodeSize = atoi(fgets(line, sizeof(line), file));
	printf("%d\n", nodeSize);
	double** Nodes = newDouble(nodeSize, 3);
	char line2[1000];
	for (int i = 0; i<nodeSize; i++)
	{
		fgets(line2, sizeof(line2), file);
		char *p;
		p = strtok(line2, " ");
		for (int j = 0; j<4; j++)
		{
			if (p)
				if (j != 0)
					Nodes[i][j - 1] = atof(p);
			p = strtok(NULL, " ");
		}
	}
	for (int i = 0; i < 2; i++)
	{
		printf("%s", fgets(line, sizeof(line), file));
	}

	int elementSize = atoi(fgets(line, sizeof(line), file));
	printf("%d\n", elementSize);
	int** Elements = newInt(elementSize, 4);
	int* elementSurface = newIntV(elementSize);
	int* elementType = newIntV(elementSize);
	
	int triangleSize = 0;
	int tetahedraSize = 0;
	for (int i = 0; i<elementSize; i++)
	{
		fgets(line2, sizeof(line2), file);
		char *p;
		p = strtok(line2, " ");
		for (int j = 0; j<9; j++)
		{
			if (p)
			{
				if (j == 1)
				{
					elementType[i] = atoi(p);
					if (elementType[i] == 2)
						triangleSize++;
					if (elementType[i] == 4)
						tetahedraSize++;
				}
					
				if (j == 4)
					elementSurface[i] = atoi(p);
				if (j >4)
					Elements[i][j - 5] = atoi(p) - 1;
			}
			
			p = strtok(NULL, " ");
		}
	}

	int** Triangles = newInt(triangleSize, 3);
	int** Tetrahedras = newInt(tetahedraSize, 4);
	int* triangleSurface = newIntV(triangleSize);
	int* tetrahedraSurface = newIntV(tetahedraSize);
	int trianglesIndex = 0;
	int tetrahedrasIndex = 0;
	for (int i = 0; i < elementSize; i++)
	{
		if (elementType[i] == 2)
		{
			triangleSurface[trianglesIndex] = elementSurface[i];
			for (int j = 0; j < 3; j++)
				Triangles[trianglesIndex][j] = Elements[i][j];
			trianglesIndex++;
			
		}
		if (elementType[i] == 4)
		{
			tetrahedraSurface[tetrahedrasIndex] = elementSurface[i];
			for (int j = 0; j < 4; j++)
				Tetrahedras[tetrahedrasIndex][j] = Elements[i][j];
			tetrahedrasIndex++;
		}
	}

	/*printMatrixInt(Triangles, triangleSize, 3,"Triangles.txt");
	printMatrixInt(Tetrahedras, tetahedraSize, 4, "Tetrahedras.txt");*/
	int* nodeSurface = newIntV(nodeSize);
	for (int i = 0; i < tetahedraSize; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			nodeSurface[Tetrahedras[i][j]] = -tetrahedraSurface[i];
		}
	}
	for (int i = 0; i < triangleSize; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			nodeSurface[Triangles[i][j]] = triangleSurface[i];
		}
	}

	freeMatrixInt(Elements, elementSize);
	free(elementSurface);
	free(elementType);

	struct geometryVolume outGeometry;
	outGeometry.nodeSize = nodeSize;
	outGeometry.Nodes = Nodes;
	outGeometry.nodeSurface = nodeSurface;

	outGeometry.triangleSize = triangleSize;
	outGeometry.Triangles = Triangles;
	outGeometry.triangleSurface = triangleSurface;

	outGeometry.tetrahedraSize = tetahedraSize;
	outGeometry.Tetrahedras = Tetrahedras;
	outGeometry.tetrahedraSurface = tetrahedraSurface;

	return outGeometry;
}

void printVolume(char* fileName, double** Nodes, int nodeSize, int* nodeSurface,
	int** Triangles, int triangleSize, int* triangleSurface,
	int** Tetrahedras, int tetrahedraSize, int* tetrahedraSurface)
{
	FILE *f = fopen(fileName, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f, "$MeshFormat\n");
	fprintf(f, "2.2 0 8\n");
	fprintf(f, "$EndMeshFormat\n");
	fprintf(f, "$Nodes\n");
	fprintf(f, "%d\n", nodeSize);
	for (int i = 0; i < nodeSize; i++)
	{
		fprintf(f, "%d ", i + 1);
		for (int j = 0; j < 3; j++)
			fprintf(f, "%.18f ", Nodes[i][j]);
		fprintf(f, "\n");
	}
	fprintf(f, "$EndNodes\n");
	fprintf(f, "$Elements\n");
	fprintf(f, "%d\n", triangleSize+tetrahedraSize);
	for (int i = 0; i < triangleSize; i++)
	{
		fprintf(f, "%d ", i + 1);
		fprintf(f, "2 2 0 ");
		fprintf(f, "%d ", triangleSurface[i]);
		for (int j = 0; j < 3; j++)
			fprintf(f, "%d ", Triangles[i][j] + 1);
		fprintf(f, "\n");
	}
	for (int i = 0; i < tetrahedraSize; i++)
	{
		fprintf(f, "%d ", i + 1+ triangleSize);
		fprintf(f, "4 2 0 ");
		fprintf(f, "%d ", tetrahedraSurface[i]);
		for (int j = 0; j < 4; j++)
			fprintf(f, "%d ", Tetrahedras[i][j] + 1);
		fprintf(f, "\n");
	}
	fprintf(f, "$EndElements\n");
	fclose(f);
}

void printVolumeValues(char* fileName, double** Nodes, int nodeSize, int* nodeSurface,
	int** Triangles, int triangleSize, int* triangleSurface,
	int** Tetrahedras, int tetrahedraSize, int* tetrahedraSurface, double** Values, int valuesSize, char* dataName)
{
	FILE *f = fopen(fileName, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f, "$MeshFormat\n");
	fprintf(f, "2.2 0 8\n");
	fprintf(f, "$EndMeshFormat\n");
	fprintf(f, "$Nodes\n");
	fprintf(f, "%d\n", nodeSize);
	for (int i = 0; i < nodeSize; i++)
	{
		fprintf(f, "%d ", i + 1);
		for (int j = 0; j < 3; j++)
			fprintf(f, "%.18f ", Nodes[i][j]);
		fprintf(f, "\n");
	}
	fprintf(f, "$EndNodes\n");
	fprintf(f, "$Elements\n");
	fprintf(f, "%d\n", triangleSize + tetrahedraSize);
	for (int i = 0; i < triangleSize; i++)
	{
		fprintf(f, "%d ", i + 1);
		fprintf(f, "2 2 0 ");
		fprintf(f, "%d ", triangleSurface[i]);
		for (int j = 0; j < 3; j++)
			fprintf(f, "%d ", Triangles[i][j] + 1);
		fprintf(f, "\n");
	}
	for (int i = 0; i < tetrahedraSize; i++)
	{
		fprintf(f, "%d ", i + 1 + triangleSize);
		fprintf(f, "4 2 0 ");
		fprintf(f, "%d ", tetrahedraSurface[i]);
		for (int j = 0; j < 4; j++)
			fprintf(f, "%d ", Tetrahedras[i][j] + 1);
		fprintf(f, "\n");
	}
	fprintf(f, "$EndElements\n");
	for (int j = 0; j < valuesSize; j++)
	{
		fprintf(f,"$NodeData\n");
		fprintf(f, "1\n");
		fprintf(f,"%s\n", dataName);
		fprintf(f, "1\n");
		fprintf(f, "%d\n",j);
		fprintf(f, "3\n");
		fprintf(f, "%d\n", j);
		fprintf(f, "1\n");
		fprintf(f, "%d\n", nodeSize);

		for (int i = 0; i < nodeSize; i++)
		{
			fprintf(f, "%d %0.10f \n", (i + 1), Values[i][j]);

		}
		fprintf(f, "$EndNodeData\n");
	}

	fclose(f);
}

void freeGeometryVolume(geometryVolume volume)
{
	free(volume.nodeSurface);
	freeMatrix(volume.Nodes, volume.nodeSize);
	free(volume.triangleSurface);
	freeMatrixInt(volume.Triangles, volume.triangleSize);
	free(volume.tetrahedraSurface);
	freeMatrixInt(volume.Tetrahedras, volume.tetrahedraSize);

}