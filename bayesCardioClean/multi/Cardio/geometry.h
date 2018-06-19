#ifndef GEOMETRY_H
#define GEOMETRY_H

typedef struct geometrySurface
{
	int nodeSize;
	int triangleSize;
	double** Nodes;
	int** Triangles;
	int* triangleSurface;
	int* nodeSurface;
}geometrySurface;

typedef struct geometryVolume
{
	int nodeSize;
	int* nodeSurface;
	double** Nodes;

	int triangleSize;
	int* triangleSurface;
	int** Triangles;

	int tetrahedraSize;
	int* tetrahedraSurface;
	int** Tetrahedras;
}geometryVolume;

geometrySurface openSurface(char* fileName);

void printSurface(char* fileName, double** Nodes, int nodeSize, int* nodeSurface,
	int** Triangles, int triangleSize, int* triangleSurface);

void freeGeometrySurface(geometrySurface surface);

geometryVolume openVolume(char* fileName);

void printVolume(char* fileName, double** Nodes, int nodeSize, int* nodeSurface,
	int** Triangles, int triangleSize, int* triangleSurface,
	int** Tetrahedras, int tetrahedraSize, int* tetrahedraSurface);

void printVolumeValues(char* fileName, double** Nodes, int nodeSize, int* nodeSurface,
	int** Triangles, int triangleSize, int* triangleSurface,
	int** Tetrahedras, int tetrahedraSize, int* tetrahedraSurface, double** Values,
	int valuesSize, char* dataName);

void freeGeometryVolume(geometryVolume volume);
#endif