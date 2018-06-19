#ifndef FEM_H
#define FEM_H
double** Stiffness3DforFEM(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface);

double** Mass3DforFEM(double** Nodes, int** Triangles, int** Tetrahedra,
	int** NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface);

int neumanTriangles(int* triangleSurface, int triangleSize, int neumannSurface);

void getNeumanData3DNEW(double** Nodes, int** Triangles, int nodeSize, int triangleSize,
	int* triangleSurface, int neumannSurface, double** NeumanData, int** Neuman, int ntriangles, 
	double t);

double* getDirichletData(double** Nodes, int nodeSize, int constrainedSurface, int* nodeSurface
	, double t);

double** Load3D(double** Nodes, int** Tetrahedra, int* NodeSurface, int nodeSize,
	int tetrahedraSize, int constrainedsurface, double* nodeDirichlet, double** NeumanData,
	int** Neuman, int ntriangles, double t);

int getPositions(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int*** Positions);

double* getStiffnessValues(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int** Positions, int sizePositions);

double* getMassValues(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int** Positions, int sizePositions);

double* getStiffnessValuesTensor(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int** Positions, int sizePositions,
	double alongfibers_cond, double perpenfibers_cond, double normalfibers_cond, double** fibers);

double* getStiffnessValuesTensorLesion(double** Nodes, int** Triangles, int** Tetrahedra,
	int* NodeSurface, int nodeSize, int triangleSize, int tetrahedraSize,
	int constrainedsurface, int** Positions, int sizePositions,
	double alongfibers_cond, double perpenfibers_cond, double normalfibers_cond, double** fibers,
	int lesionPosition, double reductionCoeff);
#endif