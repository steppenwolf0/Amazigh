#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include "mathMatrix.h"
#include "supportC.h"
void plot2D(double* x, double* y, int dim, const char* title, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f,"\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");
	fprintf(f, "		var trace1 = { \n");
	fprintf(f, "		x: [");
	for (i = 0; i < dim; i++)
	{
		fprintf(f,"%.4f, ", x[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		y: [");
	for (i = 0; i < dim; i++)
	{
		fprintf(f, "%.4f, ", y[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		mode : 'lines+markers' \n");
	fprintf(f, "		}; \n");

	fprintf(f,	"		var data = [ trace1 ]; \n");
	
	fprintf(f, "		var layout = { \n");
	fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");
		
	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");
		


	fclose(f);
}

void plot2DMat(double* x, double** y,  int dim_1, int dim_2, const char* title, const char* name )
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i, j;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");

	for (j = 0; j < dim_1; j++)
	{
		fprintf(f, "<div id = \"myDiv");
		fprintf(f, "%d", j);
		fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
		fprintf(f, "	<script>\n");
		fprintf(f, "		var trace1 = { \n");
		fprintf(f, "		x: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", x[i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		y: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", y[j][i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		mode : 'lines+markers' \n");
		fprintf(f, "		}; \n");

		fprintf(f, "		var data = [ trace1 ]; \n");

		fprintf(f, "		var layout = { \n");
		fprintf(f, "		title:'");
		fprintf(f, title);
		fprintf(f, " %d",j);
		fprintf(f, "' \n");

		fprintf(f, "		}; \n");

		fprintf(f, "		Plotly.newPlot('myDiv");
		fprintf(f, "%d", j);
		fprintf(f, "', data, layout); \n");
		fprintf(f, "	</script> \n");

		fprintf(f, "</body> \n");
	}

	
	fclose(f);
}

void plot2DComparisonNames(double* x, double** y, int dim_1, int dim_2, const char* title, const char* name, char** variables)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i, j;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");

	for (j = 0; j < dim_1; j++)
	{
		fprintf(f, "		var trace");
		fprintf(f, "%d = { \n", j);
		fprintf(f, "		x: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", x[i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		y: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", y[j][i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		name : '%s', \n",variables[j]);
		fprintf(f, "		mode : 'lines' \n");
		fprintf(f, "		}; \n");
	}




	fprintf(f, "		var data = [ ");
	for (j = 0; j < dim_1; j++)
	{
		fprintf(f, "trace%d,", j);
	}
	fprintf(f, " ]; \n");

	fprintf(f, "		var layout = { \n");
	/*fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");*/
	//X AXIS
	//fprintf(f, "		,xaxis: { \n");
	fprintf(f, "		xaxis: { \n");
	fprintf(f, "		title: 'iterations'");
	fprintf(f, "		},");
	//Y AXIS
	fprintf(f, "		yaxis: { \n");
	fprintf(f, "		title: 'log10(x)'");
	fprintf(f, "		}");

	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");



	fclose(f);
}


void plot2DComparison(double* x, double** y, int dim_1, int dim_2, const char* title, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i, j;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");

	for (j = 0; j < dim_1; j++)
	{
		fprintf(f, "		var trace");
		fprintf(f,"%d = { \n",j);
		fprintf(f, "		x: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", x[i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		y: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", y[j][i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		mode : 'lines' \n");
		fprintf(f, "		}; \n");
	}
	

	

	fprintf(f, "		var data = [ ");
	for (j = 0; j < dim_1; j++)
	{
		fprintf(f, "trace%d,", j);
	}
	fprintf(f," ]; \n");

	fprintf(f, "		var layout = { \n");
	fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");
	//X AXIS
	fprintf(f, "		,xaxis: { \n");
	fprintf(f, "		title: 'iterations'");
	fprintf(f, "		},");
	//Y AXIS
	fprintf(f, "		yaxis: { \n");
	fprintf(f, "		title: 'log10(x)'");
	fprintf(f, "		}");

	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");



	fclose(f);
}

void plot2DExtraDots(double* x, double* y, double** dots,  int dim_2,int dim_3, const char* title, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i, j;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");


	 j = 0;
		fprintf(f, "		var trace");
		fprintf(f, "%d = { \n", j);
		fprintf(f, "		x: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", x[i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		y: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", y[i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		mode : 'lines+markers' \n");
		fprintf(f, "		}; \n");
	
		j = 1;
		fprintf(f, "		var trace");
		fprintf(f, "%d = { \n", j);
		fprintf(f, "		x: [");
		for (i = 0; i < dim_3; i++)
		{
			fprintf(f, "%.4f, ", dots[0][i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		y: [");
		for (i = 0; i < dim_3; i++)
		{
			fprintf(f, "%.4f, ", dots[1][i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		mode : 'lines+markers' \n");
		fprintf(f, "		}; \n");



	fprintf(f, "		var data = [ ");
	for (j = 0; j < 2; j++)
	{
		fprintf(f, "trace%d,", j);
	}
	fprintf(f, " ]; \n");

	fprintf(f, "		var layout = { \n");
	fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");

	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");



	fclose(f);
}

void plotContour(double* x, double* y, double** z, int dim_1,int dim_2, const char* title, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i, j;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");
	fprintf(f, "		var trace1 = { \n");
	fprintf(f, "		x: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", x[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		y: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", y[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		z: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "[");
		for (j = 0; j < dim_2; j++)
		{
			fprintf(f, "%.4f, ", z[i][j]);
		}
		fprintf(f, "],\n");
	}
	fprintf(f, "], \n");
	fprintf(f, "		type: 'contour' \n");
	fprintf(f, "		}; \n");

	fprintf(f, "		var data = [ trace1 ]; \n");

	fprintf(f, "		var layout = { \n");
	fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");

	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");



	fclose(f);
}

void plotPoints(double* x, double* y, double* z, int dim_1, const char* title, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");
	fprintf(f, "		var trace1 = { \n");
	fprintf(f, "		x: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", x[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		y: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", y[i]);
	}
	fprintf(f, "], \n");

	

	fprintf(f, "		mode: 'markers', \n");

	fprintf(f, "			marker: { \n");
	fprintf(f, "				size: 10, \n");
	fprintf(f, "				color: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", z[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "				} \n");


	fprintf(f, "		}; \n");

	fprintf(f, "		var data = [ trace1 ]; \n");

	fprintf(f, "		var layout = { \n");
	fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");

	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");



	fclose(f);
}

void plotPointsExtraDots(double* x, double* y, double* z, double** dots, int dim_1, int dim_2, const char* title, const char* name)
{
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");
	fprintf(f, "		var trace1 = { \n");
	fprintf(f, "		x: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", x[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		y: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", y[i]);
	}
	fprintf(f, "], \n");



	fprintf(f, "		mode: 'markers', \n");

	fprintf(f, "			marker: { \n");
	fprintf(f, "				size: 10, \n");
	fprintf(f, "				color: [");
	for (i = 0; i < dim_1; i++)
	{
		fprintf(f, "%.4f, ", z[i]);
	}
	fprintf(f, "], \n");
	fprintf(f, "				} \n");


	fprintf(f, "		}; \n");


	fprintf(f, "		var trace");
	fprintf(f, "%d = { \n", 2);
	fprintf(f, "		x: [");
	for (i = 0; i < dim_2; i++)
	{
		fprintf(f, "%.4f, ", dots[i][0]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		y: [");
	for (i = 0; i < dim_2; i++)
	{
		fprintf(f, "%.4f, ", dots[i][1]);
	}
	fprintf(f, "], \n");
	fprintf(f, "		mode : 'markers' \n");
	fprintf(f, "		}; \n");

	fprintf(f, "		var data = [ trace1,trace2 ]; \n");

	fprintf(f, "		var layout = { \n");
	fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");

	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");



	fclose(f);
}

//*********************************************************************************NEW

void plot2DComparisonNamesNew(double** y, int dim_1, int dim_2, const char* title, const char* name, char** variables, const char* Yaxis)
{
	double** yT = transposeMatrix(y, dim_2, dim_1);
	double* x = newDoubleV(dim_2);

	for (int i = 0; i < dim_2; i++)
	{
		 x[i]=i;

	}
	FILE *f = fopen(name, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	int i, j;
	fprintf(f, "<head> \n");
	fprintf(f, "<script src = \"plotly-latest.min.js\"> \n");
	fprintf(f, "</script> \n");
	fprintf(f, "</head> \n");
	fprintf(f, "<body> \n");
	fprintf(f, "<div id = \"myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "\"><!--Plotly chart will be drawn inside this DIV--></div> \n");
	fprintf(f, "	<script>\n");

	for (j = 0; j < dim_1; j++)
	{
		fprintf(f, "		var trace");
		fprintf(f, "%d = { \n", j);
		fprintf(f, "		x: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", x[i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		y: [");
		for (i = 0; i < dim_2; i++)
		{
			fprintf(f, "%.4f, ", yT[j][i]);
		}
		fprintf(f, "], \n");
		fprintf(f, "		name : '%s', \n", variables[j]);
		fprintf(f, "		mode : 'lines' \n");
		fprintf(f, "		}; \n");
	}




	fprintf(f, "		var data = [ ");
	for (j = 0; j < dim_1; j++)
	{
		fprintf(f, "trace%d,", j);
	}
	fprintf(f, " ]; \n");

	fprintf(f, "		var layout = { \n");
	fprintf(f, "		title:'");
	fprintf(f, title);
	fprintf(f, "' \n");
	//X AXIS
	fprintf(f, "		,xaxis: { \n");
	//fprintf(f, "		xaxis: { \n");
	fprintf(f, "		title: 'iterations'");
	fprintf(f, "		},");
	//Y AXIS
	fprintf(f, "		yaxis: { \n");
	fprintf(f, "		title: '");
	fprintf(f, Yaxis);
	fprintf(f, "'");
	fprintf(f, "		}");

	fprintf(f, "		}; \n");

	fprintf(f, "		Plotly.newPlot('myDiv");
	fprintf(f, "%d", 0);
	fprintf(f, "', data, layout); \n");
	fprintf(f, "	</script> \n");



	fprintf(f, "</body> \n");



	fclose(f);

	free(x);
	freeMatrix(yT, dim_1);
}

