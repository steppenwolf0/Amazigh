#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openFile.h"
#include "supportC.h"

double** openMatrix(const char* name)
{
	FILE* file = fopen(name, "r"); /* should check the result */
	int n=8192;
    char line[60000];
	int dim_1=atoi(fgets(line, sizeof(line), file));
	int dim_2=atoi(fgets(line, sizeof(line), file));
	double** matrix=newDouble(dim_1,dim_2);
    //while (fgets(line, sizeof(line), file)) 
	for (int i=0;i<dim_1;i++)
	{
		fgets(line, sizeof(line), file);
		char *p;
		p = strtok(line, " ");
       for (int j=0;j<dim_2;j++)
       {
       		if(p)
		    //printf("%s ", p);
		    matrix[i][j]=atof(p);
		    p = strtok(NULL, " ");
	   }
    	//printf("\n");
	}
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
    /* may check feof here to make a difference between eof and io failure -- network
       timeout for instance */
   return matrix;
}

struct matDouble openMatrixDouble(const char* name)
{
	FILE* file = fopen(name, "r"); /* should check the result */
	int n=100;
    char line[100];
	int dim_1=atoi(fgets(line, sizeof(line), file));
	int dim_2=atoi(fgets(line, sizeof(line), file));
	char line2[60000];
	double** matrix=newDouble(dim_1,dim_2);
    //while (fgets(line, sizeof(line), file)) 
	for (int i=0;i<dim_1;i++)
	{
		fgets(line2, sizeof(line2), file);
		char *p;
		p = strtok(line2, " ");
       for (int j=0;j<dim_2;j++)
       {
       		if(p)
		    //printf("%s ", p);
		    matrix[i][j]=atof(p);
		    p = strtok(NULL, " ");
	   }
    	//printf("\n");
	}
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
    /* may check feof here to make a difference between eof and io failure -- network
       timeout for instance */
       
       struct matDouble outMatrix;
       outMatrix.dim_1=dim_1;
       outMatrix.dim_2=dim_2;
       outMatrix.matrix=matrix;
       
       return outMatrix;
}

struct matDouble openMatrixDouble2(const char* name, int dim_1, int dim_2)
{
	FILE* file = fopen(name, "r"); /* should check the result */

	char line2[60000];
	double** matrix=newDouble(dim_1,dim_2);
    //while (fgets(line, sizeof(line), file)) 
	for (int i=0;i<dim_1;i++)
	{
		fgets(line2, sizeof(line2), file);
		char *p;
		p = strtok(line2, " ");
       for (int j=0;j<dim_2;j++)
       {
       		if(p)
		    //printf("%s ", p);
		    matrix[i][j]=atof(p);
		    p = strtok(NULL, " ");
	   }
    	//printf("\n");
	}
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
    /* may check feof here to make a difference between eof and io failure -- network
       timeout for instance */
       
       struct matDouble outMatrix;
       outMatrix.dim_1=dim_1;
       outMatrix.dim_2=dim_2;
       outMatrix.matrix=matrix;
       
       return outMatrix;

	   fclose(file);
}

double* openVector(const char* name, int dim_1)
{
	FILE* file = fopen(name, "r"); /* should check the result */

	char line2[60000];
	double* matrix = newDoubleV(dim_1);
	//while (fgets(line, sizeof(line), file)) 
	for (int i = 0; i<dim_1; i++)
	{
		fgets(line2, sizeof(line2), file);
		
				matrix[i] = atof(line2);
			
		
		//printf("\n");
	}
	/* note that fgets don't strip the terminating \n, checking its
	presence would allow to handle lines longer that sizeof(line) */
	/* may check feof here to make a difference between eof and io failure -- network
	timeout for instance */

	fclose(file);

	return matrix;
}

int** openMatrixInt(const char* name)
{
	FILE* file = fopen(name, "r"); /* should check the result */
	int n = 8192;
	char line[60000];
	int dim_1 = atoi(fgets(line, sizeof(line), file));
	int dim_2 = atoi(fgets(line, sizeof(line), file));
	int** matrix = newInt(dim_1, dim_2);
	//while (fgets(line, sizeof(line), file)) 
	for (int i = 0; i<dim_1; i++)
	{
		fgets(line, sizeof(line), file);
		char *p;
		p = strtok(line, " ");
		for (int j = 0; j<dim_2; j++)
		{
			if (p)
				//printf("%s ", p);
				matrix[i][j] = atoi(p);
			p = strtok(NULL, " ");
		}
		//printf("\n");
	}
	/* note that fgets don't strip the terminating \n, checking its
	presence would allow to handle lines longer that sizeof(line) */
	/* may check feof here to make a difference between eof and io failure -- network
	timeout for instance */
	return matrix;
}