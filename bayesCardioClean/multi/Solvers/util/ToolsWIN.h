#ifndef TOOLSWIN_H
#define TOOLSWIN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>
//#include <sys/time.h>

void *Malloc( long size );
double **matrixNew( int n, int m );
double vectorDotProduct( double *vector0, double *vector1, int n0 );
double vectorNorm( double *vector0, int n0 );
double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
void blasDSCAL( int n, double sa, double x[], int incx );
int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );
double **choleskyDecomposition( double **matrix, int n );
int linpackDTRDI( double t[], int ldt, int n );
double **matrixLowerTriangularInverse( double **matrix, int n );
void eigenDecomposition( double **matrix, int n, double **D, double **Q );
void eigenDecompositionQLalgo2( int n, double **V, double *d, double *e );
double myhypot( double a, double b );
void eigenDecompositionHouseholder2( int n, double **V, double *d, double *e );
void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 );

int *mergeSort( double *array, int array_size );
void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q );
void mergeSortWithinBoundsInt( int *array, int *sorted, int *tosort, int p, int q );

void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q );
int *mergeSortInt( int *array, int array_size );
void mergeSortMergeInt( int *array, int *sorted, int *tosort, int p, int r, int q );

int *getRanks(double *array, int array_size );
int *getRanksFromSorted(int *sorted, int array_size );

double randomRealUniform01( void );
int randomInt( int maximum );
double random1DNormalUnit( void );
double random1DNormalParameterized( double mean, double variance );
void initializeRandomNumberGenerator( void );
int *randomPermutation( int n );
int **allPermutations( int length, int *numberOfPermutations );
int **allPermutationsSubroutine( int from, int length, int *numberOfPermutations );

long getMilliSecondsRunning();
long getMilliSecondsRunningAfterInit();
long getMilliSecondsRunningSinceTimeStamp( long timestamp );
long getCurrentTimeStampInMilliSeconds();

void startTimer( void );
double getTimer( void );
void printTimer( void );

double maxA( double x, double y );
double minA( double x, double y );
double distanceEuclidean( double *solution_a, double *solution_b, int n );
double distanceEuclidean2D( double x1, double y1, double x2, double y2 );

double *matrixVectorPartialMultiplication( double **matrix, double *vector, int n0, int n1, int number_of_elements, int *element_indices );

int64_t    random_seed,                      /* The seed used for the random-number generator. */
           random_seed_changing;             /* Internally used variable for randomly setting a random seed. */

long  timestamp_start,                       /* The time stamp in milliseconds for when the program was started. */
      timestamp_start_after_init;            /* The time stamp in milliseconds for when the algorithm was started */

double haveNextNextGaussian,             /* Internally used variable for sampling the normal distribution. */
       nextNextGaussian;                     /* Internally used variable for sampling the normal distribution. */

#endif
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Constants -=-=-=-=-=-=-=-=-=-=-=-=-=-*/
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798
#define FALSE 0
#define TRUE 1

