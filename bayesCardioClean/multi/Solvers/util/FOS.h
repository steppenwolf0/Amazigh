
#ifndef FOS_H
#define FOS_H

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "ToolsWIN.h"
#include "SO_optimization.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

typedef struct FOS {
    int length;
    int **sets;
    int *set_length;
} FOS;

/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
void printFOS( FOS *fos );
FOS *readFOSFromFile( FILE *file );
FOS *copyFOS( FOS *f );
FOS *learnLinkageTree( double **covariance_matrix );
int determineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length );
double getSimilarity( int a, int b );
double **computeMIMatrix( double **covariance_matrix, int n );
int *matchFOSElements( FOS *new_FOS, FOS *prev_FOS );
int *hungarianAlgorithm( int** similarity_matrix, int dim );
void hungarianAlgorithmAddToTree(int x, int prevx, short *S, int *prev, int *slack, int *slackx, int* lx, int *ly, int** similarity_matrix, int dim);
int determineNearestNeighbour(int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length );
void ezilaitiniFOS(FOS *lm );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
int      *mpm_number_of_indices,
          FOS_element_ub,                       /* Cut-off value for bounded fixed linkage tree (BFLT). */
          use_univariate_FOS,                   /* Whether a univariate FOS is used. */
          learn_linkage_tree,                   /* Whether the FOS is learned at the start of each generation. */
          static_linkage_tree,                  /* Whether the FOS is fixed throughout optimization. */
          random_linkage_tree,                  /* Whether the fixed linkage tree is learned based on a random distance measure. */
          FOS_element_size;                     /* If positive, the size of blocks of consecutive variables in the FOS. If negative, determines specific kind of linkage tree FOS. */
//double ***MI_matrices;
double  **S_matrix,
         *S_vector;                             /* Avoids quadratic memory requirements when a linkage tree is learned based on a random distance measure. */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

#endif
