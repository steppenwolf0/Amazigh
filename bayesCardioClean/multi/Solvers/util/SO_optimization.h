#ifndef SO_OPTIMIZATION_H
#define SO_OPTIMIZATION_H

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "FOS.h"


#ifdef _WIN32
#include "ToolsWIN.h"
#else
#include "Tools.h"
#endif

#include "../problemDef.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/



/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
double installedProblemLowerRangeBound(  int dimension );
double installedProblemUpperRangeBound(  int dimension );
void initializeParameterRangeBounds( void );
short isParameterInRangeBounds( double parameter, int dimension );
void installedProblemEvaluation( Individual* ind, int number_of_touched_parameters, int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before );
void installedProblemEvaluationWithoutRotation(Individual* ind,
	int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before );


/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
double elitist_objective_value,
       elitist_constraint_value;
//Individual **populations;
//Individual **selections;
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
short  black_box_evaluations,                         /* Whether full (black-box) evaluations must always be performed. */
       use_vtr,                                       /* Whether to terminate at the value-to-reach (VTR) (0 = no). */
       vtr_hit_status,                                /* Whether the VTR has been reached. */
      *populations_terminated,                        /* Which populations have been terminated. */
       evaluations_for_statistics_hit,                /* Can be used to write statistics after a certain number of evaluations. */
       write_generational_statistics,                 /* Whether to compute and write statistics every generation (0 = no). */
       write_generational_solutions;                  /* Whether to write the population every generation (0 = no). */
int 
       //number_of_parameters,                          /* The number of parameters to be optimized. */
       number_of_populations,                         /* The number of parallel populations that initially partition the search space. */
       block_size,                                    /* The number of variables in one block of the 'sum of rotated ellipsoid blocks' function. */
       number_of_blocks,                              /* The number of blocks the 'sum of rotated ellipsoid blocks' function. */
       block_start,                                   /* The index at which the first block starts of the 'sum of rotated ellipsoid blocks' function. */
      *number_of_generations,                         /* The current generation count of a subgeneration in the interleaved multi-start scheme. */
       total_number_of_generations,                   /* The overarching generation count in the interleaved multi-start scheme. */
      *population_sizes;                              /* The size of the population. */
double number_of_evaluations,                         /* The current number of times a function evaluation was performed. */
       vtr,                                           /* The value-to-reach (function value of best solution that is feasible). */
      *lower_range_bounds,                            /* The respected lower bounds on parameters. */
      *upper_range_bounds,                            /* The respected upper bounds on parameters. */
      *lower_init_ranges,                             /* The initialization range lower bound. */
      *upper_init_ranges,                             /* The initialization range upper bound */
       lower_user_range,                              /* The initial lower range-bound indicated by the user (same for all dimensions). */
       upper_user_range;                              /* The initial upper range-bound indicated by the user (same for all dimensions). */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

#endif
