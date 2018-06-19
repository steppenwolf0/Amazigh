/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include "SO_optimization.h"
#include "../problemDef.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/**
 * Returns 1 if x is better than y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x has a smaller objective value than y
 */
short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
    short result;

    result = 0;

    if( constraint_value_x > 0 ) /* x is infeasible */
    {
        if( constraint_value_y > 0 ) /* Both are infeasible */
        {
            if( constraint_value_x < constraint_value_y )
                result = 1;
        }
    }
    else /* x is feasible */
    {
        if( constraint_value_y > 0 ) /* x is feasible and y is not */
            result = 1;
        else /* Both are feasible */
        {
            if( objective_value_x < objective_value_y )
                result = 1;
        }
    }

    return( result );
}



/**
 * Returns the lower-range bound of an installed problem.
 */
double installedProblemLowerRangeBound( int dimension )
{
   return( functionProblemLowerRangeBound( dimension ) );

}

/**
 * Returns the upper-range bound of an installed problem.
 */
double installedProblemUpperRangeBound( int dimension )
{
     return( functionProblemUpperRangeBound( dimension ) );

}

/**
 * Returns whether a parameter is inside the range bound of
 * every problem.
 */
short isParameterInRangeBounds( double parameter, int dimension )
{
    if( parameter < installedProblemLowerRangeBound(  dimension ) ||
            parameter > installedProblemUpperRangeBound( dimension ) ||
            isnan( parameter ) )
    {
        return( 0 );
    }

    return( 1 );
}

/**
 * Initializes the parameter range bounds.
 */
void initializeParameterRangeBounds( void )
{
    int i;

    lower_range_bounds = (double *) Malloc( number_of_parameters*sizeof( double ) );
    upper_range_bounds = (double *) Malloc( number_of_parameters*sizeof( double ) );
    lower_init_ranges  = (double *) Malloc( number_of_parameters*sizeof( double ) );
    upper_init_ranges  = (double *) Malloc( number_of_parameters*sizeof( double ) );

    for( i = 0; i < number_of_parameters; i++ )
    {
        lower_range_bounds[i] = installedProblemLowerRangeBound(  i );
        upper_range_bounds[i] = installedProblemUpperRangeBound(  i );
    }

    for( i = 0; i < number_of_parameters; i++ )
    {
        lower_init_ranges[i] = lower_user_range;
        if( lower_user_range < lower_range_bounds[i] )
            lower_init_ranges[i] = lower_range_bounds[i];
        if( lower_user_range > upper_range_bounds[i] )
            lower_init_ranges[i] = lower_range_bounds[i];

        upper_init_ranges[i] = upper_user_range;
        if( upper_user_range > upper_range_bounds[i] )
            upper_init_ranges[i] = upper_range_bounds[i];
        if( upper_user_range < lower_range_bounds[i] )
            upper_init_ranges[i] = upper_range_bounds[i];
    }
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * function after rotating the parameter vector.
 * Both are returned using pointer variables.
 * Number of evaluations is increased by the
 * ratio ([0..1]) of new parameters that have
 * been changed.
 */
void installedProblemEvaluation(  Individual* ind, 	int number_of_touched_parameters, 
	int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before )
{
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //double *touched_parameters, *touched_parameters_before, last_obj, last_cons;
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	double *touched_parameters, last_obj, last_cons;
    int i, j, c, prev_block, cur_block, *block_indices;


    touched_parameters = NULL;
    if( !(touched_parameters_indices == NULL || black_box_evaluations) )
    {
        touched_parameters = (double*) Malloc( number_of_touched_parameters*sizeof( double ) );
        for( i = 0; i < number_of_touched_parameters; i++ )
            touched_parameters[i] = ind->parameters[touched_parameters_indices[i]];
		
    }


   
        if( touched_parameters_indices == NULL || black_box_evaluations ) 
		installedProblemEvaluationWithoutRotation(ind, number_of_parameters, NULL, NULL, NULL, 0, 0 );
        else 
		installedProblemEvaluationWithoutRotation( ind, number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, 
		objective_value_before, constraint_value_before );
    


    if( use_vtr && !vtr_hit_status && ind->constrained_value == 0 && ind->objective_value <= vtr  )
    {
        if( touched_parameters_indices != NULL )
            installedProblemEvaluation(ind, number_of_parameters, NULL, NULL, 0, 0 );
        if(ind->constrained_value == 0 && ind->objective_value <= vtr  )
        {
            vtr_hit_status = 1;
            elitist_objective_value = ind->objective_value;
            elitist_constraint_value = ind->constrained_value;
        }
    }

    if( !vtr_hit_status && betterFitness(ind->objective_value, ind->constrained_value, elitist_objective_value, elitist_constraint_value) )
    {
        elitist_objective_value = ind->objective_value;
        elitist_constraint_value = ind->constrained_value;
    }

    if( touched_parameters_indices != NULL )
        free( touched_parameters );

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * without rotating the parameter vector.
 * Both are returned using pointer variables.
 */
void installedProblemEvaluationWithoutRotation( Individual* ind,
	int number_of_touched_parameters, int *touched_parameters_indices, double *touched_parameters, double *parameters_before, double objective_value_before, double constraint_value_before )
{
    ind->objective_value  = 0.0;
	ind->constrained_value = 0.0;



    if( black_box_evaluations || touched_parameters_indices == NULL )
    {
		
         FunctionProblemEvaluation(ind, number_of_parameters );

        
        number_of_evaluations++;
    }
    else
    {
         FunctionPartialProblemEvaluation(ind,number_of_touched_parameters, touched_parameters_indices, touched_parameters, parameters_before, 
		 objective_value_before, constraint_value_before ); 

        
        number_of_evaluations += number_of_touched_parameters/(double)number_of_parameters;
    }
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
