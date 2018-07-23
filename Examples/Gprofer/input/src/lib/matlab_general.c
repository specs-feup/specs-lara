/* Implementation file for lib/matlab_general */

#include "matlab_general.h"


/**
 *  max implementation for Y = max(X,Z), when X and Z are both scalars 
 */
int max_scalars_dec_ii(int scalar1, int scalar2)
{

	/* return the correct element */
	if( scalar1 > scalar2 )

		return scalar1;

	return scalar2;
}
