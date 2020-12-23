#include <stdio.h>


int relaxation_matrix_set_ex(hw_params_t* hw_params,
                             double* mat,
                             uint32_t ln,  /* # local elements including ghost in the x dimension */
                             uint32_t lm,  /* # local elements including ghost in the y dimension */
                             int rowpos,   /* my row position on the global matrix */
                             int colpos)   /* my col position on the global matrix */
{
    uint32_t i, j;
    double dist;
    int N = hw_params->resolution;

    for( i = 0; i < hw_params->num_sources; i++ ) {
        /**
         * The heat dissipate linearly in circles around the central
         * point up to the defined range. It only affects the interface
         * between the mediums, so it only has an impact on the boundaries.
         */
        if( 0 == rowpos ) {
            for( j = 1; j < ln-1; j++ ) {  /* initialize the top row */
                dist = sqrt( pow((double)(colpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].x, 2) +
                             pow(hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[rowpos*ln+j] += ((hw_params->heat_sources[i].range - dist) /
                               hw_params->heat_sources[i].range *
                               hw_params->heat_sources[i].temp);
                }
            }
        }
        if( rowpos >= (N - (lm - 2)) ) {
            for( j = 1; j < ln-1; j++ ) {  /* initialize the bottom row */
                dist = sqrt( pow((double)(colpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].x, 2) +
                             pow(1-hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[(lm-1)*ln+j] += ((hw_params->heat_sources[i].range - dist) /
                                         hw_params->heat_sources[i].range *
                                         hw_params->heat_sources[i].temp);
                }
            }
        }
        if( 0 == colpos ) {
            for( j = 1; j < lm-1; j++ ) {  /* left-most column */
                dist = sqrt( pow(hw_params->heat_sources[i].x, 2) +
                             pow((double)(rowpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[j*ln+0] += ((hw_params->heat_sources[i].range - dist) /
                                  hw_params->heat_sources[i].range *
                                  hw_params->heat_sources[i].temp);
                }
            }
        }
        if( colpos >= (N - (ln - 2)) ) {
            for( j = 1; j < lm-1; j++ ) {  /* right-most column */
                dist = sqrt( pow(1-hw_params->heat_sources[i].x, 2) +
                             pow((double)(rowpos+j)/(double)(N+1) -
                                 hw_params->heat_sources[i].y, 2));
                if( dist <= hw_params->heat_sources[i].range ) {
                    mat[j*ln+(ln-1)] += ((hw_params->heat_sources[i].range - dist) /
                                         hw_params->heat_sources[i].range *
                                         hw_params->heat_sources[i].temp);
                }
            }
        }
    }
    return 0;
}