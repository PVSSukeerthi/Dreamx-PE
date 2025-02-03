#ifndef CONDITIONALP_H
#define CONDITIONALP_H

#include <stdbool.h> 
#include <math.h>    
#include <stdlib.h>  


void calculate_conditional_matrix_and_bit_probabilities(
    const double *dist, int bw, int cl, double **conditional_matrix, double *bit_probabilities);

#endif 
