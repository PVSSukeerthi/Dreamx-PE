#ifndef LV_H
#define LV_H

#include <stdbool.h> 
#include <stdlib.h>   
#include <math.h>   

void computeLinkVector(double **Conditional_Matrix, double *Bit_Probabilities, int BW, double alpha, int *LV, int *LV_size);

#endif 
