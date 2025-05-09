#ifndef CONDITIONALP_H
#define CONDITIONALP_H

#include <stdbool.h> 
#include <math.h>    
#include <stdlib.h>  


// void calculate_conditional_matrix_and_bit_probabilities(
//     const double *dist, int bw, int cl, double **conditional_matrix, double *bit_probabilities);

void calculate_all_CPMs(
    const double *dist, int bw, int cl,
    double **CPM_00, double **CPM_01, double **CPM_10, double **CPM_11,
    double *bit_probabilities);

    void print_CPM(const char* name, double** CPM, int bw);

#endif 
