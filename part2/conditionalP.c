#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include "conditionalP.h"

#define MAX_BITS 32

void calculate_conditional_matrix_and_bit_probabilities(
    const double *dist, int bw, int cl, double **conditional_matrix, double *bit_probabilities) {

    int num_values = (int)pow(2, bw);
    int *dec_values = (int *)malloc(num_values * sizeof(int));
    
    if (!dec_values) {
        printf("Error: Memory allocation failed for dec_values!\n");
        exit(1);
    }

    for (int i = 0; i < num_values; i++) {
        dec_values[i] = i;
    }

    // **🚀 Allocate Memory for Conditional Matrix Correctly**
    for (int i = 0; i < bw; i++) {
        for (int j = 0; j < bw; j++) {
            conditional_matrix[i][j] = 0.0; // Initialize to avoid garbage values
        }
    }

    // **🚀 Calculating Bit Probabilities**
    for (int i = 0; i < bw; i++) {
        double sum_prob = 0;
        for (int j = 0; j < num_values; j++) {
            if ((dec_values[j] % (int)pow(2, i + 1)) >= (int)pow(2, i)) {
                sum_prob += dist[j];
            }
        }
        bit_probabilities[i] = sum_prob;
    }

    // **🚀 Calculating Conditional Probabilities**
    for (int i = bw - 1; i >= 1; i--) {
        for (int j = i - 1; j >= fmax(0, i - cl); j--) {

            for (int bit_i_val = 0; bit_i_val < 2; bit_i_val++) {
                for (int bit_j_val = 0; bit_j_val < 2; bit_j_val++) {

                    bool *inds_for_i = (bool *)malloc(num_values * sizeof(bool));
                    bool *inds_for_j = (bool *)malloc(num_values * sizeof(bool));

                    if (!inds_for_i || !inds_for_j) {
                        printf("Error: Memory allocation failed for inds_for_i/j!\n");
                        exit(1);
                    }
                    
                    for (int k = 0; k < num_values; k++) {
                        inds_for_i[k] = ((dec_values[k] % (int)pow(2, i + 1)) >= (int)pow(2, i)) ^ bit_i_val;
                        inds_for_j[k] = ((dec_values[k] % (int)pow(2, j + 1)) >= (int)pow(2, j)) ^ bit_j_val;
                    }

                    double num = 0;
                    double denom = 0;
                    for (int k = 0; k < num_values; k++) {
                        if (inds_for_j[k]) {
                            denom += dist[k];
                            if (inds_for_i[k]) {
                                num += dist[k];
                            }
                        }
                    }

                    if (denom > 0) {
                        conditional_matrix[i][j] = num / denom;
                    } else {
                        conditional_matrix[i][j] = 0;
                    }

                    free(inds_for_i);
                    free(inds_for_j);
                }
            }
        }
    }

    free(dec_values);
}

// int main() {
//     int bw = 4; // Example bitwidth
//     int cl = 2; // Example coverage length

//     // Example input distribution (normalized probabilities for 2^bw values)
//     double dist[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

//     double *conditional_matrix[MAX_BITS * MAX_BITS];
//     for (int i = 0; i < MAX_BITS * MAX_BITS; i++) {
//         conditional_matrix[i] = (double *)calloc(4, sizeof(double)); // 4 possible combinations of bit values
//     }

//     double bit_probabilities[MAX_BITS] = {0};

//     calculate_conditional_matrix_and_bit_probabilities(dist, bw, cl, conditional_matrix, bit_probabilities);

//     // Print results
//     printf("Bit Probabilities:\n");
//     for (int i = 0; i < bw; i++) {
//         printf("Bit %d: %.3f\n", i, bit_probabilities[i]);
//     }

//     printf("\nConditional Matrix:\n");
//     for (int i = 0; i < bw; i++) {
//         for (int j = 0; j < bw; j++) {
//             if (i * bw + j < MAX_BITS * MAX_BITS) {
//                 printf("(%d,%d): [%.3f %.3f %.3f %.3f]\n", i, j, conditional_matrix[i * bw + j][0],
//                        conditional_matrix[i * bw + j][1],
//                        conditional_matrix[i * bw + j][2],
//                        conditional_matrix[i * bw + j][3]);
//             }
//         }
//     }

//     for (int i = 0; i < MAX_BITS * MAX_BITS; i++) {
//         free(conditional_matrix[i]);
//     }

//     return 0;
// }
