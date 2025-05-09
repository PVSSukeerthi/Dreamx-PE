#define MAX_BITS 32

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// Function to allocate memory for CPMs
double** allocate_CPM(int bw) {
    double** CPM = (double **)malloc(bw * sizeof(double *));
    if (!CPM) {
        printf("Error: Memory allocation failed for CPM!\n");
        exit(1);
    }
    for (int i = 0; i < bw; i++) {
        CPM[i] = (double *)calloc(bw, sizeof(double)); // Use calloc to initialize to zero
        if (!CPM[i]) {
            printf("Error: Memory allocation failed for CPM[%d]!\n", i);
            exit(1);
        }
    }
    return CPM;
}

// Function to calculate all four Conditional Probability Matrices (CPMs)
void calculate_all_CPMs(
    const double *dist, int bw, int cl,
    double **CPM_00, double **CPM_01, double **CPM_10, double **CPM_11,
    double *bit_probabilities) {

    int num_values = (int)pow(2, bw);
    int *dec_values = (int *)malloc(num_values * sizeof(int));
    
    if (!dec_values) {
        printf("Error: Memory allocation failed for dec_values!\n");
        exit(1);
    }

    for (int i = 0; i < num_values; i++) {
        dec_values[i] = i;
    }

    // 🚀 Initialize all four CPM matrices
    for (int i = 0; i < bw; i++) {
        for (int j = 0; j < bw; j++) {
            CPM_00[i][j] = 0.0;
            CPM_01[i][j] = 0.0;
            CPM_10[i][j] = 0.0;
            CPM_11[i][j] = 0.0;
        }
    }

    // 🚀 Calculate bit probabilities
    for (int i = 0; i < bw; i++) {
        double sum_prob = 0;
        for (int j = 0; j < num_values; j++) {
            if ((dec_values[j] % (int)pow(2, i + 1)) >= (int)pow(2, i)) {
                sum_prob += dist[j];
            }
        }
        bit_probabilities[i] = sum_prob;
    }

    // 🚀 Compute all four CPMs
    for (int i = bw - 1; i >= 1; i--) {
        for (int j = i - 1; j >= fmax(0, i - cl); j--) {
            double num_00 = 0, denom_00 = 0;
            double num_01 = 0, denom_01 = 0;
            double num_10 = 0, denom_10 = 0;
            double num_11 = 0, denom_11 = 0;

            for (int k = 0; k < num_values; k++) {
                int bit_i = (dec_values[k] >> i) & 1;  // Extract bit i
                int bit_j = (dec_values[k] >> j) & 1;  // Extract bit j

                if (bit_j == 0) {  // Conditioning on aj = 0
                    denom_00 += dist[k];  
                    denom_10 += dist[k];  
                    if (bit_i == 0) num_00 += dist[k];  
                    if (bit_i == 1) num_10 += dist[k];  
                } 
                else {  // Conditioning on aj = 1
                    denom_01 += dist[k];  
                    denom_11 += dist[k];  
                    if (bit_i == 0) num_01 += dist[k];  
                    if (bit_i == 1) num_11 += dist[k];  
                }
            }

            // 🚀 Store computed probabilities in CPMs
            CPM_00[i][j] = (denom_00 > 0) ? (num_00 / denom_00) : 0;
            CPM_01[i][j] = (denom_01 > 0) ? (num_01 / denom_01) : 0;
            CPM_10[i][j] = (denom_10 > 0) ? (num_10 / denom_10) : 0;
            CPM_11[i][j] = (denom_11 > 0) ? (num_11 / denom_11) : 0;
        }
    }

    free(dec_values);
}

// Function to print CPM matrices
void print_CPM(const char* name, double** CPM, int bw) {
    printf("\n%s:\n", name);
    printf("    ");
    for (int j = 0; j < bw; j++) printf("%4d ", j);
    printf("\n");

    for (int i = 0; i < bw; i++) {
        printf("%2d: ", i);
        for (int j = 0; j < bw; j++) {
            printf("%4.2f ", CPM[i][j]);
        }
        printf("\n");
    }
}

// int main() {
//     int BW = 8;  
//     int CL = 3;  
//     double alpha = 0.05;  

//     // Example: Probability distribution (Replace with real data)
//     double dist[256]; 
//     for (int i = 0; i < 256; i++) dist[i] = 1.0 / 256;  

//     // 🚀 Allocate CPM matrices
//     double **CPM_00 = allocate_CPM(BW);
//     double **CPM_01 = allocate_CPM(BW);
//     double **CPM_10 = allocate_CPM(BW);
//     double **CPM_11 = allocate_CPM(BW);

//     // Allocate memory for bit probabilities
//     double bit_probabilities[BW];

//     // 🚀 Compute the CPMs
//     calculate_all_CPMs(dist, BW, CL, CPM_00, CPM_01, CPM_10, CPM_11, bit_probabilities);

//     // 🚀 Print all CPM matrices
//     print_CPM("CPM{0,0}", CPM_00, BW);
//     print_CPM("CPM{0,1}", CPM_01, BW);
//     print_CPM("CPM{1,0}", CPM_10, BW);
//     print_CPM("CPM{1,1}", CPM_11, BW);

//     // 🚀 Free allocated memory
//     for (int i = 0; i < BW; i++) {
//         free(CPM_00[i]);
//         free(CPM_01[i]);
//         free(CPM_10[i]);
//         free(CPM_11[i]);
//     }
//     free(CPM_00);
//     free(CPM_01);
//     free(CPM_10);
//     free(CPM_11);

//     return 0;
// }

// int main() {
//     int bw = 8;  // Example: 8-bit adder
//     int cl = 3;  // Coverage length (how many previous bits to check)
//     double dist[256]; // Example probability distribution (replace with real data)
    
//     // Fill `dist` with some example data
//     for (int i = 0; i < 256; i++) dist[i] = 1.0 / 256;  // Uniform distribution

//     // Allocate memory for the four CPM matrices
//     double **CPM_00 = (double **)malloc(bw * sizeof(double *));
//     double **CPM_01 = (double **)malloc(bw * sizeof(double *));
//     double **CPM_10 = (double **)malloc(bw * sizeof(double *));
//     double **CPM_11 = (double **)malloc(bw * sizeof(double *));
//     for (int i = 0; i < bw; i++) {
//         CPM_00[i] = (double *)malloc(bw * sizeof(double));
//         CPM_01[i] = (double *)malloc(bw * sizeof(double));
//         CPM_10[i] = (double *)malloc(bw * sizeof(double));
//         CPM_11[i] = (double *)malloc(bw * sizeof(double));
//     }

//     // Allocate memory for bit probabilities
//     double bit_probabilities[bw];

//     // Call function to compute all four CPMs
//     calculate_all_CPMs(dist, bw, cl, CPM_00, CPM_01, CPM_10, CPM_11, bit_probabilities);

//     // Print all CPM matrices
//     print_CPM("CPM{0,0}", CPM_00, bw);
//     print_CPM("CPM{0,1}", CPM_01, bw);
//     print_CPM("CPM{1,0}", CPM_10, bw);
//     print_CPM("CPM{1,1}", CPM_11, bw);

//     // Free allocated memory
//     for (int i = 0; i < bw; i++) {
//         free(CPM_00[i]);
//         free(CPM_01[i]);
//         free(CPM_10[i]);
//         free(CPM_11[i]);
//     }
//     free(CPM_00);
//     free(CPM_01);
//     free(CPM_10);
//     free(CPM_11);

//     return 0;
// }


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
