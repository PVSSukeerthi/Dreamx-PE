#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "conditionalP.h"

#define MAX_BITS 32  // Adjust based on bit-width

void computeLinkVector(
    double** CPM_00, double** CPM_01, double** CPM_10, double** CPM_11, 
    double* Bit_Probabilities, int BW, double alpha, int* LV, int* LV_size) {

    int i, j, idx = 0;
    
    // 🚀 Allocate Dependence Matrix (BW x BW)
    bool** Dependence = (bool**)malloc(BW * sizeof(bool*));
    for (i = 0; i < BW; i++) {
        Dependence[i] = (bool*)calloc(BW, sizeof(bool));
    }

    // 🚀 Step 1: Compute Dependence Matrix
    for (i = 1; i < BW; i++) { // Start from bit 1 (MSB)
        for (j = 0; j < i; j++) { // Compare with previous bits
            
            double delta_00 = fabs(CPM_00[i][j] - Bit_Probabilities[i]);
            double delta_01 = fabs(CPM_01[i][j] - Bit_Probabilities[i]);
            double delta_10 = fabs(CPM_10[i][j] - Bit_Probabilities[i]);
            double delta_11 = fabs(CPM_11[i][j] - Bit_Probabilities[i]);

            // If ANY of the 4 CPM values indicate strong dependence, mark as dependent
            if (delta_00 >= alpha || delta_01 >= alpha || delta_10 >= alpha || delta_11 >= alpha) {
                Dependence[i][j] = true;
            }
        }
    }

    // 🚀 Step 2: Identify Dependent Ranges
    int ranges_min[MAX_BITS] = {0}, ranges_max[MAX_BITS] = {0};
    int ranges_count = 0;

    for (i = 0; i < BW; i++) {
        bool hasDependence = false;
        for (j = 0; j < BW; j++) {
            if (Dependence[i][j]) {
                hasDependence = true;
                break;
            }
        }
        if (hasDependence) {
            if (ranges_count == 0 || i > ranges_max[ranges_count - 1] + 1) {
                ranges_min[ranges_count] = i;
                ranges_max[ranges_count] = i;
                ranges_count++;
            } else {
                ranges_max[ranges_count - 1] = i;
            }
        }
    }

    // 🚀 Step 3: Generate Link Vector (LV)
    idx = 0;
    int current_bit = 0;

    for (i = 0; i < ranges_count; i++) {
        // Add independent bits before the current dependent range
        while (current_bit < ranges_min[i]) {
            LV[idx++] = current_bit++;
            LV[idx++] = -1; // Separator for independent bit
        }

        // Add the current dependent range
        for (j = ranges_min[i]; j <= ranges_max[i]; j++) {
            LV[idx++] = j;
        }
        LV[idx++] = -1; // Separator for dependent range

        current_bit = ranges_max[i] + 1;
    }

    // Add remaining independent bits after the last dependent range
    while (current_bit < BW) {
        LV[idx++] = current_bit++;
        LV[idx++] = -1; // Separator for independent bit
    }

    *LV_size = idx;

    // 🚀 Step 4: Free allocated memory
    for (i = 0; i < BW; i++) {
        free(Dependence[i]);
    }
    free(Dependence);
}

// int main() {
//     int BW = 8;  // Example: 8-bit adder
//     int CL = 3;  // Coverage length (how many previous bits to check)
//     double alpha = 0.05; // Dependence threshold

//     // Example: Probability distribution (Replace with real data)
//     double dist[256]; // For 8-bit, 2^8 = 256
//     for (int i = 0; i < 256; i++) dist[i] = 1.0 / 256;  // Uniform distribution

//     // Allocate memory for 4 CPMs
//     double **CPM_00, **CPM_01, **CPM_10, **CPM_11;
//     CPM_00 = (double **)malloc(BW * sizeof(double *));
//     CPM_01 = (double **)malloc(BW * sizeof(double *));
//     CPM_10 = (double **)malloc(BW * sizeof(double *));
//     CPM_11 = (double **)malloc(BW * sizeof(double *));
//     for (int i = 0; i < BW; i++) {
//         CPM_00[i] = (double *)calloc(BW, sizeof(double));
//         CPM_01[i] = (double *)calloc(BW, sizeof(double));
//         CPM_10[i] = (double *)calloc(BW, sizeof(double));
//         CPM_11[i] = (double *)calloc(BW, sizeof(double));
//     }

//     // Allocate memory for bit probabilities
//     double bit_probabilities[BW];

//     // 🚀 Compute the CPMs using the actual probability distribution
//     calculate_all_CPMs(dist, BW, CL, CPM_00, CPM_01, CPM_10, CPM_11, bit_probabilities);

//     // 🚀 Print all CPM matrices for verification
//     print_CPM("CPM{0,0}", CPM_00, BW);
//     print_CPM("CPM{0,1}", CPM_01, BW);
//     print_CPM("CPM{1,0}", CPM_10, BW);
//     print_CPM("CPM{1,1}", CPM_11, BW);

//     // Allocate memory for Link Vector
//     int LV[MAX_BITS], LV_size = 0;

//     // 🚀 Compute Link Vector using the computed CPMs
//     computeLinkVector(CPM_00, CPM_01, CPM_10, CPM_11, bit_probabilities, BW, alpha, LV, &LV_size);

//     // 🚀 Print the Link Vector
//     printf("\nLink Vector (LV):\n");
//     for (int i = 0; i < LV_size; i++) {
//         if (LV[i] == -1) printf("-1 "); // Separator
//         else printf("%d ", LV[i]);
//     }
//     printf("\n");

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
//     // Example bitwidth
//     int BW = 8;  // Number of bits (you can change this as needed)
//     double alpha = 0.5;  // Threshold for dependence

//     // Allocate and initialize Conditional_Matrix
//     double** Conditional_Matrix = (double**)malloc(BW * sizeof(double*));
//     for (int i = 0; i < BW; i++) {
//         Conditional_Matrix[i] = (double*)malloc(BW * sizeof(double));
//     }

//     // Example initialization of Conditional_Matrix
//     for (int i = 0; i < BW; i++) {
//         for (int j = 0; j < BW; j++) {
//             Conditional_Matrix[i][j] = (double)(rand() % 100) / 100.0;  // Random values between 0 and 1
//         }
//     }

//     // Allocate and initialize Bit_Probabilities
//     double* Bit_Probabilities = (double*)malloc(BW * sizeof(double));

//     // Example initialization of Bit_Probabilities
//     for (int i = 0; i < BW; i++) {
//         Bit_Probabilities[i] = (double)(rand() % 100) / 100.0;  // Random values between 0 and 1
//     }

//     // Allocate memory for the Link Vector (LV)
//     int* LV = (int*)malloc(MAX_BITS * sizeof(int));  // Maximum possible size
//     int LV_size = 0;

//     // Call the computeLinkVector function
//     computeLinkVector(Conditional_Matrix, Bit_Probabilities, BW, alpha, LV, &LV_size);

//     // Print the resulting Link Vector
//     printf("Link Vector: ");
//     for (int i = 0; i < LV_size; i++) {
//         printf("%d ", LV[i]);
//     }
//     printf("\n");

//     // Free allocated memory
//     for (int i = 0; i < BW; i++) {
//         free(Conditional_Matrix[i]);
//     }
//     free(Conditional_Matrix);
//     free(Bit_Probabilities);
//     free(LV);

//     return 0;
// }
