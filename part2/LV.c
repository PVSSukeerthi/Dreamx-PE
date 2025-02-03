#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define MAX_BITS 32

void computeLinkVector(double** Conditional_Matrix, double* Bit_Probabilities, int BW, double alpha, int* LV, int* LV_size) {
    int i, j, idx = 0;
    double** Dependence_Matrix = (double**)malloc(BW * sizeof(double*));
    bool** Dependence = (bool**)malloc(BW * sizeof(bool*));

    for (i = 0; i < BW; i++) {
        Dependence_Matrix[i] = (double*)calloc(BW, sizeof(double));
        Dependence[i] = (bool*)calloc(BW, sizeof(bool));
    }

    // Dependence Matrix
    for (i = 1; i < BW; i++) { // Start from the second bit
        for (j = 0; j < i; j++) { // Compare with previous bits
            Dependence_Matrix[i][j] = fabs(Conditional_Matrix[i][j] - Bit_Probabilities[i]);
            if (Dependence_Matrix[i][j] >= alpha) {
                Dependence[i][j] = true;
            }
        }
    }

    int ranges_min[MAX_BITS] = {0}, ranges_max[MAX_BITS] = {0};
    int ranges_count = 0;

    // Identify dependent ranges
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

    idx = 0;
    int current_bit = 0;

    // Merge dependent and independent bits in order
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

    // Free allocated memory
    for (i = 0; i < BW; i++) {
        free(Dependence_Matrix[i]);
        free(Dependence[i]);
    }
    free(Dependence_Matrix);
    free(Dependence);
}

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
