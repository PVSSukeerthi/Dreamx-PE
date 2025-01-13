#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define MAX_BITS 32 // Maximum bitwidth

// Function to compute the Link Vector (LV)
void computeLinkVector(double** Conditional_Matrix, double* Bit_Probabilities, int BW, double alpha, int* LV, int* LV_size) {
    int i, j, k, l, idx = 0;
    double** Dependence_Matrix = (double**)malloc(BW * sizeof(double*));
    bool** Dependence = (bool**)malloc(BW * sizeof(bool*));

    for (i = 0; i < BW; i++) {
        Dependence_Matrix[i] = (double*)calloc(BW, sizeof(double));
        Dependence[i] = (bool*)calloc(BW, sizeof(bool));
    }

    // Compute Dependence Matrix
    for (i = 1; i < BW; i++) { // Start from second bit
        for (j = 0; j < i; j++) { // Compare with previous bits
            Dependence_Matrix[i][j] = fabs(Conditional_Matrix[i][j] - Bit_Probabilities[i]);
            if (Dependence_Matrix[i][j] >= alpha) {
                Dependence[i][j] = true; // Mark bits as dependent
            }
        }
    }

    // Generate ranges for Link Vector
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
            if (ranges_count == 0 || i > ranges_max[ranges_count - 1]) {
                ranges_min[ranges_count] = i;
                ranges_max[ranges_count] = i;
                ranges_count++;
            } else {
                ranges_max[ranges_count - 1] = i;
            }
        }
    }

    // Merge overlapping ranges
    for (i = 0; i < ranges_count - 1; i++) {
        if (ranges_max[i] >= ranges_min[i + 1]) {
            ranges_max[i] = fmax(ranges_max[i], ranges_max[i + 1]);
            for (j = i + 1; j < ranges_count - 1; j++) {
                ranges_min[j] = ranges_min[j + 1];
                ranges_max[j] = ranges_max[j + 1];
            }
            ranges_count--;
            i--; // Recheck merged range
        }
    }

    // Construct Link Vector (LV)
    idx = 0;
    for (i = 0; i < ranges_count; i++) {
        for (j = ranges_min[i]; j <= ranges_max[i]; j++) {
            LV[idx++] = j;
        }
        if (i < ranges_count - 1 && ranges_min[i + 1] > ranges_max[i] + 1) {
            LV[idx++] = -1; // Separator for independent bits
        }
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
