#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_BITS 32 // Maximum bitwidth for inputs
#define MAX_VALUES (1 << MAX_BITS) // Maximum possible values for distribution

// Function to compute conditional probabilities
void conditionalProbabilityEstimation(double* Dist, int BW, int CL, double** Conditional_Matrix, double* Bit_Probabilities) {
    int i, j, k, l;
    unsigned int Dec_Values[MAX_VALUES]; // Array of decimal values
    int numValues = (1 << BW); // Total number of values based on bitwidth

    // Initialize decimal values
    for (i = 0; i < numValues; i++) {
        Dec_Values[i] = i;
    }

    // Compute bit-level probabilities
    for (i = 0; i < BW; i++) {
        double bitSum = 0.0;
        for (j = 0; j < numValues; j++) {
            if ((Dec_Values[j] & (1 << i)) != 0) { // Check if ith bit is 1
                bitSum += Dist[j];
            }
        }
        Bit_Probabilities[i] = bitSum;
    }

    // Compute conditional probabilities
    for (i = BW - 1; i >= 1; i--) { // Start from most significant bit
        for (j = i - 1; j >= fmax(0, i - CL); j--) { // Check coverage length
            for (k = 0; k <= 1; k++) { // Possible values for bit_i
                for (l = 0; l <= 1; l++) { // Possible values for bit_j
                    double numerator = 0.0, denominator = 0.0;
                    for (int m = 0; m < numValues; m++) {
                        int bit_i_val = (Dec_Values[m] & (1 << i)) != 0;
                        int bit_j_val = (Dec_Values[m] & (1 << j)) != 0;

                        if (bit_j_val == l) {
                            denominator += Dist[m];
                            if (bit_i_val == k) {
                                numerator += Dist[m];
                            }
                        }
                    }
                    if (denominator > 0) {
                        Conditional_Matrix[i][j] = numerator / denominator;
                    } else {
                        Conditional_Matrix[i][j] = 0.0; // Handle division by zero
                    }
                }
            }
        }
    }
}
