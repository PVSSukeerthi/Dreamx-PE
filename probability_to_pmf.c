#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void Probability_to_PMF(const int *Indexes, int Indexes_len, const double *IPM_full, const double *Err_full, double Sig, double *PMF) {
    // Initialize PMF array
    int PMF_size = (int)(6 * Sig + 1);
    memset(PMF, 0, PMF_size * sizeof(double));

    // Iterate over Indexes to compute PMF
    for (int i = 0; i < Indexes_len; i++) {
        int idx = Indexes[i];
        double Err = Err_full[idx];      // Get corresponding error value
        double IPM = IPM_full[idx];      // Get corresponding probability

        // Calculate PMF index: (Err + 4) * Sig - Sig + 1
        int PMF_index = (int)((Err + 4) * Sig - Sig + 1);

        // Update PMF
        if (PMF_index >= 0 && PMF_index < PMF_size) {
            PMF[PMF_index] += IPM;
        }
    }
}

int main() {
    // Example usage
    int Indexes[] = {0, 1, 2};
    int Indexes_len = 3;
    double IPM_full[] = {0.1, 0.2, 0.3, 0.4}; // Full IPM array
    double Err_full[] = {-1, 0, 1, 2};        // Full Err array
    double Sig = 2;                           // Signal weight

    int PMF_size = (int)(6 * Sig + 1);
    double *PMF = (double *)calloc(PMF_size, sizeof(double)); // Allocate PMF array

    Probability_to_PMF(Indexes, Indexes_len, IPM_full, Err_full, Sig, PMF);

    // Print PMF
    for (int i = 0; i < PMF_size; i++) {
        printf("PMF[%d] = %.2f\n", i, PMF[i]);
    }

    free(PMF); // Free allocated memory
    return 0;
}
