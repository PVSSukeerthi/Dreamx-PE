#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Function prototypes
void Validate_Inputs(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length);
void Partial_PMF_Computation(int I, double A_p, double B_p, double C_in_p, double Sig, double *C_out_p, double *PMF_00, double *PMF_01, double *PMF_10, double *PMF_11);
void Convolve(const double *x, int x_len, const double *y, int y_len, double *result);

void PEMACx(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length, double *MSE, double *MED) {
    // Validate inputs
    Validate_Inputs(Adder_Config, Probability_A_bits, Probability_B_bits, Probability_C_in, length);

    double PMF_0[100] = {0}, PMF_1[100] = {0}, PMF_0_new[200] = {0};
    double PMF_00[10] = {0}, PMF_01[10] = {0}, PMF_10[10] = {0}, PMF_11[10] = {0};
    double C_out_p = 0, C_in_p = Probability_C_in;
    double PMF[200] = {0};
    int Min_val = -3;
    int PMF_0_len = 1, PMF_1_len = 1;

    for (int i = 0; i < length; i++) {
        double A_p = Probability_A_bits[i];
        double B_p = Probability_B_bits[i];
        double Sig = pow(2, i);
        int I = Adder_Config[i];

        // Compute partial PMFs
        Partial_PMF_Computation(I, A_p, B_p, C_in_p, Sig, &C_out_p, PMF_00, PMF_01, PMF_10, PMF_11);

        if (i == 0) {
            // Initialize PMF_0 and PMF_1
            for (int j = 0; j < 10; j++) {
                PMF_0[j] = PMF_00[j] + PMF_10[j];
                PMF_1[j] = PMF_01[j] + PMF_11[j];
            }
        } else {
            // Convolve PMFs
            Convolve(PMF_00, 10, PMF_0, PMF_0_len, PMF_0_new);
            Convolve(PMF_10, 10, PMF_1, PMF_1_len, PMF_0);
            Convolve(PMF_01, 10, PMF_0, PMF_0_len, PMF_1);
            Convolve(PMF_11, 10, PMF_1, PMF_1_len, PMF_1);
            Min_val += (-3) * Sig;
        }

        // Normalize PMFs
        if (i != length - 1) {
            double sum_PMF_0 = 0, sum_PMF_1 = 0;
            for (int j = 0; j < PMF_0_len; j++) sum_PMF_0 += PMF_0[j];
            for (int j = 0; j < PMF_1_len; j++) sum_PMF_1 += PMF_1[j];
            for (int j = 0; j < PMF_0_len; j++) PMF_0[j] /= sum_PMF_0;
            for (int j = 0; j < PMF_1_len; j++) PMF_1[j] /= sum_PMF_1;
        }
    }

    // Combine PMF_0 and PMF_1
    for (int i = 0; i < PMF_0_len; i++) {
        PMF[i] = PMF_0[i] + PMF_1[i];
    }

    // Compute MSE and MED
    double mse = 0, med = 0;
    for (int i = Min_val; i < Min_val + PMF_0_len; i++) {
        mse += pow(i, 2) * PMF[i - Min_val];
        med += fabs(i) * PMF[i - Min_val];
    }
    *MSE = mse;
    *MED = med;
}

// Placeholder for validation function
void Validate_Inputs(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length) {
    // Add input validation logic
}

// Placeholder for partial PMF computation
void Partial_PMF_Computation(int I, double A_p, double B_p, double C_in_p, double Sig, double *C_out_p, double *PMF_00, double *PMF_01, double *PMF_10, double *PMF_11) {
    // Add logic for PMF computation
}

// Function to perform convolution
void Convolve(const double *x, int x_len, const double *y, int y_len, double *result) {
    memset(result, 0, (x_len + y_len - 1) * sizeof(double));
    for (int i = 0; i < x_len; i++) {
        for (int j = 0; j < y_len; j++) {
            result[i + j] += x[i] * y[j];
        }
    }
}
