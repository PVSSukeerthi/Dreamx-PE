#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Function prototypes
void Probability_to_PMF(const int *Indexes, int Indexes_len, const double *IPM, const double *Err, double Sig, double *PMF);

void Partial_PMF_Computation(int I, double A_p, double B_p, double C_in_p, double Sig,
                             double *C_out_p, double *PMF_00, double *PMF_01, double *PMF_10, double *PMF_11) {
    // Define all possible combinations of A, B, and C_in
    int A[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    int B[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    int C_in[8] = {0, 1, 0, 1, 0, 1, 0, 1};
    double A_v_p[8], B_v_p[8], C_in_v_p[8], IPM[8];
    double C_out[8], Err[8];

    // Compute individual probabilities
    for (int i = 0; i < 8; i++) {
        A_v_p[i] = (A[i] ? A_p : 1 - A_p);
        B_v_p[i] = (B[i] ? B_p : 1 - B_p);
        C_in_v_p[i] = (C_in[i] ? C_in_p : 1 - C_in_p);
        IPM[i] = A_v_p[i] * B_v_p[i] * C_in_v_p[i];
    }

    // Assign C_out and Err based on the configuration
    switch (I) {
        case 0: // Accurate Full Adder
            memcpy(C_out, (double[]){0, 0, 0, 1, 0, 1, 1, 1}, sizeof(C_out));
            memcpy(Err, (double[]){0, 0, 0, 0, 0, 0, 0, 0}, sizeof(Err));
            break;
        case 1: // Approximate Full Adder Type 1
            memcpy(C_out, (double[]){0, 0, 1, 1, 0, 1, 1, 1}, sizeof(C_out));
            memcpy(Err, (double[]){0, 0, 1, 0, -1, 0, 0, 0}, sizeof(Err));
            break;
        case 2: // Approximate Full Adder Type 2
            memcpy(C_out, (double[]){0, 0, 0, 1, 0, 1, 1, 1}, sizeof(C_out));
            memcpy(Err, (double[]){1, 0, 0, 0, 0, 0, 0, -1}, sizeof(Err));
            break;
        case 3: // Approximate Full Adder Type 3
            memcpy(C_out, (double[]){0, 0, 1, 1, 0, 1, 1, 1}, sizeof(C_out));
            memcpy(Err, (double[]){1, 0, 1, 0, 0, 0, 0, -1}, sizeof(Err));
            break;
        case 4: // Approximate Full Adder Type 4
            memcpy(C_out, (double[]){0, 0, 0, 0, 1, 1, 1, 1}, sizeof(C_out));
            memcpy(Err, (double[]){0, 0, -1, -1, 1, 0, 0, 0}, sizeof(Err));
            break;
        case 5: // Approximate Full Adder Type 5
            memcpy(C_out, (double[]){0, 0, 0, 0, 1, 1, 1, 1}, sizeof(C_out));
            memcpy(Err, (double[]){0, -1, 0, -1, 1, 0, 1, 0}, sizeof(Err));
            break;
        case 6: // Approximate Full Adder Type 6
            memcpy(C_out, (double[]){0, 1, 0, 1, 0, 1, 0, 1}, sizeof(C_out));
            memcpy(Err, (double[]){0, 2, 0, 0, 0, 0, -2, 0}, sizeof(Err));
            break;
        case 7: // Approximate Full Adder Type 7
            memcpy(C_out, (double[]){0, 0, 0, 1, 0, 1, 1, 1}, sizeof(C_out));
            memcpy(Err, (double[]){0, 0, 0, 1, 0, 1, 0, 0}, sizeof(Err));
            break;
    }

    // Compute carry-out probability
    *C_out_p = 0;
    for (int i = 0; i < 8; i++) {
        *C_out_p += C_out[i] * IPM[i];
    }

    // Compute partial error PMFs
    int Indexes[8];
    int count;

    // PMF_00: C_in = 0 and C_out = 0
    count = 0;
    for (int i = 0; i < 8; i++) {
        if (C_in[i] == 0 && C_out[i] == 0) Indexes[count++] = i;
    }
    Probability_to_PMF(Indexes, count, IPM, Err, Sig, PMF_00);

    // PMF_01: C_in = 0 and C_out = 1
    count = 0;
    for (int i = 0; i < 8; i++) {
        if (C_in[i] == 0 && C_out[i] == 1) Indexes[count++] = i;
    }
    Probability_to_PMF(Indexes, count, IPM, Err, Sig, PMF_01);

    // PMF_10: C_in = 1 and C_out = 0
    count = 0;
    for (int i = 0; i < 8; i++) {
        if (C_in[i] == 1 && C_out[i] == 0) Indexes[count++] = i;
    }
    Probability_to_PMF(Indexes, count, IPM, Err, Sig, PMF_10);

    // PMF_11: C_in = 1 and C_out = 1
    count = 0;
    for (int i = 0; i < 8; i++) {
        if (C_in[i] == 1 && C_out[i] == 1) Indexes[count++] = i;
    }
    Probability_to_PMF(Indexes, count, IPM, Err, Sig, PMF_11);
}

// Placeholder for Probability_to_PMF function
void Probability_to_PMF(const int *Indexes, int Indexes_len, const double *IPM, const double *Err, double Sig, double *PMF) {
    // Implement logic to compute the PMF for given indexes
    // (Placeholder: PMF computation logic should go here)
}
