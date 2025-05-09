#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "conditionalP.h"
#include "LV.h"
#include "merging.h"

#define MAX_BITS 32




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Function prototypes
// void Validate_Inputs(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length);
void Partial_PMF_Computation(int I, double A_p, double B_p, double C_in_p, double Sig, double *C_out_p, double *PMF_00, double *PMF_01, double *PMF_10, double *PMF_11);
// void Convolve(const double *x, int x_len, const double *y, int y_len, double *result);
void Probability_to_PMF(const int *Indexes, int Indexes_len, const double *IPM_full, const double *Err_full, double Sig, double *PMF) ;
void Convolve(const double *x, int x_len, const double *y, int y_len, double **result, int *result_len) ;
void Validate_Inputs(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length);
double* Convolution(double* a, int len_a, double* b, int len_b) ;


void PEMACx(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length, double *MSE, double *MED) {
    // Validate inputs
 Validate_Inputs(Adder_Config,Probability_A_bits, Probability_B_bits, Probability_C_in,  length);

 // Initialize variables
    int PMF_0_length = 7;  // Initial size

    double* PMF_0 = (double*)calloc(PMF_0_length, sizeof(double));  // Adjust size as needed
    double* PMF_1 = (double*)calloc(PMF_0_length, sizeof(double));  // Adjust size as needed
    
    
        double C_out_p;
        double Min_val = -3;
        double sig;
        double final_length;
        int new_length =PMF_0_length;
    
    // Loop through the bits in the adder configuration
    for (int i = 0; i < length; i++) {
        sig = pow(2, i);
         int PMF_length = (int)(6 * sig + 1);
        double* PMF_00 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
        double* PMF_01 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
        double* PMF_10 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
         double* PMF_11 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
        double A_p = Probability_A_bits[i];
        double B_p = Probability_B_bits[i];
        int I = Adder_Config[i];
        

        double C_in_p;
        if (i == 0) {
            C_in_p = Probability_C_in;
        } else {
            C_in_p = C_out_p;
        }

        // Compute the partial PMF values
        Partial_PMF_Computation(I, A_p, B_p, C_in_p, sig, &C_out_p, PMF_00, PMF_01, PMF_10, PMF_11);

        // Update PMFs for current bit
        if (i == 0) {
            for (int j = 0; j < PMF_length; j++) {  // Adjust size as needed
                PMF_0[j] = PMF_00[j] + PMF_10[j];  // Convolution would be done here
                PMF_1[j] = PMF_01[j] + PMF_11[j];  // Convolution would be done here

                Min_val = -3;
            }
        }
    else{        // Normalize PMFs if not the last bit

       
        double* temp_PMF_1 = Convolution(PMF_00, PMF_length, PMF_0, new_length);
        double* temp_PMF_2 = Convolution(PMF_10, PMF_length, PMF_1, new_length);
        double* temp_PMF_3 = Convolution(PMF_01, PMF_length, PMF_0, new_length);
        double* temp_PMF_4 = Convolution(PMF_11,  PMF_length, PMF_1, new_length);
                                    

        new_length = (PMF_length + new_length )- 1;
        final_length = new_length;

        PMF_0 = (double*)realloc(PMF_0, new_length * sizeof(double));
        PMF_1 = (double*)realloc(PMF_1, new_length * sizeof(double));
        if (!PMF_0 || !PMF_1) {
            printf("Memory allocation failed during reallocation!\n");
            exit(1);
        }
                              
        
    //Combine results to update PMF_1
        for (int i = 0; i < new_length; i++) {
            PMF_0[i] = temp_PMF_1[i] + temp_PMF_2[i];
            PMF_1[i] = temp_PMF_3[i] + temp_PMF_4[i];
        }
        

        Min_val = Min_val + ((-3)*sig);

       

    }
     if (i != length - 1) {
            double sum_0 = 0, sum_1 = 0;
            
            for (int j = 0; j < new_length; j++) {  // Adjust size as needed
                sum_0 += PMF_0[j];
                sum_1 += PMF_1[j];
            }
            for (int j = 0; j < new_length; j++) {  // Adjust size as needed
                PMF_0[j] /= sum_0;
                PMF_1[j] /= sum_1;
            }
    }
    free(PMF_00);
    free(PMF_01);
    free(PMF_10);
    free(PMF_11);
    }

    // Combine the final PMFs
    double* PMF = (double*)calloc(final_length, sizeof(double)); // Adjust size as needed
    for (int i = 0; i < final_length; i++) {  // Adjust size as needed
        PMF[i] = PMF_0[i] + PMF_1[i];

    }

    // Compute MSE
    *MSE = 0;
    for (int i = 0; i <final_length ; i++) {  // Adjust size as needed
        *MSE += (pow(Min_val + i, 2)) * PMF[i];
    }

    // Compute MED
    *MED = 0;
    for (int i = 0; i < final_length; i++) {  // Adjust size as needed
        *MED += fabs(Min_val + i) * PMF[i];
    }

    // Free dynamically allocated memory
   
    free(PMF);
}





void Validate_Inputs(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length) {
    if (length <= 0) {
        fprintf(stderr, "Error: Invalid length.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < length; i++) {
        if (Adder_Config[i] < 0 || Adder_Config[i] > 1000) {
            fprintf(stderr, "Error: Invalid adder configuration.\n");
            exit(EXIT_FAILURE);
        }
        if (Probability_A_bits[i] < 0 || Probability_A_bits[i] > 1 || Probability_B_bits[i] < 0 || Probability_B_bits[i] > 1) {
            fprintf(stderr, "Error: Invalid probability value.\n");
            exit(EXIT_FAILURE);
        }
    }
    if (Probability_C_in < 0 || Probability_C_in > 1) {
        fprintf(stderr, "Error: Invalid C_in probability.\n");
        exit(EXIT_FAILURE);
    }
}



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
          A_v_p[i] = (A[i] * A_p) + ((!A[i]) * (1 - A_p));  // Bitwise NOT operation on A
        B_v_p[i] = (B[i] * B_p) + ((!B[i]) * (1 - B_p));  // Bitwise NOT operation on B
        C_in_v_p[i] = (C_in[i] * C_in_p) + ((!C_in[i]) * (1 - C_in_p)); 
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
        case 7: // Approximate Full Adder Type 1000
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
void Probability_to_PMF(const int *Indexes, int Indexes_len, const double *IPM_full, const double *Err_full, double Sig, double *PMF) {
    
    int PMF_size = (int)(6 * Sig + 1);

    for (int i = 0; i < PMF_size; i++) {
        PMF[i] = 0.0; // Explicit zero initialization
    }

    // Extract Err and IPM values based on Indexes
    double Err[Indexes_len], IPM[Indexes_len];
    for (int i = 0; i < Indexes_len; i++) {
        Err[i] = Err_full[Indexes[i]];
        IPM[i] = IPM_full[Indexes[i]];
    }

    // Populate PMF array
    for (int i = 0; i < Indexes_len; i++) {
        int PMF_index = (int)(((Err[i] + 4) * Sig )- Sig ); // Compute index
        if (PMF_index >= 0 && PMF_index < PMF_size) { // Check bounds
            PMF[PMF_index] += IPM[i]; // Add corresponding probability
        }
    }
}

// Function to perform convolution

double* Convolution(double* a, int len_a, double* b, int len_b) {
    int len_result = len_a + len_b - 1;  // Length of the result array
    double* result = (double*)malloc(len_result * sizeof(double));  // Allocate memory for the result array

    // Initialize the result array to zero
    for (int i = 0; i < len_result; i++) {
        result[i] = 0.0;
    }

    // Perform the convolution
    for (int i = 0; i < len_result; i++) {
        for (int j = 0; j < len_a; j++) {
            if (i - j >= 0 && i - j < len_b) {
                result[i] += a[j] * b[i - j];
            }
        }
    }

    return result;  // Return the result array
}








int main() {
    int BW = 5;  // Bit Width
    double alpha = 0.2;  // Dependence threshold

    // Example configuration array (user-specified adders for each bit position)
    int configuration[MAX_BITS] = {5, 2, 3, 0, 1};

    // Allocate probability distribution
    double *dist = (double *)malloc((1 << BW) * sizeof(double));
    if (!dist) {
        printf("Memory allocation failed for probability distribution!\n");
        return 1;
    }
    for (int i = 0; i < (1 << BW); i++) {
        dist[i] = 1.0 / (1 << BW);  // Uniform probability distribution
    }

    // Allocate matrices
    double **Conditional_Matrix = (double **)malloc(BW * sizeof(double *));
    if (!Conditional_Matrix) {
        printf("Memory allocation failed for Conditional_Matrix!\n");
        return 1;
    }
    for (int i = 0; i < BW; i++) {
        Conditional_Matrix[i] = (double *)malloc(BW * sizeof(double));
        if (!Conditional_Matrix[i]) {
            printf("Memory allocation failed for Conditional_Matrix[%d]!\n", i);
            return 1;
        }
        // Initialize to prevent garbage values
        for (int j = 0; j < BW; j++) {
            Conditional_Matrix[i][j] = 0.0;
        }
    }

    double *Bit_Probabilities = (double *)malloc(BW * sizeof(double));
    if (!Bit_Probabilities) {
        printf("Memory allocation failed for Bit_Probabilities!\n");
        return 1;
    }

    // Compute Conditional Probability Matrix and Bit Probabilities
    calculate_conditional_matrix_and_bit_probabilities(dist, BW, 2, Conditional_Matrix, Bit_Probabilities);


    // Compute Link Vector
    int *LV = (int *)malloc(MAX_BITS * sizeof(int));
    if (!LV) {
        printf("Memory allocation failed for Link Vector!\n");
        return 1;
    }
    int LV_size = 0;
    computeLinkVector(Conditional_Matrix, Bit_Probabilities, BW, alpha, LV, &LV_size);

    // Print Link Vector
    printf("\nLink Vector: ");
    for (int i = 0; i < LV_size; i++) {
        printf("%d ", LV[i]);
    }
    printf("\n");

    // Example setup
    int length = 5; // Length of the configuration and probability arrays
    double Probability_C_in = 0.5;
    // Output variables for PEMACx
    double MSE = 0.0, MED = 0.0;
    
    // Output variables for Exhaustive Simulations
    double MSE_Ex = 0.0, MED_Ex = 0.0;

    // Call PEMACx function
    PEMACx(configuration, Bit_Probabilities, Bit_Probabilities,Probability_C_in, length, &MSE, &MED);


    printf("PEMACx Results:\n");
    printf("MSE: %f\n", MSE);
    printf("MED: %f\n\n", MED);






    // Initialize truth tables
    int ***outputs = NULL;
    int **accurateAdder = NULL;
    int rows, output_cols, numAdders;
    initializeTruthTables(&outputs, &rows, &output_cols, &numAdders, &accurateAdder);

    // Process and merge truth tables based on Link Vector groups
    processAndMerge(outputs, numAdders, rows, output_cols, accurateAdder, LV, LV_size, configuration);

    // Free allocated memory
    free(dist);
    free(Bit_Probabilities);
    free(LV);
    for (int i = 0; i < BW; i++) free(Conditional_Matrix[i]);
    free(Conditional_Matrix);
    
    for (int k = 0; k < numAdders; k++) {
        for (int i = 0; i < rows; i++) free(outputs[k][i]);
        free(outputs[k]);
    }
    free(outputs);

    for (int i = 0; i < rows; i++) free(accurateAdder[i]);
    free(accurateAdder);

    return 0;
}
