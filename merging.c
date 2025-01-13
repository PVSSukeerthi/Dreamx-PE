#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define MAX_BITS 32      // Maximum bitwidth
#define MAX_VALUES (1 << MAX_BITS) // Maximum possible input values

// Function Prototypes
void conditionalProbabilityEstimation(double* Dist, int BW, int CL, double** Conditional_Matrix, double* Bit_Probabilities);
void computeLinkVector(double** Conditional_Matrix, double* Bit_Probabilities, int BW, double alpha, int* LV, int* LV_size);
double** mergeTruthTables(double** T1, int T1_rows, int T1_cols, double** T2, int T2_rows, int T2_cols, int* T_rows, int* T_cols);

int main() {
    // Input Parameters
    int BW = 8;                 // Bitwidth of the input operands
    int CL = 2;                 // Coverage length for dependence checking
    double alpha = 0.1;         // Dependence threshold
    int T1_rows = 4, T1_cols = 3; // Example dimensions for T1 truth table
    int T2_rows = 8, T2_cols = 2; // Example dimensions for T2 truth table

    // Example probability distribution (uniform distribution)
    double Dist[MAX_VALUES];
    for (int i = 0; i < (1 << BW); i++) {
        Dist[i] = 1.0 / (1 << BW);
    }

    // Allocate memory for Conditional Matrix
    double** Conditional_Matrix = (double**)malloc(BW * sizeof(double*));
    for (int i = 0; i < BW; i++) {
        Conditional_Matrix[i] = (double*)calloc(BW, sizeof(double));
    }

    // Allocate memory for Bit Probabilities
    double* Bit_Probabilities = (double*)calloc(BW, sizeof(double));

    // Step 1: Compute Conditional Probabilities
    printf("Step 1: Conditional Probability Estimation\n");
    conditionalProbabilityEstimation(Dist, BW, CL, Conditional_Matrix, Bit_Probabilities);

    printf("Bit Probabilities:\n");
    for (int i = 0; i < BW; i++) {
        printf("P(bit_%d = 1) = %.4f\n", i, Bit_Probabilities[i]);
    }

    printf("\nConditional Matrix:\n");
    for (int i = 0; i < BW; i++) {
        for (int j = 0; j < BW; j++) {
            if (Conditional_Matrix[i][j] > 0) {
                printf("P(bit_%d | bit_%d) = %.4f\n", i, j, Conditional_Matrix[i][j]);
            }
        }
    }

    // Step 2: Compute Link Vector
    int LV[MAX_BITS];
    int LV_size = 0;

    printf("\nStep 2: Link Vector Computation\n");
    computeLinkVector(Conditional_Matrix, Bit_Probabilities, BW, alpha, LV, &LV_size);

    printf("Link Vector (LV):\n");
    for (int i = 0; i < LV_size; i++) {
        if (LV[i] == -1) {
            printf(" | "); // Separator for independent groups
        } else {
            printf("%d ", LV[i]);
        }
    }
    printf("\n");

    // Step 3: Merge Truth Tables
    printf("\nStep 3: Merging Truth Tables\n");

    // Example truth tables (values are illustrative)
    double T1[4][3] = {
        {0, 0, 0},
        {0, 1, 1},
        {1, 0, 0},
        {1, 1, 1}
    };
    double T2[8][2] = {
        {0, 0},
        {1, 1},
        {0, 1},
        {1, 0},
        {0, 0},
        {1, 1},
        {0, 1},
        {1, 0}
    };

    // Convert static tables to dynamically allocated arrays
    double** T1_dyn = (double**)malloc(T1_rows * sizeof(double*));
    double** T2_dyn = (double**)malloc(T2_rows * sizeof(double*));
    for (int i = 0; i < T1_rows; i++) {
        T1_dyn[i] = (double*)malloc(T1_cols * sizeof(double));
        for (int j = 0; j < T1_cols; j++) {
            T1_dyn[i][j] = T1[i][j];
        }
    }
    for (int i = 0; i < T2_rows; i++) {
        T2_dyn[i] = (double*)malloc(T2_cols * sizeof(double));
        for (int j = 0; j < T2_cols; j++) {
            T2_dyn[i][j] = T2[i][j];
        }
    }

    // Merge the truth tables
    int T_rows, T_cols;
    double** T = mergeTruthTables(T1_dyn, T1_rows, T1_cols, T2_dyn, T2_rows, T2_cols, &T_rows, &T_cols);

    printf("Merged Truth Table:\n");
    for (int i = 0; i < T_rows; i++) {
        for (int j = 0; j < T_cols; j++) {
            printf("%.1f ", T[i][j]);
        }
        printf("\n");
    }

    // Free allocated memory
    for (int i = 0; i < BW; i++) free(Conditional_Matrix[i]);
    for (int i = 0; i < T1_rows; i++) free(T1_dyn[i]);
    for (int i = 0; i < T2_rows; i++) free(T2_dyn[i]);
    for (int i = 0; i < T_rows; i++) free(T[i]);

    free(Conditional_Matrix);
    free(Bit_Probabilities);
    free(T1_dyn);
    free(T2_dyn);
    free(T);

    return 0;
}
