#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "conditionalP.h"
#include "LV.h"
#include "merging.h"

#define MAX_BITS 32



int main() {
    int BW = 8;  // Bit Width
    double alpha = 0.8;  // Dependence threshold

    double MSE = 0.0, MED = 0.0;
    
    // Output variables for Exhaustive Simulations
    double MSE_Ex = 0.0, MED_Ex = 0.0;

    // Example configuration array (user-specified adders for each bit position)
    int configuration[MAX_BITS] = {0,1, 2, 3, 4,5,6,7};

    // Allocate probability distribution
    double *dist = (double *)malloc((1 << BW) * sizeof(double));
    if (!dist) {
        printf("Memory allocation failed for probability distribution!\n");
        return 1;
    }
    for (int i = 0; i < (1 << BW); i++) {
        dist[i] = 1.0/ (1 << BW);  // Uniform probability distribution
    }

    double ** CPM_00= (double **)malloc(BW * sizeof(double *));
        for (int i = 0; i < BW; i++) {
            CPM_00[i] = (double *)malloc(BW* sizeof(double));
        }
        double ** CPM_01= (double **)malloc(BW * sizeof(double *));
        for (int i = 0; i < BW; i++) {
            CPM_01[i] = (double *)malloc(BW* sizeof(double));
        }
        double ** CPM_10= (double **)malloc(BW * sizeof(double *));
        for (int i = 0; i < BW; i++) {
            CPM_10[i] = (double *)malloc(BW* sizeof(double));
        }
        double ** CPM_11= (double **)malloc(BW * sizeof(double *));
        for (int i = 0; i < BW; i++) {
            CPM_11[i] = (double *)malloc(BW* sizeof(double));
        }
    // if (!Conditional_Matrix) {
    //     printf("Memory allocation failed for Conditional_Matrix!\n");
    //     return 1;
    // }
    // for (int i = 0; i < BW; i++) {
    //     Conditional_Matrix[i] = (double *)malloc(BW * sizeof(double));
    //     if (!Conditional_Matrix[i]) {
    //         printf("Memory allocation failed for Conditional_Matrix[%d]!\n", i);
    //         return 1;
    //     }
    //     // Initialize to prevent garbage values
    //     for (int j = 0; j < BW; j++) {
    //         Conditional_Matrix[i][j] = 0.0;
    //     }
    // }
    //Probability_A_bits, Probability_B_bits, Probability_C_in,
    double *Probability_A_bits = (double *)malloc(BW * sizeof(double));
    if (!Probability_A_bits) {
        printf("Memory allocation failed for Bit_Probabilities!\n");
        return 1;
    }

    // Compute Conditional Probability Matrix and Bit Probabilities
    calculate_all_CPMs(dist, BW, 2, CPM_00, CPM_01, CPM_10, CPM_11, Probability_A_bits);

    double *Probability_B_bits = (double *)malloc(BW * sizeof(double));
        if (!Probability_B_bits) {
            printf("Memory allocation failed for Bit_Probabilities!\n");
            return 1;
        }

        // Compute Conditional Probability Matrix and Bit Probabilities
    calculate_all_CPMs(dist, BW, 2, CPM_00, CPM_01, CPM_10, CPM_11, Probability_B_bits);

    double Probability_C_in = 0.5;
    printf("CPM_00:\n");
    for (int i = 0; i < BW; i++) {
        for (int j = 0; j < BW; j++) {
            printf("%0.4lf ", CPM_00[i][j]);
        }
        printf("\n");
    }
    printf("CPM_01:\n");
    for (int i = 0; i < BW; i++) {
        for (int j = 0; j < BW; j++) {
            printf("%0.4lf ", CPM_01[i][j]);
        }
        printf("\n");
    }
    printf("CPM_10:\n");
    for (int i = 0; i < BW; i++) {
        for (int j = 0; j < BW; j++) {
            printf("%0.4lf ", CPM_10[i][j]);
        }
        printf("\n");
    }
    printf("CPM_11:\n");
    for (int i = 0; i < BW; i++) {
        for (int j = 0; j < BW; j++) {
            printf("%0.4lf ", CPM_11[i][j]);
        }
        printf("\n");
    }
    // Compute Link Vector
    int *LV = (int *)malloc(MAX_BITS * sizeof(int));
    if (!LV) {
        printf("Memory allocation failed for Link Vector!\n");
        return 1;
    }
    int LV_size = 0;
    computeLinkVector(CPM_00, CPM_01, CPM_10, CPM_11, Probability_B_bits , BW, alpha, LV, &LV_size);

    // Print Link Vector
    printf("\nLink Vector: ");
    for (int i = 0; i < LV_size; i++) {
        printf("%d ", LV[i]);
    }
    printf("\n");

    // Initialize truth tables
    int ***outputs = NULL;
    int **accurateAdder = NULL;
    int rows, output_cols, numAdders;
    initializeTruthTables(&outputs, &rows, &output_cols, &numAdders, &accurateAdder);

    // Process and merge truth tables based on Link Vector groups
    // processAndMerge(outputs, numAdders, rows, output_cols, accurateAdder, LV, LV_size, configuration);

    processAndMerge(outputs,numAdders,rows,output_cols, 
                     accurateAdder, LV,  LV_size,configuration ,Probability_A_bits, Probability_B_bits, Probability_C_in, BW, &MSE,&MED, CPM_00, CPM_01, CPM_10, CPM_11) ;
   
        printf("PEMACx Results:\n");
    printf("MSE: %f\n", MSE);
    printf("MED: %f\n\n", MED);
    
    // Free allocated memory
    free(dist);
    free(Probability_A_bits);
    free(Probability_B_bits);
    free(LV);
    for (int i = 0; i < BW; i++) {
        free(CPM_00[i]);
        free(CPM_01[i]);
        free(CPM_10[i]);
        free(CPM_11[i]);
    }
    free(CPM_00);
    free(CPM_01);
    free(CPM_10);
    free(CPM_11);

    for (int k = 0; k < numAdders; k++) {
        for (int i = 0; i < rows; i++) free(outputs[k][i]);
        free(outputs[k]);
    }
    free(outputs);

    for (int i = 0; i < rows; i++) free(accurateAdder[i]);
    free(accurateAdder);

    return 0;
}
