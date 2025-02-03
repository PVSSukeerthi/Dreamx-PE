#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "conditionalP.h"
#include "LV.h"
#include "merging.h"

#define MAX_BITS 32

// Function to process and merge truth tables based on Link Vector groups
void processAndMerge(int ***outputs, int numAdders, int rows, int output_cols, 
                     int **accurateAdder, int *LV, int LV_size, int *configuration) {
    
    printf("\nProcessing Truth Tables Based on Link Vector Groups:\n");

    int start = 0;
    while (start < LV_size) {
        int groupSize = 0;

        // **Identify the size of the current group (skip -1 separators)**
        while (start + groupSize < LV_size && LV[start + groupSize] != -1) {
            groupSize++;
        }

        // **Store the first bit in the group for correct indexing**
        int groupStartIndex = start;

        // **Skip the separator (-1)**
        start += groupSize + 1;

        // **Skip empty groups**
        if (groupSize == 0) continue;

        printf("\nProcessing a Group of Size: %d\n", groupSize);

        // **✅ Extract the correct subset from `configuration`**
        int configurationSubset[groupSize];  
        int subsetIndex = 0;  

        for (int i = 0; i < groupSize; i++) {  
            int bitIndex = LV[groupStartIndex + i];  // ✅ **Fixed index tracking**
            
            // **Skip -1 values (separators)**
            if (bitIndex == -1) continue;

            // **Ensure valid bitIndex**
            if (bitIndex < 0 || bitIndex >= LV_size) {
                printf("Error: Invalid bit index %d in configuration\n", bitIndex);
                continue;
            }

            // **Ensure valid adder index**
            if (configuration[bitIndex] < 0 || configuration[bitIndex] >= numAdders) {
                printf("Error: Invalid adder index %d for bit %d in configuration\n", configuration[bitIndex], bitIndex);
                continue;
            }

            configurationSubset[subsetIndex++] = configuration[bitIndex];  // ✅ **Fix: Include first bit properly**
        }

        // **Ensure at least one valid entry in `configurationSubset`**
        if (subsetIndex == 0) {
            printf("Warning: No valid adders found for this group, skipping merge.\n");
            continue;
        }

        // **Print configurationSubset for debugging**
        printf("Using Adders: ");
        for (int i = 0; i < subsetIndex; i++) {
            printf("%d ", configurationSubset[i]);
        }
        printf("\n");

        // **Ensure `outputs` is not NULL**
        if (!outputs) {
            printf("Error: Outputs array is NULL\n");
            return;
        }

        // **Ensure `outputs[configurationSubset[i]]` exists before accessing**
        for (int i = 0; i < subsetIndex; i++) {
            if (!outputs[configurationSubset[i]]) {
                printf("Error: outputs[%d] is NULL\n", configurationSubset[i]);
                return;
            }
        }

        int **finalOutputs = NULL;
        int finalRows, finalOutputCols;
        mergeMultipleTruthTables(outputs, numAdders, rows, output_cols, 
                                 configurationSubset, subsetIndex, &finalOutputs, 
                                 &finalRows, &finalOutputCols);

        // **Ensure finalOutputs is not NULL before proceeding**
        if (!finalOutputs) {
            printf("Error: mergeMultipleTruthTables returned NULL\n");
            return;
        }

        // **Generate accurate merged truth table**
        int **finalAccurateOutputs = NULL;
        int finalAccurateRows, finalAccurateOutputCols;
        mergeAccurateAdder(accurateAdder, rows, output_cols, subsetIndex, 
                           &finalAccurateOutputs, &finalAccurateRows, &finalAccurateOutputCols);

        // **Ensure finalAccurateOutputs is not NULL before proceeding**
        if (!finalAccurateOutputs) {
            printf("Error: mergeAccurateAdder returned NULL\n");
            return;
        }

        // **Print Headers**
        char output_headers[subsetIndex + 1][10];
        for (int i = 0; i < subsetIndex; i++) {
            sprintf(output_headers[i], "Sum%d", i);
        }
        sprintf(output_headers[subsetIndex], "Cout");

        printf("\nMerged Approximate Truth Table:\n");
        if(groupSize==1)
        printTruthTable(finalOutputs, finalRows, finalOutputCols, output_headers);

        printf("\nMerged Accurate Truth Table:\n");
        // printTruthTable(finalAccurateOutputs, finalAccurateRows, finalAccurateOutputCols, output_headers);

        // **Free allocated memory**
        for (int i = 0; i < finalRows; i++) free(finalOutputs[i]);
        free(finalOutputs);

        for (int i = 0; i < finalAccurateRows; i++) free(finalAccurateOutputs[i]);
        free(finalAccurateOutputs);
    }
}


int main() {
    int BW = 8;  // Bit Width
    double alpha = 0.5;  // Dependence threshold

    // Example configuration array (user-specified adders for each bit position)
    int configuration[MAX_BITS] = {5, 2, 3, 1, 4, 5, 3, 6};

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
