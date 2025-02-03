#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to merge two truth tables
void mergeTwoTruthTables(int **T1_outputs, int T1_rows, int T1_output_cols,
                         int **T2_outputs, int T2_rows, int T2_output_cols,
                         int ***T_outputs, int *T_rows, int *T_output_cols) {
    // Calculate the dimensions of the merged truth table
    *T_rows = (T1_rows * T2_rows) / 2;
    *T_output_cols = T1_output_cols + T2_output_cols - 1;

    // Allocate memory for the merged outputs and initialize to zero
    *T_outputs = (int **)malloc((*T_rows) * sizeof(int *));
    for (int i = 0; i < *T_rows; i++) {
        (*T_outputs)[i] = (int *)calloc(*T_output_cols, sizeof(int)); // Use calloc for zero initialization
    }

    // Step 1: Copy T1 outputs into merged outputs (excluding the last column)
    for (int i = 0; i < *T_rows; i++) {
        for (int j = 0; j < T1_output_cols - 1; j++) {
            (*T_outputs)[i][j] = T1_outputs[i % T1_rows][j];
        }
    }

    // Step 2: Generate indexes for T2 rows based on T1's Cout and repeating logic
    int *indexes = (int *)malloc((*T_rows) * sizeof(int));
    for (int i = 0; i < *T_rows; i++) {
        indexes[i] = T1_outputs[i % T1_rows][T1_output_cols - 1] + (i / T1_rows) * 2;
    }

    // Step 3: Append T2 outputs to merged outputs using the calculated indexes
    for (int i = 0; i < *T_rows; i++) {
        for (int j = 0; j < T2_output_cols; j++) {
            (*T_outputs)[i][T1_output_cols - 1 + j] = T2_outputs[indexes[i]][j];
        }
    }

    // Free the allocated memory for indexes
    free(indexes);
}

// Function to merge multiple truth tables based on the configuration
void mergeMultipleTruthTables(int ***outputs, int numAdders, int rows, int output_cols,
                              int *configuration, int configSize, int ***finalOutputs,
                              int *finalRows, int *finalOutputCols) {
    int **mergedOutputs = NULL;
    int mergedRows = rows;
    int mergedOutputCols = output_cols;

    // Allocate initial merged table for safety
    mergedOutputs = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++) {
        mergedOutputs[i] = (int *)malloc(output_cols * sizeof(int));
        memcpy(mergedOutputs[i], outputs[configuration[0]][i], output_cols * sizeof(int));
    }

    for (int i = 1; i < configSize; i++) {
        int **tempMergedOutputs = NULL;
        int tempMergedRows, tempMergedOutputCols;

        // Merge the current merged table with the next adder in the configuration
        mergeTwoTruthTables(mergedOutputs, mergedRows, mergedOutputCols,
                            outputs[configuration[i]], rows, output_cols,
                            &tempMergedOutputs, &tempMergedRows, &tempMergedOutputCols);

        // Free the previous mergedOutputs
        for (int j = 0; j < mergedRows; j++) {
            free(mergedOutputs[j]);
        }
        free(mergedOutputs);

        mergedOutputs = tempMergedOutputs;
        mergedRows = tempMergedRows;
        mergedOutputCols = tempMergedOutputCols;
    }

    *finalOutputs = mergedOutputs;
    *finalRows = mergedRows;
    *finalOutputCols = mergedOutputCols;
}

// Function to merge accurate adders repeatedly
void mergeAccurateAdder(int **accurateAdder, int rows, int output_cols, int configSize,
                        int ***finalOutputs, int *finalRows, int *finalOutputCols) {
    int **mergedOutputs = NULL;
    int mergedRows = rows;
    int mergedOutputCols = output_cols;

    mergedOutputs = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++) {
        mergedOutputs[i] = (int *)malloc(output_cols * sizeof(int));
        memcpy(mergedOutputs[i], accurateAdder[i], output_cols * sizeof(int));
    }

    for (int i = 1; i < configSize; i++) {
        int **tempOutputs = NULL;
        int tempRows, tempOutputCols;

        mergeTwoTruthTables(mergedOutputs, mergedRows, mergedOutputCols,
                            accurateAdder, rows, output_cols,
                            &tempOutputs, &tempRows, &tempOutputCols);

        for (int j = 0; j < mergedRows; j++) {
            free(mergedOutputs[j]);
        }
        free(mergedOutputs);

        mergedOutputs = tempOutputs;
        mergedRows = tempRows;
        mergedOutputCols = tempOutputCols;
    }

    *finalOutputs = mergedOutputs;
    *finalRows = mergedRows;
    *finalOutputCols = mergedOutputCols;
}


// Utility function to print a truth table
void printTruthTable(int **outputs, int rows, int output_cols, char output_headers[][10]) {
    // Print output column headers
    for (int j = 0; j < output_cols; j++) {
        printf("%s\t", output_headers[j]);
    }
    printf("\n");

    // Print table rows
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < output_cols; j++) {
            printf("%d\t", outputs[i][j]);
        }
        printf("\n");
    }
}

// Hardcoded truth tables for approximate and accurate adders
void initializeTruthTables(int ****outputs, int *rows, int *output_cols, int *numAdders,
                           int ***accurateAdder) {
    *numAdders = 7;
    *rows = 8;
    *output_cols = 2;

    *outputs = (int ***)malloc((*numAdders) * sizeof(int **));
    for (int k = 0; k < *numAdders; k++) {
        (*outputs)[k] = (int **)malloc((*rows) * sizeof(int *));
        for (int i = 0; i < *rows; i++) {
            (*outputs)[k][i] = (int *)malloc((*output_cols) * sizeof(int));
        }
    }

    // Hardcoded outputs for approximate adders
    int hardcoded_outputs[7][8][2] = {
        {{0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {0, 1}, {1, 1}},
        {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 1}, {1, 1}, {1, 1}, {0, 0}},
        {{0, 0}, {1, 0}, {1, 0}, {1, 1}, {0, 1}, {1, 1}, {1, 1}, {0, 0}},
        {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {1, 0}, {0, 1}, {1, 1}, {0, 0}},
        {{0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {1, 1}},
        {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {0, 1}, {1, 1}, {0, 0}, {1, 0}},
        {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {1, 0}, {0, 1}, {1, 1}, {1, 1}}};

    for (int k = 0; k < 7; k++) {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 2; j++) {
                (*outputs)[k][i][j] = hardcoded_outputs[k][i][j];
            }
        }
    }

    // Hardcoded accurate adder
    *accurateAdder = (int **)malloc((*rows) * sizeof(int *));
    int accurate[8][2] = {
        {0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {0, 1}, {1, 1}};

    for (int i = 0; i < *rows; i++) {
        (*accurateAdder)[i] = (int *)malloc((*output_cols) * sizeof(int));
        for (int j = 0; j < *output_cols; j++) {
            (*accurateAdder)[i][j] = accurate[i][j];
        }
    }
}

// int main() {
//     int ***outputs = NULL;
//     int **accurateAdder = NULL;
//     int rows, output_cols, numAdders;

//     // Initialize truth tables
//     initializeTruthTables(&outputs, &rows, &output_cols, &numAdders, &accurateAdder);

//     // Configuration to merge approximate adders
//     int configuration[] = {5}; // Merge adder 1, 2, and 3
//     int configSize = sizeof(configuration) / sizeof(configuration[0]);

//     // Dynamically create output headers for approximate adders
//     char output_headers[configSize + 1][10]; // configSize + 1 for Cout
//     for (int i = 0; i < configSize; i++) {
//         sprintf(output_headers[i], "Sum%d", i);
//     }
//     sprintf(output_headers[configSize], "Cout");

//     // Merge approximate truth tables
//     int **finalOutputs = NULL;
//     int finalRows, finalOutputCols;
//     mergeMultipleTruthTables(outputs, numAdders, rows, output_cols,
//                              configuration, configSize, &finalOutputs,
//                              &finalRows, &finalOutputCols);

//     // Display the merged approximate truth table
//     printf("Merged Approximate Truth Table:\n");
//     printTruthTable(finalOutputs, finalRows, finalOutputCols, output_headers);

//     // Merge accurate truth tables
//     int numBits = sizeof(configuration)/sizeof(int); // Number of bits for accurate adder
//     int **finalAccurateOutputs = NULL;
//     int finalAccurateRows, finalAccurateOutputCols;
//     mergeAccurateAdder(accurateAdder, rows, output_cols, numBits,
//                        &finalAccurateOutputs, &finalAccurateRows, &finalAccurateOutputCols);

//     // Dynamically create output headers for accurate adders
//     char accurate_headers[numBits + 1][10];
//     for (int i = 0; i < numBits; i++) {
//         sprintf(accurate_headers[i], "Sum%d", i);
//     }
//     sprintf(accurate_headers[numBits], "Cout");

//     // Display the merged accurate truth table
//     printf("Merged Accurate Truth Table:\n");
//     printTruthTable(finalAccurateOutputs, finalAccurateRows, finalAccurateOutputCols, accurate_headers);

//     // Free allocated memory
//     for (int k = 0; k < numAdders; k++) {
//         for (int i = 0; i < rows; i++) {
//             free(outputs[k][i]);
//         }
//         free(outputs[k]);
//     }
//     free(outputs);

//     for (int i = 0; i < rows; i++) {
//         free(accurateAdder[i]);
//     }
//     free(accurateAdder);

//     for (int i = 0; i < finalRows; i++) {
//         free(finalOutputs[i]);
//     }
//     free(finalOutputs);

//     for (int i = 0; i < finalAccurateRows; i++) {
//         free(finalAccurateOutputs[i]);
//     }
//     free(finalAccurateOutputs);

//     return 0;
// }
