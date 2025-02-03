#ifndef MERGING_H
#define MERGING_H

// Function Prototypes

// Function to merge two truth tables
void mergeTwoTruthTables(int **T1_outputs, int T1_rows, int T1_output_cols,
                         int **T2_outputs, int T2_rows, int T2_output_cols,
                         int ***T_outputs, int *T_rows, int *T_output_cols);

// Function to merge multiple truth tables based on configuration
void mergeMultipleTruthTables(int ***outputs, int numAdders, int rows, int output_cols,
                              int *configuration, int configSize, int ***finalOutputs,
                              int *finalRows, int *finalOutputCols);

// Function to merge accurate adders to match approximate adders' size
void mergeAccurateAdder(int **accurateAdder, int rows, int output_cols, int configSize,
                        int ***finalOutputs, int *finalRows, int *finalOutputCols);

// Function to print a truth table
void printTruthTable(int **outputs, int rows, int output_cols, char output_headers[][10]);

// Function to initialize truth tables for approximate and accurate adders
void initializeTruthTables(int ****outputs, int *rows, int *output_cols, int *numAdders,
                           int ***accurateAdder);

#endif // MERGING_H
