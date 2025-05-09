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
void printTruthTable(int **outputs, int **A_Array, int **B_Array, int * Cin_Array, int rows, int output_cols, int numBits) ;

// Function to print a truth table
// void printTruthTable(int **outputs, int rows, int output_cols, char output_headers[][10]);

// Function to initialize truth tables for approximate and accurate adders
void initializeTruthTables(int ****outputs, int *rows, int *output_cols, int *numAdders,
                           int ***accurateAdder);


// void processAndMerge(int ***outputs, int numAdders, int rows, int output_cols, 
//                      int **accurateAdder, int *LV, int LV_size, int *configuration);


// void processAndMerge(int ***outputs, int numAdders, int rows, int output_cols, 
//                      int **accurateAdder, int *LV, int LV_size, int *configurationconst , const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length, double *MSE, double *MED, double** CPM_00, double** CPM_01, double** CPM_10, double** CPM_11) ;
void processAndMerge(int ***outputs, int numAdders, int rows, int output_cols, 
    int **accurateAdder, int *LV, int LV_size, int *configuration , const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length, double *MSE, double *MED,double **CPM_00,double **CPM_10,double **CPM_01,double **CPM_11,double *output_pmf_accum ,double* PMF );

#endif // MERGING_H
