#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// Function prototypes
// void Validate_Inputs(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length);
// void Partial_PMF_Computation(int I, double A_p, double B_p, double C_in_p, double Sig, double *C_out_p, double *PMF_00, double *PMF_01, double *PMF_10, double *PMF_11);
// void Convolve(const double *x, int x_len, const double *y, int y_len, double *result);
void Probability_to_PMF(const int *Indexes, int Indexes_len, const double *IPM_full, const double *Err_full, double Sig, double *PMF);
void Convolve(const double *x, int x_len, const double *y, int y_len, double **result, int *result_len) ;
void Validate_Inputs(const int *Adder_Config, const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length);
double* Convolution(double* a, int len_a, double* b, int len_b) ;
void Partial_PMF_Computation( double *A_p, double *B_p, double C_in_p, double Sig,
    double *C_out_p, double *PMF_00, double *PMF_01, double *PMF_10, double *PMF_11,int rows,int numBits,int bitpos,int **finalOutputs, int *errorArray, double **CPM_00,double **CPM_10,double **CPM_01,double **CPM_11);

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
// Updated function to print a truth table with input bits
void printTruthTable(int **outputs, int **A_Array, int **B_Array, int * Cin_Array, int rows, int output_cols, int numBits) {
    // **1️⃣ Print Header Row**
    for (int i = 0; i < numBits; i++) {
        // printf("A%d\t", i);
    }
    for (int i = 0; i < numBits; i++) {
        // printf("B%d\t", i);
    }
    printf("Cin\t");  // **Carry input**
    
    for (int i = 0; i < numBits; i++) {
        // printf("Sum%d\t", i);
    }
    // printf("Cout\n");  // **Carry out**

     A_Array = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++) {
            A_Array[i] = (int *)malloc(numBits* sizeof(int));
        }
     B_Array = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++) {
            B_Array[i] = (int *)malloc(numBits* sizeof(int));
        }
     Cin_Array = (int *)malloc(rows * sizeof(int ));
 



    // **2️⃣ Print Truth Table Rows**
    
    for (int i = 0; i < rows; i++) { 
        // Print **A bits**
        for (int j = 0; j < numBits; j++) {
            A_Array[i][j] =  (i >> j) & 1;
           // printf("%d\t",A_Array[i][j] );  // Extract bit values for A
        }

        // Print **B bits**
        for (int j = 0; j < numBits; j++) {
             B_Array[i][j] =  (i >> (j + numBits)) & 1;
          //  printf("%d\t",B_Array[i][j]);  // Extract bit values for B
        }

          Cin_Array[i] =  (i >> (2 * numBits)) & 1;
        // Print **Cin bit** (Last bit before sum)

      // printf("%d\t", Cin_Array[i] );

        // Print **Sum and Cout from the truth table**
        for (int j = 0; j < output_cols; j++) {
          //  printf("%d\t", outputs[i][j]);
        }
        //printf("\n");
    }
}


// Hardcoded truth tables for approximate and accurate adders
void initializeTruthTables(int ****outputs, int *rows, int *output_cols, int *numAdders,
                           int ***accurateAdder) {
    *numAdders = 8;
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
    // int hardcoded_outputs[7][8][2] = {
        // {{0, 0}, {0, 0}, {0, 1}, {0, 1}, {1, 0}, {0, 1}, {0, 1}, {1, 1}},
        // {{1, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {0, 1}, {0, 1}},
        // {{1, 0}, {1, 0}, {0, 1}, {0, 1}, {1, 0}, {0, 1}, {0, 1}, {0, 1}},
        // {{0, 0}, {0, 1}, {0, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 0}, {1, 1}},
        // {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {0, 0}, {0, 1}, {1, 0}, {1, 1}},
        // {{0, 0}, {1, 0}, {1, 0}, {0, 0}, {1, 1}, {0, 1}, {0, 1}, {1, 1}},
        // {{0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {1, 1}, {1, 1}, {1, 1}}};
        int hardcoded_outputs[8][8][2] = {
            {{0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {1, 1}},
            {{0, 0}, {0, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {1, 1}},
            {{0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {0, 1}},
            {{0, 0}, {0, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {1, 1}},
            {{0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {0, 1}},
            {{0, 0}, {0, 0}, {1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}},
            {{0, 0}, {1, 0}, {1, 0}, {0, 1}, {1, 0}, {0, 1}, {1, 1}, {0, 1}},
            {{0, 0}, {0, 0}, {1, 0}, {1, 1}, {1, 0}, {1, 1}, {1, 1}, {1, 1}},};

    for (int k = 0; k < 8; k++) {
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


void processAndMerge(int ***outputs, int numAdders, int rows, int output_cols, 
                     int **accurateAdder, int *LV, int LV_size, int *configuration , const double *Probability_A_bits, const double *Probability_B_bits, double Probability_C_in, int length, double *MSE, double *MED,double **CPM_00,double **CPM_10,double **CPM_01,double **CPM_11) {
    printf("\nProcessing Truth Tables Based on Link Vector Groups:\n");

    int start = 0;
    int PMF_0_length = 7;  // Initial size
    int ** A_Array;
    int ** B_Array;
    int * Cin_Array;
    int K=0;
    int num_grps =0;
    int G =0;
    double C_out_p;
    double Min_val = -3;
    double sig;
    double final_length;
    int new_length =PMF_0_length;
    double* PMF_0 = (double*)calloc(PMF_0_length, sizeof(double));  // Adjust size as needed
    double* PMF_1 = (double*)calloc(PMF_0_length, sizeof(double));  // Adjust size as needed
    
    for(int i=0; i<LV_size;i++){
        if(LV[i]==-1){
            num_grps++;
        }
    }
        
    
    while (start < LV_size) {
        int groupSize = 0;

        while (start + groupSize < LV_size && LV[start + groupSize] != -1) {
            groupSize++;
        }
        
        int groupStartIndex = start;
        start += groupSize + 1;
        if (groupSize == 0) continue;
        printf("\ngrp start index: %d\n", groupStartIndex);
        printf("\nProcessing a Group of Size: %d\n", groupSize);

        int configurationSubset[groupSize];  
        int subsetIndex = 0;
        for (int i = 0; i < groupSize; i++) {  
            int bitIndex = LV[groupStartIndex + i];
            if (bitIndex == -1) continue;
            configurationSubset[subsetIndex++] = configuration[bitIndex];
        }

        if (subsetIndex == 0) {
            printf("Warning: No valid adders found for this group, skipping merge.\n");
            continue;
        }

        int **finalOutputs = NULL;
        int finalRows, finalOutputCols;
        mergeMultipleTruthTables(outputs, numAdders, rows, output_cols, 
                                 configurationSubset, subsetIndex, &finalOutputs, 
                                 &finalRows, &finalOutputCols);

        int **finalAccurateOutputs = NULL;
        int finalAccurateRows, finalAccurateOutputCols;
        mergeAccurateAdder(accurateAdder, rows, output_cols, subsetIndex, 
                           &finalAccurateOutputs, &finalAccurateRows, &finalAccurateOutputCols);

        printf("\nMerged Approximate Truth Table:\n");
        printTruthTable(finalOutputs,A_Array, B_Array , Cin_Array, finalRows, finalOutputCols, groupSize);

        printf("\nMerged Accurate Truth Table:\n");

        printTruthTable(finalAccurateOutputs,A_Array, B_Array , Cin_Array, finalAccurateRows, finalAccurateOutputCols, groupSize);

        // Compute and Store Error Array


        // 
        
        // printf("\nError Array:\n");
        // for (int i = 0; i < finalRows; i++) {
        //     int approx_decimal = 0, accurate_decimal = 0;
        //     int weight = 1;

        //     // Convert the group's binary values to decimal
        //     for (int j = 0; j < groupSize; j++) {
        //         approx_decimal += finalOutputs[i][j] * weight;
        //         accurate_decimal += finalAccurateOutputs[i][j] * weight;
        //         weight *= 2; // Increment weight for next bit
        //     }

        //     // Include carry-out in the decimal computation
        //     approx_decimal += finalOutputs[i][finalOutputCols - 1] * weight;
        //     accurate_decimal += finalAccurateOutputs[i][finalAccurateOutputCols - 1] * weight;

        //     // Compute the decimal error
        //     int error_decimal = approx_decimal - accurate_decimal;

        //     // Convert error back to binary
        //     for (int j = 0; j < groupSize; j++) {
        //         errorArray[i][j] = (error_decimal >> j) & 1; // Extract each bit of error
        //     }
        // }

        int *errorArray = (int *)malloc(finalRows * sizeof(int));

        printf("\nError Array:\n");
        for (int i = 0; i < finalRows; i++) {
            int approx_decimal = 0, accurate_decimal = 0;
            int weight = 1;

            // Convert the group's binary values to decimal
            for (int j = 0; j < groupSize; j++) {
                approx_decimal += finalOutputs[i][j] * weight;
                accurate_decimal += finalAccurateOutputs[i][j] * weight;
                weight *= 2; // Increment weight for next bit
            }

            // Include carry-out in the decimal computation
            approx_decimal += finalOutputs[i][finalOutputCols - 1] * weight;
            accurate_decimal += finalAccurateOutputs[i][finalAccurateOutputCols - 1] * weight;

            // Compute the decimal error
            int error_decimal = approx_decimal - accurate_decimal;
            errorArray[i] = error_decimal;
            // printf("error of  %d\n",errorArray[i]);
            // Convert error back to binary
            // for (int j = 0; j < groupSize; j++) {
            //     errorArray[i][j] = (error_decimal >> j) & 1; // Extract each bit of error
            // }
        }



        // printf("\nError Array:\n");
        // for (int i = 0; i < finalRows; i++) {
        //   //  printf("Row %d: ", i);
        //     for (int j = 0; j < groupSize; j++) {
        //         errorArray[i][j] = (finalOutputs[i][j] - finalAccurateOutputs[i][j]) + 
        //                            2 * (finalOutputs[i][finalOutputCols - 1] - finalAccurateOutputs[i][finalAccurateOutputCols - 1]);
        //        // printf("%d\t", errorArray[i][j]);
        //     }
        //    // printf("\n");
        // }

///////////////////////////////////////////////////////////////////
         // Initialize variables
    

    // Loop through the bits in the adder configuration

        sig = pow(2, K);
         int PMF_length = (int)(6 * sig + 1);
        double* PMF_00 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
        double* PMF_01 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
        double* PMF_10 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
         double* PMF_11 = (double*)calloc(PMF_length, sizeof(double)); // Adjust size as needed
        
         double *A_p = (double *)malloc(groupSize * sizeof(double ));
         double *B_p = (double *)malloc(groupSize * sizeof(double ));
        
         for (int i = 0; i < groupSize; i++) {
         A_p[i] = Probability_A_bits[i+G];
         B_p[i]= Probability_B_bits[i+G];
        
        }

        double C_in_p;
        if (K == 0) {
            C_in_p = Probability_C_in;
                           // printf("iiiiiiiiiiiiiiiiiiii %f\n",C_in_p);

                                  
        } else {
            C_in_p = C_out_p;
            //printf("iiiiiiiiiiiiiiiiiiii %f\n",C_in_p);
        }
       C_out_p = 0.0;

        // Compute the partial PMF values
        Partial_PMF_Computation(A_p, B_p, C_in_p, sig, &C_out_p, PMF_00, PMF_01, PMF_10, PMF_11,finalRows,groupSize,0,finalOutputs,errorArray,CPM_00,CPM_10,CPM_01,CPM_11);

        
        // Update PMFs for current bit
        if (K == 0) {
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
            // for(int i = 0; i < PMF_length; i++){
            //         printf("iiiiiiiiiiiiiiiiiiii %f\n",PMF_0[i]);
            // }

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
    // for (int i = 0; i < new_length; i++) {
    //         printf("PMF_0= %.10f, PMF_1=%.10f\n",PMF_0[i],PMF_1[i]);
    // }
     K++;  

     if (K != num_grps ) {
             printf("hhhhhhhhhhhhhhhhhaaaaaaaaaaaaaaaaaaa  %d\n",K);
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
        
   // }

/////////////////////////////////////////////////////////////////////////////////////////

        for (int i = 0; i < finalRows; i++) free(finalOutputs[i]);
        free(finalOutputs);
        for (int i = 0; i < finalAccurateRows; i++) free(finalAccurateOutputs[i]);
        free(finalAccurateOutputs);
        free(errorArray);

    G = G + groupSize;
    }
    
    // Combine the final PMFs
    double* PMF = (double*)calloc(final_length, sizeof(double)); // Adjust size as needed
    for (int i = 0; i < final_length; i++) {  // Adjust size as needed
        PMF[i] = PMF_0[i] + PMF_1[i];
        // printf("iiiiiiiiiiiiiiiiiiii %f\n",PMF_1[i]);

    }
     printf("data pointa %f\n",final_length);
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


//  #define NUM_BITS 5
// #define NUM_COMBINATIONS (1 << NUM_BITS) // 2^5 = 32

// Function to compute conditional probability P(a_i | a_j)
double compute_conditional_probability(int a_i, int a_j, double **CPM_11, double **CPM_01, double **CPM_10, double** CPM_00,int i,int j) {
    if (a_i == 1 && a_j == 1) {
        return CPM_11[i][j];
    } else if (a_i == 0 && a_j == 1) {
        return CPM_01[i][j];
    } else if (a_i == 1 && a_j == 0) {
        return CPM_10[i][j];
    } else {
        return CPM_00[i][j];
    }
}

// Function to compute joint probability P(a0, a1, ..., a4)
double compute_joint_probability(int * bits, double *bit_probabilities, double **CPM_11, double **CPM_01, double **CPM_10, double** CPM_00,int numBits) {
    double joint_prob = 1.0;
    double total_prob = 0;
    for (int i = 0; i < numBits; i++) {
        if (i == 0) {
            // P(a0) is given directly
            joint_prob *= (bits[i] == 1) ? bit_probabilities[i] : (1 - bit_probabilities[i]);
        } else {
            // P(ai | a0, a1, ..., a_{i-1}) is computed using the CPMs
            int a_i = bits[i];
            int a_j = bits[i - 1];
            double conditional_prob = compute_conditional_probability(a_i, a_j, CPM_11, CPM_01, CPM_10, CPM_00, i,i-1);
            // printf("The joint prob: %d %d\n",i,i-1);
            joint_prob *= conditional_prob;
            
            // if(i>1){
            //     conditional_prob = compute_conditional_probability(a_i, a_j, CPM_11, CPM_01, CPM_10, CPM_00, i,i-2);
            //     joint_prob /= conditional_prob;
            // }
        }
    }
    // for(int i = 0; i < numBits; i++) {
    //     total_prob = total_prob + bit_probabilities[i];
    // }
    // if (numBits >2){
    //     joint_prob = total_prob - joint_prob; 
    // }
    return joint_prob;
}

// Function to compute the IPM
void compute_ipm(double *IPM, double *bit_probabilities, double **CPM_11, double **CPM_01, double **CPM_10, double **CPM_00,int **A,int numBits, int rows) {
    // Generate all possible combinations of bits
    for (int i = 0; i < rows; i++) {
        int *bits= (int *)malloc(numBits * sizeof(int ));
        for (int j = 0; j < numBits; j++) {
           bits[j] = A[i][j];
        }
        // Compute the joint probability for the current combination
        IPM[i] = compute_joint_probability(bits, bit_probabilities, CPM_11, CPM_01, CPM_10, CPM_00,numBits);
    }

    // // Normalize the IPM to ensure the probabilities sum to 1
    // double total_probability = 0.0;
    // for (int i = 0; i < rows; i++) {
    //     total_probability += IPM[i];
    // }
    // for (int i = 0; i < rows; i++) {
    //     IPM[i] /= total_probability;
    // }
}




void Partial_PMF_Computation( double *A_p, double *B_p, double C_in_p, double Sig,
                             double *C_out_p, double *PMF_00, double *PMF_01, double *PMF_10, double *PMF_11,int rows,int numBits,int bitpos,int **finalOutputs, int *errorArray, double **CPM_00,double **CPM_10,double **CPM_01,double **CPM_11) {
    // Define all possible combinations of A, B, and C_in

    int ** A= (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++) {
            A[i] = (int *)malloc(numBits* sizeof(int));
        }
    int ** B = (int **)malloc(rows * sizeof(int *));
        for (int i = 0; i < rows; i++) {
            B[i] = (int *)malloc(numBits* sizeof(int));
        }
     int *C_in = (int *)malloc(rows * sizeof(int ));
 
    // **2️⃣ Print Truth Table Rows**
    for (int i = 0; i < rows; i++) { 
        // Print **A bits**
        for (int j = 0; j < numBits; j++) {
            A[i][j] =  (i >> j) & 1;
            // printf("%d\t",A[i][j] );  // Extract bit values for A
        }

        // Print **B bits**
        for (int j = 0; j < numBits; j++) {
             B[i][j] =  (i >> (j + numBits)) & 1;
             //printf("%d\t",B[i][j]);  // Extract bit values for B
        }

          C_in[i] =  (i >> (2 * numBits)) & 1;
        // Print **Cin bit** (Last bit before sum)
                //   printf(" %d\t", C_in[i] );

                //   printf("\n" );  // Extract bit values for A

    
    }

    // int A[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    // int B[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    // int C_in[8] = {0, 1, 0, 1, 0, 1, 0, 1};
    double *A_v_p = (double *)malloc(rows * sizeof(double ));
    double *B_v_p = (double *)malloc(rows * sizeof(double ));

    double C_in_v_p[rows], IPM[rows];
    double C_out[rows];
    double Err[rows];
    

    
        
    // if(numBits ==1){
    // // Compute individual probabilities
    //     for (int i = 0; i < rows; i++) {
    //         A_v_p[i] = (A[i][bitpos] * A_p[0]) + ((!A[i][bitpos]) * (1 - A_p[0]));  // Bitwise NOT operation on A
    //         B_v_p[i] = (B[i][bitpos] * B_p[0]) + ((!B[i][bitpos]) * (1 - B_p[0]));  // Bitwise NOT operation on B
    //         C_in_v_p[i] = (C_in[i] * C_in_p) + ((!C_in[i]) * (1 - C_in_p)); 
    //         IPM[i] = A_v_p[i] * B_v_p[i] * C_in_v_p[i];

    //     }   
    // }   
    // else{
        compute_ipm(A_v_p, A_p, CPM_11, CPM_01, CPM_10, CPM_00,A,numBits,rows);             // Bitwise NOT operation on A
        compute_ipm(B_v_p, B_p, CPM_11, CPM_01, CPM_10, CPM_00,B,numBits,rows); // Bitwise NOT operation on B
       
         for (int i = 0; i < rows; i++) {
            C_in_v_p[i] = (C_in[i] * C_in_p) + ((!C_in[i]) * (1 - C_in_p)); 
            IPM[i] = A_v_p[i] * B_v_p[i] * C_in_v_p[i];
        } 
    // }

       // Output: IPM array


    for (int i = 0; i < rows; i++) { 
        C_out[i]=finalOutputs[i][numBits] ;
         Err[i]=errorArray[i];
                

    }

    // Compute carry-out probability
    *C_out_p = 0;
    for (int i = 0; i < rows; i++) {
        *C_out_p += C_out[i] * IPM[i];
                
    }


    
    for(int i = 0; i <numBits; i++){
        printf("B_p[%d] = %.10f\n",
       i, B_p[i]);
    }

    printf("Before Partial_PMF_Computation: K=%d, C_in_p=%.10f\n", bitpos, C_in_p);
    for(int i = 0; i <rows; i++){
        printf("IPM[%d] = %.10f, A_v_p=%.10f, B_v_p=%.10f, C_in_v_p=%.10f,C_out=%.10f\n",
       i, IPM[i], A_v_p[i], B_v_p[i], C_in_v_p[i],C_out[i]);
    }

    printf("C_out_p after update: %.10f\n", *C_out_p);


    // Compute partial error PMFs
    int Indexes[rows];
    int count;

    // PMF_00: C_in = 0 and C_out = 0
    count = 0;
    for (int i = 0; i < rows; i++) {
        if (C_in[i] == 0 && C_out[i] == 0) Indexes[count++] = i;

    }
    Probability_to_PMF(Indexes, count, IPM, Err, Sig, PMF_00);

    // PMF_01: C_in = 0 and C_out = 1
    count = 0;
    for (int i = 0; i < rows; i++) {
        if (C_in[i] == 0 && C_out[i] == 1) Indexes[count++] = i;
    }
    Probability_to_PMF(Indexes, count, IPM, Err, Sig, PMF_01);

    // PMF_10: C_in = 1 and C_out = 0
    count = 0;
    for (int i = 0; i < rows; i++) {
        if (C_in[i] == 1 && C_out[i] == 0) Indexes[count++] = i;
    }
    Probability_to_PMF(Indexes, count, IPM, Err, Sig, PMF_10);

    // PMF_11: C_in = 1 and C_out = 1
    count = 0;
    for (int i = 0; i < rows; i++) {
        if (C_in[i] == 1 && C_out[i] == 1) {
            Indexes[count++] = i;

        }
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





