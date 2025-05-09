#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_BITS 32
#define N 8
#define SIZE 256
#define MAX_NODES 1024

typedef struct {
    int bit_position;               // Output bit index this adder contributes to
    double error_pmf[SIZE];         // Error distribution for this adder
} AdderNode;


AdderNode adder_nodes[MAX_NODES];
int adder_node_count = 0;

void compute_bit_probabilities(const double* P, double* P_bits);
void compute_partial_product_row_distribution(int a_bit_idx, const double* P_A_bits, const double* P_B, double* PP_row_dist);
void convolve_distributions(const double* dist1, const double* dist2, double* result);
void dreamx_estimate_error(const double* A_dist, const double* B_dist, double* output_dist, double* error_pmf);
double compute_expected_error(const double* error_pmf);
double compute_expected_value(const double* dist);
double femto_compute_total_weighted_error();
void print_distribution(const char* label, const double* dist);
void femto_convolve_all_error_pmfs(double* final_error_pmf);
int main() {
    // Uniform input distributions
    double P_A[SIZE], P_B[SIZE];
    for (int i = 0; i < SIZE; ++i) {
        P_A[i] = 1.0 / SIZE;
        P_B[i] = 1.0 / SIZE;
    }

    // Bit-level probabilities
    double P_A_bits[N], P_B_bits[N];
    compute_bit_probabilities(P_A, P_A_bits);
    compute_bit_probabilities(P_B, P_B_bits);

    // Partial product rows (from AND gates)
    double PP_rows[N][SIZE] = {0};
    for (int i = 0; i < N; ++i) {
        compute_partial_product_row_distribution(i, P_A_bits, P_B, PP_rows[i]);
    }
    double total_error = 0.0;
    double current_output[SIZE] = {0};
    for (int i = 0; i < SIZE; ++i)
        current_output[i] = PP_rows[0][i];

        for (int i = 0; i < SIZE; ++i) current_output[i] = PP_rows[0][i];

        for (int i = 1; i < N; ++i) {
            double shifted_pp[SIZE] = {0};
            double new_output[SIZE] = {0};
            double error_pmf[SIZE] = {0};
        
            // Shift PP_rows[i] left by i bits
            for (int j = 0; j < SIZE; ++j) {
                if (j << i < SIZE)
                    shifted_pp[j << i] = PP_rows[i][j];
            }
        
            // DREAMx adds current_output + shifted partial product
            dreamx_estimate_error(current_output, shifted_pp, new_output, error_pmf);
        
            // Register adder node for FEMTO traversal
            adder_nodes[adder_node_count].bit_position = i;
            for (int j = 0; j < SIZE; ++j)
                adder_nodes[adder_node_count].error_pmf[j] = error_pmf[j];
            adder_node_count++;
        
            // Accumulate expected error
            double err = compute_expected_error(error_pmf);
            total_error += err;
        
            // Update current output
            for (int j = 0; j < SIZE; ++j)
                current_output[j] = new_output[j];
        }
        
        

    // Output results
    print_distribution("Final Output Distribution", current_output);
    printf("\nExpected Value: %.6f\n", compute_expected_value(current_output));
    printf("FEMTO Estimated Weighted Error: %.6f\n", femto_compute_total_weighted_error());

    double femto_error_pmf[SIZE] = {0};
    femto_convolve_all_error_pmfs(femto_error_pmf);
    print_distribution("FEMTO Final Error PMF", femto_error_pmf);


    return 0;
}

void compute_bit_probabilities(const double* P, double* P_bits) {
    for (int bit = 0; bit < N; ++bit) {
        P_bits[bit] = 0.0;
        for (int i = 0; i < SIZE; ++i) {
            if ((i >> bit) & 1)
                P_bits[bit] += P[i];
        }
    }
}

void compute_partial_product_row_distribution(int a_bit_idx, const double* P_A_bits, const double* P_B, double* PP_row_dist) {
    double prob = P_A_bits[a_bit_idx];
    for (int i = 0; i < SIZE; ++i) {
        int shifted = i << a_bit_idx;
        if (shifted < SIZE)
            PP_row_dist[shifted] = P_B[i] * prob;
    }
}

void convolve_distributions(const double* dist1, const double* dist2, double* result) {
    for (int i = 0; i < SIZE; ++i)
        result[i] = 0;
    for (int i = 0; i < SIZE; ++i)
        for (int j = 0; j < SIZE; ++j)
            if (i + j < SIZE)
                result[i + j] += dist1[i] * dist2[j];
}

// 💡 You will implement your DREAMx logic here
void dreamx_estimate_error(const double* A_dist, const double* B_dist, double* output_dist, double* error_pmf) {
    convolve_distributions(A_dist, B_dist, output_dist);

    // For now: dummy error model (change this!)
    for (int i = 0; i < SIZE; ++i)
        error_pmf[i] = 0.0;
}

// FEMTO weighted error = ∑ error_pmf[i] × i × 2^bit_position
double femto_compute_total_weighted_error() {
    double total = 0.0;
    for (int i = 0; i < adder_node_count; ++i) {
        int pos = adder_nodes[i].bit_position;
        for (int j = 0; j < SIZE; ++j) {
            total += fabs(j) * adder_nodes[i].error_pmf[j] * (1 << pos);
        }
    }
    return total;
}

void femto_convolve_all_error_pmfs(double* final_error_pmf) {
    // Initialize final_error_pmf as delta distribution (identity for convolution)
    for (int i = 0; i < SIZE; ++i)
        final_error_pmf[i] = (i == 0) ? 1.0 : 0.0;

    double temp[SIZE];

    for (int i = 0; i < adder_node_count; ++i) {
        // Scale error PMF by 2^bit_position
        double scaled_pmf[SIZE] = {0};

        for (int j = 0; j < SIZE; ++j) {
            int shifted_idx = j << adder_nodes[i].bit_position;
            if (shifted_idx < SIZE)
                scaled_pmf[shifted_idx] = adder_nodes[i].error_pmf[j];
        }

        // Convolve with accumulated result
        convolve_distributions(final_error_pmf, scaled_pmf, temp);

        // Update final_error_pmf with result
        for (int k = 0; k < SIZE; ++k)
            final_error_pmf[k] = temp[k];
    }
}

double compute_expected_error(const double* error_pmf) {
    double error = 0.0;
    for (int i = 0; i < SIZE; ++i)
        error += fabs(i) * error_pmf[i];
    return error;
}

double compute_expected_value(const double* dist) {
    double expected = 0.0;
    for (int i = 0; i < SIZE; ++i)
        expected += i * dist[i];
    return expected;
}

void print_distribution(const char* label, const double* dist) {
    printf("\n%s:\n", label);
    for (int i = 0; i < SIZE; ++i) {
        if (dist[i] > 1e-6)
            printf("  %3d: %.6f\n", i, dist[i]);
    }
}
