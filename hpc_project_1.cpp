#include <iostream>
#include <vector>
#include <omp.h>
#include <cstdlib>
#include <ctime>

#define SIZE 100  // Matrix size

using namespace std;

// Ordinary matrix multiplication
void matrixMultiplication(const vector<vector<int>>& X, const vector<vector<int>>& Y, vector<vector<int>>& Z) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            Z[i][j] = 0;
            for (int k = 0; k < SIZE; k++) {
                Z[i][j] += X[i][k] * Y[k][j];
            }
        }
    }
}

// Matrix multiplication using OpenMP
void matrixMultiplicationParallel(const vector<vector<int>>& X, const vector<vector<int>>& Y, vector<vector<int>>& Z, int num_threads) {
    #pragma omp parallel for collapse(2) num_threads(num_threads)
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            Z[i][j] = 0;
            for (int k = 0; k < SIZE; k++) {
                Z[i][j] += X[i][k] * Y[k][j];
            }
        }
    }
}

// Generate a sparse matrix
vector<vector<int>> generateSparseMatrix(double probability) {
    vector<vector<int>> matrix(SIZE, vector<int>(SIZE, 0));
    srand(time(0));
    
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if ((rand() % 100) < (probability * 100)) {
                matrix[i][j] = rand() % 10 + 1;  // Random integer between 1 and 10
            }
        }
    }
    return matrix;
}

// Row compression for sparse matrices
void rowCompression(const vector<vector<int>>& matrix, vector<vector<int>>& B, vector<vector<int>>& C) {
    for (int i = 0; i < SIZE; i++) {
        vector<int> b_row;
        vector<int> c_row;
        
        for (int j = 0; j < SIZE; j++) {
            if (matrix[i][j] != 0) {
                b_row.push_back(matrix[i][j]);
                c_row.push_back(j);
            }
        }
        
        if (b_row.empty()) {  // Handle empty row
            b_row.push_back(0);
            c_row.push_back(-1);  // -1 to indicate an empty row
        }
        
        B.push_back(b_row);
        C.push_back(c_row);
    }
}

// Evaluate performance of parallel matrix multiplication
void evaluatePerformance(const vector<vector<int>>& X, const vector<vector<int>>& Y, vector<vector<int>>& Z) {
    double best_time = 1e9;
    int best_threads = 1;

    for (int threads = 1; threads <= 64; threads++) {
        double start_time = omp_get_wtime();

        matrixMultiplicationParallel(X, Y, Z, threads);

        double end_time = omp_get_wtime();
        double elapsed_time = end_time - start_time;

        cout << "Threads: " << threads << ", Time: " << elapsed_time << " seconds" << endl;

        if (elapsed_time < best_time) {
            best_time = elapsed_time;
            best_threads = threads;
        }
    }

    cout << "Best performance with " << best_threads << " threads: " << best_time << " seconds" << endl;
}

int main() {
    srand(time(0));

    // Task 1: Ordinary Matrix Multiplication
    vector<vector<int>> X(SIZE, vector<int>(SIZE, rand() % 10 + 1));  // Example matrix X
    vector<vector<int>> Y(SIZE, vector<int>(SIZE, rand() % 10 + 1));  // Example matrix Y
    vector<vector<int>> Z(SIZE, vector<int>(SIZE, 0));  // Result matrix

    cout << "Ordinary Matrix Multiplication:" << endl;
    matrixMultiplication(X, Y, Z);
    
    // Task 2: Generate B and C matrices
    double probability = 0.05;  // Example sparsity probability
    vector<vector<int>> sparseX = generateSparseMatrix(probability);
    vector<vector<int>> XB, XC, YB, YC;
    
    rowCompression(sparseX, XB, XC);
    rowCompression(sparseX, YB, YC);  // Just an example, could use different matrices

    cout << "Generated B and C matrices for sparse matrix X." << endl;

    // Task 3: Evaluate performance with OpenMP
    cout << "Evaluating performance with varying threads:" << endl;
    evaluatePerformance(X, Y, Z);

    return 0;
}
