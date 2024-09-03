#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip> // for std::setprecision
#include <cstdlib> // for rand() and srand()
#include <ctime>   // for time()
#include <omp.h>

using namespace std;

// Function to generate a random sparse matrix
void generateSparseMatrix(int size, double probability, vector<vector<int>>& matrix) {
    matrix.resize(size, vector<int>(size, 0));
    srand(static_cast<unsigned>(time(0))); // Seed for random number generation

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (static_cast<double>(rand()) / RAND_MAX < probability) {
                matrix[i][j] = rand() % 10 + 1; // Random value between 1 and 10
            }
        }
    }
}

// Function to break down a matrix into row compression format
void breakdownMatrix(const vector<vector<int>>& A, vector<vector<int>>& B, vector<vector<int>>& C) {
    int size = A.size();

    B.resize(size);
    C.resize(size);

    for (int i = 0; i < size; ++i) {
        vector<int> values;
        vector<int> columns;

        for (int j = 0; j < size; ++j) {
            if (A[i][j] != 0) {
                values.push_back(A[i][j]);
                columns.push_back(j);
            }
        }

        if (values.empty()) {
            // If no non-zero elements, store two consecutive 0s
            B[i] = {0, 0};
            C[i] = {0, 0};
        } else {
            // Store non-zero elements and their column indices
            B[i] = values;
            C[i] = columns;
        }
    }
}

// Function to write a matrix to a file
void writeMatrixToFile(const vector<vector<int>>& matrix, const string& filename) {
    ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    for (const auto& row : matrix) {
        for (int val : row) {
            file << val << " ";
        }
        file << endl;
    }

    file.close();
}
 

vector<vector<int>> multiplySparseMatrices(
    const vector<vector<int>>& B1, const vector<vector<int>>& C1,
    const vector<vector<int>>& B2, const vector<vector<int>>& C2, int size) {

    vector<vector<int>> result(size, vector<int>(size, 0));

    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            int sum = 0;
            for (int k = 0; k < B1[i].size(); ++k) {
                int value1 = B1[i][k];
                int col1 = C1[i][k];

                // Efficient lookup with std::find can be slow; consider using unordered_map for large matrices
                auto it = find(C2[col1].begin(), C2[col1].end(), j);
                if (it != C2[col1].end()) {
                    int index2 = it - C2[col1].begin();
                    int value2 = B2[col1][index2];
                    sum += value1 * value2;
                }
            }
            result[i][j] = sum;
        }
    }

    return result;
}





// Function to multiply two dense matrices
vector<vector<int>> multiplyDenseMatrices(const vector<vector<int>>& A, const vector<vector<int>>& B) {
    int size = A.size();
    vector<vector<int>> result(size, vector<int>(size, 0));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            int sum = 0;
            for (int k = 0; k < size; ++k) {
                sum += A[i][k] * B[k][j];
            }
            result[i][j] = sum;
        }
    }

    return result;
}

// Function to check the results of sparse matrix multiplication against dense multiplication
void checkResults(const vector<vector<int>>& sparseResult, const vector<vector<int>>& denseResult) {
    int size = sparseResult.size();
    bool match = true;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (sparseResult[i][j] != denseResult[i][j]) {
                cout << "Mismatch at (" << i << ", " << j << "): "
                     << "Sparse result = " << sparseResult[i][j] << ", "
                     << "Dense result = " << denseResult[i][j] << endl;
                match = false;
            }
        }
    }

    if (match) {
        cout << "Results match!" << endl;
    }
}

int main() {
    int size = 100; // Adjust size as needed
    double probabilities[] = {0.01, 0.02, 0.05}; // Probabilities to test

    for (double probability : probabilities) {
        vector<vector<int>> X, Y;
        vector<vector<int>> XB, XC, YB, YC;

        // Generate random sparse matrices X and Y
        generateSparseMatrix(size, probability, X);
        generateSparseMatrix(size, probability, Y);

        // Break down matrices X and Y into matrices XB, XC, YB, and YC
        breakdownMatrix(X, XB, XC);
        breakdownMatrix(Y, YB, YC);

        // Write matrices XB, XC, YB, and YC to files
        string probStr = to_string(probability).substr(2); // Get string after "0."
        writeMatrixToFile(XB, "XB_" + probStr + ".txt");
        writeMatrixToFile(XC, "XC_" + probStr + ".txt");
        writeMatrixToFile(YB, "YB_" + probStr + ".txt");
        writeMatrixToFile(YC, "YC_" + probStr + ".txt");

        cout << "Matrices for probability " << probability << " written to files." << endl;

        // Time the sparse matrix multiplication
        auto start = chrono::high_resolution_clock::now();

        vector<vector<int>> sparse_result = multiplySparseMatrices(XB, XC, YB, YC, size);

        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "Sparse matrix multiplication took " << fixed << setprecision(6) << elapsed.count() << " seconds." << endl;


        // Perform dense matrix multiplications for comparison
        vector<vector<int>> dense_result = multiplyDenseMatrices(X, Y);


        // Check results
        checkResults(sparse_result, dense_result);
    }

    return 0;
}
