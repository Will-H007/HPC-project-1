#include <iostream>
#include <vector>
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
            b_row.push_back(0);
            c_row.push_back(0);
            c_row.push_back(0);
        }
        
        B.push_back(b_row);
        C.push_back(c_row);
    }
}

int main() {
    double probabilities[] = {0.01, 0.02, 0.05};

    for (double prob : probabilities) {
        // Generate sparse matrices X and Y
        vector<vector<int>> X = generateSparseMatrix(prob);
        vector<vector<int>> Y = generateSparseMatrix(prob);

        // Compress matrices X and Y
        vector<vector<int>> XB, XC, YB, YC;
        rowCompression(X, XB, XC);
        rowCompression(Y, YB, YC);

        // Print some details for verification
        cout << "Probability: " << prob << endl;
        cout << "XB size: " << XB.size() << ", XC size: " << XC.size() << endl;
        cout << "YB size: " << YB.size() << ", YC size: " << YC.size() << endl;

        // Optional: Perform matrix multiplication to verify
        vector<vector<int>> Z(SIZE, vector<int>(SIZE, 0));
        matrixMultiplication(X, Y, Z);

        cout << "Compressed matrices generated successfully for probability " << prob << endl;
    }

    return 0;
}
