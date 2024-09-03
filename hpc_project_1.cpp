#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>

#define SIZE 100 // Define the matrix size

using namespace std;

struct CRSMatrix {
    vector<int> values;
    vector<int> columns;
    vector<int> rowPointer;
    int nonZeros;
    CRSMatrix() : nonZeros(0) {
        rowPointer.resize(SIZE + 1, 0);
    }
};

// Function to generate a sparse matrix in CRS format
CRSMatrix generateSparseMatrixCRS(int size, double probability) {
    CRSMatrix crsMatrix;
    vector<int> valuesTemp, columnsTemp;
    crsMatrix.rowPointer[0] = 0;

    srand(time(0));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if ((rand() % 100) < (probability * 100)) {
                valuesTemp.push_back(rand() % 10 + 1);
                columnsTemp.push_back(j);
            }
        }
        crsMatrix.rowPointer[i + 1] = valuesTemp.size();
    }

    crsMatrix.nonZeros = valuesTemp.size();
    crsMatrix.values = valuesTemp;
    crsMatrix.columns = columnsTemp;



    return crsMatrix;
}

void generateRowCompressedMatrices(const CRSMatrix& X, CRSMatrix& B, CRSMatrix& C) {
    vector<int> valuesB, columnsB;
    vector<int> valuesC, columnsC;

    B.rowPointer.push_back(0);
    C.rowPointer.push_back(0);

    for (int i = 0; i < SIZE; i++) {
        int rowStart = X.rowPointer[i];
        int rowEnd = X.rowPointer[i + 1];
        int numNonZeros = rowEnd - rowStart;

        if (numNonZeros > 0) {
            for (int j = rowStart; j < rowEnd; j++) {
                valuesB.push_back(X.values[j]);
                columnsB.push_back(X.columns[j]);
                valuesC.push_back(X.values[j]); // Copy values for C as well
                columnsC.push_back(X.columns[j]); // Copy columns for C
            }
        } else {
            // Handle rows with no non-zero values if necessary
            // Here you should add your logic if needed
        }

        B.rowPointer.push_back(valuesB.size());
        C.rowPointer.push_back(valuesC.size());
    }

    B.nonZeros = valuesB.size();
    C.nonZeros = valuesC.size();
    B.values = valuesB;
    B.columns = columnsB; // Ensure B.columns is correctly populated
    C.values = valuesC;
    C.columns = columnsC;
}



void saveCRSToFile(const CRSMatrix& matrix, const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    file << "RowPointers: ";
    for (size_t i = 0; i < matrix.rowPointer.size(); ++i) {
        file << matrix.rowPointer[i] << " ";
    }
    file << std::endl;

    file << "Values: ";
    for (size_t i = 0; i < matrix.values.size(); ++i) {
        file << matrix.values[i] << " ";
    }
    file << std::endl;

    file << "Columns: ";
    for (size_t i = 0; i < matrix.columns.size(); ++i) {
        file << matrix.columns[i] << " ";
    }
    file << std::endl;

    file.close();
}


void printMatrixInfo(const CRSMatrix& matrix, const string& name) {
    cout << name << " nonZeros: " << matrix.nonZeros << endl;

    cout << name << " RowPointers: ";
    for (size_t i = 0; i < matrix.rowPointer.size(); ++i) {
        cout << matrix.rowPointer[i] << " ";
    }
    cout << endl;

    cout << name << " Values: ";
    for (size_t i = 0; i < matrix.values.size(); ++i) {
        cout << matrix.values[i] << " ";
    }
    cout << endl;

    cout << name << " Columns: ";
    for (size_t i = 0; i < matrix.columns.size(); ++i) {
        cout << matrix.columns[i] << " ";
    }
    cout << endl;
}

// Function to evaluate performance and save results
void evaluatePerformance(int size, const vector<double>& probabilities) {
    for (size_t i = 0; i < probabilities.size(); ++i) {
        double probability = probabilities[i];
        cout << "Generating matrices for set " << (i + 1) << " with probability " << probability << endl;

        CRSMatrix X = generateSparseMatrixCRS(size, probability);
        CRSMatrix Y = generateSparseMatrixCRS(size, probability);

        CRSMatrix XB, XC, YB, YC;
        generateRowCompressedMatrices(X, XB, XC);
        generateRowCompressedMatrices(Y, YB, YC);

        // printMatrixInfo(XB, "MatrixX_B");
        // printMatrixInfo(XC, "MatrixX_C");
        // printMatrixInfo(YB, "MatrixY_B");
        // printMatrixInfo(YC, "MatrixY_C");

        // Save matrices to files
        saveCRSToFile(XB, "MatrixX_B_" + to_string(i + 1) + ".txt");
        saveCRSToFile(XC, "MatrixX_C_" + to_string(i + 1) + ".txt");
        saveCRSToFile(YB, "MatrixY_B_" + to_string(i + 1) + ".txt");
        saveCRSToFile(YC, "MatrixY_C_" + to_string(i + 1) + ".txt");
    }
}

int main() {
    vector<double> probabilities = {0.01, 0.02, 0.05};
    evaluatePerformance(SIZE, probabilities);

    return 0;
}
