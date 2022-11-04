#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

#define DELIMETER '\t'

vector<vector<double>> read_matrix_from_tsv(const string &file_name) {
    string line;
    ifstream myfile("matrix.tsv");

    vector<vector<double>> matrix;

    vector<double> row;
    if (myfile.is_open()) {
        while (getline(myfile, line)) {

            size_t pos = 0;
            string word;
            while ((pos = line.find(DELIMETER)) != std::string::npos) {
                word = line.substr(0, pos);
                row.push_back(stod(word));
                line.erase(0, pos + 1);
            }
            row.push_back(stod(line));

            matrix.push_back(row);
            row.clear();
        }
        myfile.close();
    } else cout << "Unable to open file";

    return matrix;
}

void print_matrix(const vector<vector<double>> &matrix) {
    for (const auto &row: matrix) {
        for (auto element: row) cout << element << " ";
        cout << endl;
    }
}

void solve_gauss(vector<vector<double>> matrix) {
    vector<double> res(matrix.size());

    for (int i = 0; i < matrix.size() - 1; i++) {
        for (int j = i + 1; j < matrix.size(); j++) {
            double coeff_here = matrix[j][i] / matrix[i][i];
            for (int k = 0; k < matrix[0].size(); k++) {
                matrix[j][k] = matrix[j][k] - coeff_here * matrix[i][k];
            }
        }
    }

    for (int i = matrix.size() - 1; i >= 0; --i) {
        res[i] = matrix[i][matrix.size()];

        for (int j = i + 1; j < matrix.size(); j++) {
            if (i != j) {
                res[i] = res[i] - matrix[i][j] * res[j];
            }
        }
        res[i] = res[i] / matrix[i][i];


    }
    cout << "RESULT:" << endl;
    for (int i = 0; i < matrix.size(); ++i) cout << "x" << i << "=" << res[i] << "; ";
    cout << endl;
}

vector<vector<double>> minor(const vector<vector<double>> &matrix, int i0, int j0) {
    vector<vector<double>> result(matrix.size() - 1, vector<double>(matrix.size() - 1));
    for (int i = 0; i < matrix.size(); ++i) {
        if (i != i0) {
            for (int j = 0; j < matrix[0].size(); ++j) {
                if (j != j0) result[i].push_back(matrix[i][j]);
            }
        }
    }

    return result;
}

double det(const vector<vector<double>> &matrix) {
    if (matrix.size() == 1) return matrix[0][0];
    if (matrix.size() == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double result = 0;
    for (int j = 0; j < matrix[0].size(); ++j) {
        result += pow(-1, j) * matrix[0][j] * det(minor(matrix, 0, j));
    }
    return result;
}

void solve_cramer(vector<vector<double>> matrix) {
    vector<vector<double>> main_matrix(matrix.size(), vector<double>(matrix.size()));
    for (int i = 0; i < matrix.size(); ++i)
        for (int j = 0; j < matrix.size(); ++j)
            main_matrix[i][j] = matrix[i][j];

    vector<double> constant_terms(matrix.size());
    for (int i = 0; i < matrix.size(); ++i) constant_terms[i] = matrix[i][matrix.size()];

    double main_det = det(main_matrix);
    vector<double> res(matrix.size());
    for (int j = 0; j < matrix.size(); ++j) {
        vector<vector<double>> j_matrix = main_matrix;
        for (int i = 0; i < matrix.size(); ++i) {
            j_matrix[i][j] = constant_terms[j];
        }
        res[j] = det(j_matrix) / main_det;
    }


    cout << "RESULT:" << endl;
    for (int i = 0; i < matrix.size(); ++i) cout << "x" << i << "=" << res[i] << "; ";
    cout << endl;
}

int main() {
    auto matrix = read_matrix_from_tsv("matrix.tsv");

    print_matrix(matrix);

    solve_gauss(matrix);

    solve_cramer(matrix);

    return 0;
}
