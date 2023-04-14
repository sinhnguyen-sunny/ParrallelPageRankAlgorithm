#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

const double damping_factor = 0.85; 
const double threshold = 1e-3;
const int max_iterations = 1000;

vector<double> pagerank(const SparseMatrix<double>& links_matrix) {
    int n = links_matrix.cols();
    VectorXd pr(n);
    pr.setOnes();
    pr = pr / n;

    SparseMatrix<double> A;
    A = damping_factor * links_matrix.transpose();

    MatrixXd B = (1 - damping_factor) / n * MatrixXd::Ones(n, n);
    SparseMatrix<double> C = B.sparseView();

    A += C;

    BiCGSTAB<SparseMatrix<double>> solver(A);
    solver.setTolerance(threshold);

    for (int iter = 0; iter < max_iterations; iter++) {
        VectorXd pr_new = solver.solve(pr);
        double diff = (pr_new - pr).cwiseAbs().sum();
        if (diff < threshold) {
            break;
        }
        pr = pr_new;
    }

    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = pr(i);
    }
    return result;
}

int main() {
    int n = 700;//web
    int m = 1000000;// Số liên kết
    vector<Triplet<double>> triplets;

    for (int i = 0; i < m; i++) {
        int from = rand() % n;
        int to = rand() % n;
        triplets.emplace_back(from, to, 1.0);
    }

    SparseMatrix<double> links_matrix(n, n);

    links_matrix.setFromTriplets(triplets.begin(), triplets.end());

    vector<double> pr = pagerank(links_matrix);

    for (int i = 0; i < n; i++) {
        cout << "Page " << i << ": " << pr[i] << endl;
    }

    return 0;
}