#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cstdlib>

using namespace std;
using namespace Eigen;

const double damping_factor = 0.95; // he so dan truyen
const double threshold = 1e-3; // nguong hoi tu
const int max_iterations = 1000; // so lan lap toi da

vector<double> pagerank(const SparseMatrix<double>& links_matrix) {
    int n = links_matrix.cols(); // so trang web
    VectorXd pr(n); // vector ban dau
    pr.setOnes(); // khoi tao gia tri 1 cho vector ban dau
    pr = pr / n; // chuan hoa vector ban dau

    SparseMatrix<double> A; // Khai bao ma tran A
    A = damping_factor * links_matrix.transpose(); // Gan gia tri cho ma tran A

    //Tao ma tran xac dinh muc do phan phoi xac suat chuyen tiep giua cac trang web
    MatrixXd B = (1 - damping_factor) / n * MatrixXd::Ones(n, n); 
    SparseMatrix<double> C = B.sparseView(); //Tao ma tran thua C tu ma tran B

    A += C; // Them gia tri (1 - damping_factor) / n vao cac phan tu cua ma tran A

    //Khoi tao solver de tinh lai vector pagerank moi
    // BiCGSTAB<SparseMatrix<double>> solver(A);
    // solver.setTolerance(threshold);

    for (int iter = 0; iter < max_iterations; iter++) {
        // tinh lai vector pagerank
        VectorXd pr_new = A * pr;
        // VectorXd pr_new = solver.solve(pr); 

        // tính khoang cach Euclidean giua vector pagerank ban dau va vector pagerank dc tinh toan
        double diff = (pr_new - pr).cwiseAbs().sum(); 
        if (diff < threshold) {
            break;  
        }
        // neu chua hoi tu thi gan lai vector pr va tiep tuc
        pr = pr_new;
    }

    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = pr(i);
    }
    return result;
}

int main() {
    int n = 5000; // Số trang web
    int m = 5000; // Số lien ket
    vector<Triplet<double>> triplets;

    //random data
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