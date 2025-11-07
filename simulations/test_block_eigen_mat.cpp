#include "include/utils.h"

using namespace Eigen;
int main(){
    
    int m = 8;
    int n = 2;
    std::vector<SparseMatrix<double>> blocks;
    blocks.reserve(m);
    
    for (int k = 0; k < m; ++k) {
        SparseMatrix<double> A(n, n);
        std::vector<Triplet<double>> trip;
        // Fill with "simple" numbers depending on k
        trip.emplace_back(0, 0, 10 * (k + 1) + 1); // e.g. 11, 21, 31
        trip.emplace_back(0, 1, 10 * (k + 1) + 2);
        trip.emplace_back(1, 0, 10 * (k + 1) + 3);
        trip.emplace_back(1, 1, 10 * (k + 1) + 4);
        A.setFromTriplets(trip.begin(), trip.end());
        blocks.push_back(A);
    }
    SparseMatrix<double> M = blockDiag(blocks); 

    std::cout << M << std::endl;

    M.leftCols(n) = 0.5*M.leftCols(n);
    M.rightCols(n) = 0.5*M.rightCols(n);

    std::cout << M << std::endl;

}

