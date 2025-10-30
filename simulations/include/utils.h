#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <chrono> 
#include <filesystem>
#include <limits>
#include <algorithm>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
read_mtx(const std::string& file_name) {
    Eigen::SparseMatrix<T> buff;
    Eigen::loadMarket(buff, file_name);
    return buff;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> 
read_TXT(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<std::vector<T>> data;
    std::string line;

    // Read the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<T> row;
        T value;

        // Extract numbers from each line
        while (ss >> value) {
            row.push_back(value);
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    file.close();

    if (data.empty()) {
        throw std::runtime_error("The file is empty or has invalid format.");
    }

    // Determine matrix dimensions
    size_t rows = data.size();
    size_t cols = data[0].size();

    // Ensure all rows have the same number of columns
    for (const auto& row : data) {
        if (row.size() != cols) {
            throw std::runtime_error("Inconsistent number of columns in the file.");
        }
    }

    // Create Eigen::MatrixXd and populate it
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix(rows, cols);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matrix(i, j) = data[i][j];
        }
    }

    return matrix;
}

template<typename T> void eigen2ext(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M, const std::string& sep, const std::string& filename, bool append = false){
    std::ofstream file;

    if(!append) 
        file.open(filename);
    else
        file.open(filename, std::ios_base::app); 
    
    for(std::size_t i = 0; i < M.rows(); ++i) {
            for(std::size_t j=0; j < M.cols()-1; ++j) file << M(i,j) << sep;
            file << M(i, M.cols()-1) <<  "\n";  
    }
    file.close();
}

template<typename T> void eigen2txt(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M, const std::string& filename = "mat.txt", bool append = false){
    eigen2ext<T>(M, " ", filename, append);
}

template<typename T> void eigen2csv(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& M, const std::string& filename = "mat.csv", bool append = false){
    eigen2ext<T>(M, ",", filename, append);
}

template< typename T> void vector2ext(const std::vector<T>& V, const std::string& sep, const std::string& filename, bool append = false){
    std::ofstream file;

    if(!append) 
        file.open(filename);
    else
        file.open(filename, std::ios_base::app);
    
    for(std::size_t i = 0; i < V.size()-1; ++i) file << V[i] << sep;
    
    file << V[V.size()-1] << "\n";  
    
    file.close();
}

template< typename T> void vector2txt(const std::vector<T>& V, const std::string& filename = "vec.txt", bool append = false){
   vector2ext<T>(V, " ", filename, append);
}

template< typename T> void vector2csv(const std::vector<T>& V, const std::string& filename = "vec.csv", bool append = false){
   vector2ext<T>(V, ",", filename, append);
}

void write_table(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M, const std::vector<std::string>& header = {}, const std::string& filename = "data.txt"){

    std::ofstream file(filename);

    if(header.empty() || header.size() != M.cols()){
        std::vector<std::string> head(M.cols());
        for(std::size_t i = 0; i < M.cols(); ++i)
                head[i] =  "V" + std::to_string(i);
        vector2txt<std::string>(head, filename);    
    }else vector2txt<std::string>(header, filename);
    
    eigen2txt<double>(M, filename, true);
}

void write_csv(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M, const std::vector<std::string>& header = {}, const std::string& filename = "data.csv"){
    std::ofstream file(filename);

    if(header.empty() || header.size() != M.cols()){
        std::vector<std::string> head(M.cols());
        for(std::size_t i = 0; i < M.cols(); ++i)
                head[i] =  "V" + std::to_string(i);
        vector2csv(head, filename);    
    }else vector2csv(header, filename);
    
    eigen2csv<double>(M, filename, true);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i + 1 != v.size()) os << " ";
    }
    return os;
}

auto noise(std::size_t nrow, std::size_t ncol, double sigma, std::mt19937 gen) {
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(nrow, ncol);
    std::normal_distribution<> __noise(0.0, sigma);
    for (std::size_t i = 0; i < nrow; ++i){
      for( std::size_t j = 0; j < ncol; ++j ){
       res(i, j) = __noise(gen); 
      }
    }
    return res;
}


template<class Scalar>
Eigen::SparseMatrix<Scalar> blockDiag(const std::vector<Eigen::SparseMatrix<Scalar>>& blocks) {
    const int m = static_cast<int>(blocks.size());
    const int n = blocks.empty() ? 0 : blocks[0].rows(); 
    Eigen::SparseMatrix<Scalar> M(m*n, m*n);

    std::vector<Eigen::Triplet<Scalar>> T;
    T.reserve([&]{
        size_t s = 0;
        for (auto& A : blocks) s += A.nonZeros();
        return s;
    }());

    for (int k = 0; k < m; ++k) {
        const auto& A = blocks[k];             
        for (int outer = 0; outer < A.outerSize(); ++outer) {
            for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(A, outer); it; ++it) {
                const int i = it.row();           
                const int j = it.col();     
                T.emplace_back(k*n + i, k*n + j, it.value());
            }
        }
    }
    M.setFromTriplets(T.begin(), T.end());
    return M;
}

template<typename T>
void range(const Eigen::SparseMatrix<T>& A){
    auto max = A.coeffs().maxCoeff();
    auto min = A.coeffs().minCoeff();
    std::cout << "matrix values: " << min << " " << max << std::endl; 
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
expand_grid(const std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>>& factors) {
    using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    const size_t K = factors.size();
    if (K == 0) return Mat(0, 0);

    // total rows = product of sizes
    size_t nrows = 1;
    std::vector<size_t> sz(K);
    for (size_t j = 0; j < K; ++j) {
        sz[j] = static_cast<size_t>(factors[j].size());
        if (sz[j] == 0) return Mat(0, static_cast<Eigen::Index>(K));
        nrows *= sz[j];
    }

    Mat M(static_cast<Eigen::Index>(nrows), static_cast<Eigen::Index>(K));

    // prefix products: prod_{t < j} sz[t]
    std::vector<size_t> prefix(K, 1);
    for (size_t j = 1; j < K; ++j) prefix[j] = prefix[j - 1] * sz[j - 1];

    for (size_t j = 0; j < K; ++j) {
        const auto& f = factors[j];
        const size_t nj = sz[j];
        const size_t repeat_len = prefix[j];                 // first factor repeats 1, next repeats sz[0], etc.
        const size_t block = repeat_len * nj;                // size of one full cycle for column j
        const size_t cycles = nrows / block;

        Eigen::Index row = 0;
        for (size_t c = 0; c < cycles; ++c) {
            for (size_t i = 0; i < nj; ++i) {
                for (size_t r = 0; r < repeat_len; ++r) {
                    M(row++, static_cast<Eigen::Index>(j)) = f(static_cast<Eigen::Index>(i));
                }
            }
        }
    }
    return M;
}
