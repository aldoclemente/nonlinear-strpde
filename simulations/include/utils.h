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

