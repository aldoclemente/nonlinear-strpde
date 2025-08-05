#include <Eigen/Dense>
#include <iostream>

#include "include/kFoldCV.h"

int main() {
    // Toy data: 10 samples, 3 features; targets with 1 column
    Eigen::MatrixXd X = Eigen::MatrixXd::Random(11, 3);
    
    Eigen::MatrixXd Y = Eigen::MatrixXd::Random(11, 1);

    std::cout << X << std::endl;
    std::cout << Y << std::endl;
    

    KFoldCV kfold(X, Y, /*K=*/3, /*shuffle=*/true, /*seed=*/123);

    for (int k = 0; k < kfold.k(); ++k) {
        std::cout << "Fold " << k << ":\n";
        std::cout << "  X_train: " << kfold.X_train(k).rows() << " x " << kfold.X_train(k).cols() << "\n";
        std::cout << kfold.X_train(k) << std::endl;

        
        std::cout << "  X_test : "  << kfold.X_test(k).rows()  << " x " << kfold.X_test(k).cols()  << "\n";
        std::cout << kfold.X_test(k) << std::endl;
        
        std::cout << "  Y_train: " << kfold.Y_train(k).rows() << " x " << kfold.Y_train(k).cols() << "\n";
        std::cout << kfold.Y_train(k) << std::endl;
        
        std::cout << "  Y_test : "  << kfold.Y_test(k).rows()  << " x " << kfold.Y_test(k).cols()  << "\n";
        std::cout << kfold.Y_test(k) << std::endl;
        
    }

    {   
        std::cout << "\n" << std::endl;
        std::cout << "\n" << std::endl;
        
        Eigen::MatrixXd X = Eigen::MatrixXd::Ones(11, 3);
        X.block(3,0,3,3) = 2* X.block(3,0,3,3);
        
        X.block(6,0,3,3) = 3* X.block(6,0,3,3);
        
        Eigen::MatrixXd Y = Eigen::MatrixXd::Random(11, 1);

        std::cout << X << std::endl;
        std::cout << Y << std::endl;
    
        KFoldCV kfold(X, Y, /*K=*/5, /*shuffle=*/false, /*seed=*/123);

    for (int k = 0; k < kfold.k(); ++k) {
        std::cout << "Fold " << k << ":\n";
        std::cout << "  X_train: " << kfold.X_train(k).rows() << " x " << kfold.X_train(k).cols() << "\n";
        std::cout << kfold.X_train(k) << std::endl;

        
        std::cout << "  X_test : "  << kfold.X_test(k).rows()  << " x " << kfold.X_test(k).cols()  << "\n";
        std::cout << kfold.X_test(k) << std::endl;
        
        std::cout << "  Y_train: " << kfold.Y_train(k).rows() << " x " << kfold.Y_train(k).cols() << "\n";
        std::cout << kfold.Y_train(k) << std::endl;
        
        std::cout << "  Y_test : "  << kfold.Y_test(k).rows()  << " x " << kfold.Y_test(k).cols()  << "\n";
        std::cout << kfold.Y_test(k) << std::endl;
        
    }
    }

    {   
        
        std::cout << "\n" << std::endl;
        std::cout << "\n" << std::endl;
        std::cout << "LOOCV " << std::endl;
        std::cout << X << std::endl;
        std::cout << Y << std::endl;

        
        KFoldCV kfold(X, Y, /*K=*/10, /*shuffle=*/false, /*seed=*/123);

    for (int k = 0; k < kfold.k(); ++k) {
        std::cout << "Fold " << k << ":\n";
        std::cout << "  X_train: " << kfold.X_train(k).rows() << " x " << kfold.X_train(k).cols() << "\n";
        std::cout << kfold.X_train(k) << std::endl;

        
        std::cout << "  X_test : "  << kfold.X_test(k).rows()  << " x " << kfold.X_test(k).cols()  << "\n";
        std::cout << kfold.X_test(k) << std::endl;
        
        std::cout << "  Y_train: " << kfold.Y_train(k).rows() << " x " << kfold.Y_train(k).cols() << "\n";
        std::cout << kfold.Y_train(k) << std::endl;
        
        std::cout << "  Y_test : "  << kfold.Y_test(k).rows()  << " x " << kfold.Y_test(k).cols()  << "\n";
        std::cout << kfold.Y_test(k) << std::endl;
        
    }
    }
}
