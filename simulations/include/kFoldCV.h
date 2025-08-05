#ifndef __KFOLD_CROSS_VALIDATION_H_
#define __KFOLD_CROSS_VALIDATION_H_

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <numeric>
#include <stdexcept>
#include <algorithm>

class KFoldCV {
public:
    using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using index_t  = Eigen::Index;

    // Construct K folds from X (n x p) and Y (n x q).
    // If shuffle=true, rows are randomly permuted with the given seed before splitting.
    KFoldCV(const matrix_t& X,
            const matrix_t& Y,
            int K,
            bool shuffle = true,
            uint64_t seed = 42)
        : X_(X), Y_(Y), K_(K)
    {
        const index_t n = X_.rows();
        if (K_ < 2 || K_ > n) {
            throw std::invalid_argument("K must be in [2, n_samples].");
        }
        if (Y_.rows() != n) {
            throw std::invalid_argument("X and Y must have the same number of rows (samples).");
        }

        // Build index vector [0, 1, ..., n-1]
        std::vector<index_t> idx(static_cast<std::size_t>(n));
        std::iota(idx.begin(), idx.end(), index_t{0});

        // Optional shuffle
        if (shuffle) {
            std::mt19937_64 rng(seed);
            std::shuffle(idx.begin(), idx.end(), rng);
        }

        // Compute fold sizes
        const index_t base = n / K_;
        const index_t rem  = n % K_; // first 'rem' folds get one extra sample

        // Partition indices into K folds
        folds_.reserve(K_);
        index_t start = 0;
        for (int k = 0; k < K_; ++k) {
            const index_t sz = base + (k < rem ? 1 : 0);
            folds_.emplace_back();
            auto& f = folds_.back();
            f.reserve(static_cast<std::size_t>(sz));
            for (index_t i = 0; i < sz; ++i) {
                f.push_back(static_cast<int>(idx[static_cast<std::size_t>(start + i)]));
            }
            start += sz;
        }
    }

    // Number of folds
    int k() const noexcept { return K_; }

    // Access helpers for indices (optional but handy)
    const std::vector<int>& test_indices(int k) const { return at_(folds_, k); }
    std::vector<int> train_indices(int k) const {
        check_k_(k);
        const auto n = static_cast<int>(X_.rows());
        std::vector<int> tr;
        tr.reserve(n - static_cast<int>(folds_[k].size()));
        for (int i = 0; i < K_; ++i) if (i != k) {
            tr.insert(tr.end(), folds_[i].begin(), folds_[i].end());
        }
        return tr;
    }

    // Materialize matrices for a specific fold k (0-based)
    matrix_t X_train(int k) const { return take_rows_(X_, train_indices(k)); }
    matrix_t Y_train(int k) const { return take_rows_(Y_, train_indices(k)); }
    matrix_t X_test (int k) const { return take_rows_(X_, at_(folds_, k)); }
    matrix_t Y_test (int k) const { return take_rows_(Y_, at_(folds_, k)); }

private:
    int K_;
    const matrix_t& X_;
    const matrix_t& Y_;
    std::vector<std::vector<int>> folds_; // test folds only

    void check_k_(int k) const {
        if (k < 0 || k >= K_) throw std::out_of_range("Fold index out of range.");
    }

    static const std::vector<int>& at_(const std::vector<std::vector<int>>& v, int k) {
        if (k < 0 || k >= static_cast<int>(v.size())) throw std::out_of_range("Fold index out of range.");
        return v[static_cast<std::size_t>(k)];
    }

    // Utility: gather specified rows from M into a new matrix_t
    static matrix_t take_rows_(const matrix_t& M, const std::vector<int>& rows) {
        matrix_t out(static_cast<index_t>(rows.size()), M.cols());
        for (index_t i = 0; i < static_cast<index_t>(rows.size()); ++i) {
            out.row(i) = M.row(static_cast<index_t>(rows[static_cast<std::size_t>(i)]));
        }
        return out;
    }
};

#endif // __KFOLD_CROSS_VALIDATION_H_
