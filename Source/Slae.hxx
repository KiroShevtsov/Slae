#pragma once
#include <vector>
#include <ostream>
#include <stdexcept>
#include <cmath>
#include <map>
class Vector{
    std::vector<double> v_;
public:
    std::size_t dim;
    Vector(const std::vector<double>& v) : v_(v), dim(v.size()){
        if(v.empty()){throw std::invalid_argument("Vector is empty");}
    }
    const std::vector<double>& data() const {return v_;}

    Vector operator+(const Vector& o) const {
        if(v_.size() != o.v_.size()) {throw std::invalid_argument("Impossible to sum vectors");}
        std::vector<double> res;
        res.reserve(v_.size());
        double sum = 0;
        for(std::size_t i = 0; i < v_.size(); ++i){
            sum = v_[i] + o.v_[i];
            res.push_back(sum);
        }
        return res;
    }

    Vector operator-(const Vector& o) const {
        if(v_.size() != o.v_.size()) {throw std::invalid_argument("Vectors have different dimensions");}
        std::vector<double> res;
        res.reserve(v_.size());
        for(std::size_t i = 0; i < v_.size(); ++i) {
            double r = v_[i] - o.v_[i];
            res.push_back(r);
        }
        return res;
    }

    Vector operator*(double a) const {
        std::vector<double> res;
        res.reserve(v_.size());
        for(const auto& item : v_){
            res.push_back(a * item);
        }
        return res;
    }
    inline double operator*(const Vector& o) const;

    double& operator[](std::size_t i) {return v_[i];}

    const double& operator[](std::size_t i) const {return v_[i];}
    bool operator==(const Vector& o) const {return v_ == o.v_;}
};
inline std::ostream& operator<<(std::ostream& os, const Vector& v){
    std::string whitesp = " ";
    for(const auto& item : v.data()){
        os << item << whitesp;
    }
    return os;
}
inline double EuclidNorm(const Vector& v) {
    if(v.dim == 0) {throw std::invalid_argument("vector empty");}
    return std::sqrt(v * v);
}
double Vector::operator*(const Vector& o) const {
    if(v_.size() != o.v_.size()) {throw std::invalid_argument("Impossible to mult vectors");}
    double res = 0;
    for(std::size_t i = 0; i < v_.size(); ++i){
        res += v_[i] * o.v_[i];
    }
    return res;
}
/*nx - num of columns, ny - rows*/
class Matrix{
    std::vector<double> mtx_;
public:
    const std::size_t nx_;
    const std::size_t ny_;
    Matrix(std::size_t nx, std::size_t ny, const std::vector<double>& mtx) : nx_(nx), ny_(ny){
        if (mtx.empty() or mtx.size() != nx_ * ny_){
            throw std::invalid_argument("Matrix is incorrect");
        }
        mtx_ = mtx;
    }
    const double& operator()(std::size_t i, std::size_t j) const {return mtx_[i * nx_ + j];}

    double& operator()(std::size_t i, std::size_t j) {return mtx_[i * nx_ + j];}

    const std::vector<double>& data() const {return mtx_;}

    Vector operator*(const Vector& v) const {
        std::size_t v_size = v.data().size();
        if(nx_ != v_size) {throw std::invalid_argument("Impossible to mult on vector");}
        std::vector<double> res;
        res.reserve(ny_);
        for(std::size_t i = 0; i < ny_; ++i){
            double current = 0;
            for(std::size_t k = 0; k < v_size; ++k){
                current += mtx_[i * nx_ + k] * v[k];
            }
            res.push_back(current);
        }
        return {res};
    }
};
inline std::ostream& operator<<(std::ostream& os, const Matrix& matrix){
    std::string whitesp = " ";
    for(std::size_t i = 0; i < matrix.ny_; ++i){
        for(std::size_t j = 0; j < matrix.nx_; ++j){
            os << matrix(i, j) << whitesp;
        }
        os << std::endl;
    }
    return os;
}
class SparseMatrix{
    using S = std::map<std::pair<std::size_t, std::size_t>, double>;
    std::vector<double> values_;
    std::vector<std::size_t> cols_, rows_;
public:
    const std::size_t nx_;
    const std::size_t ny_;
    
    SparseMatrix(std::size_t nx, std::size_t ny, const S& mtx) : nx_(nx), ny_(ny){
        /*rows(0) = 0*/
        rows_.push_back(0);
        std::size_t non_zero = 0;
        for(std::size_t i = 0; i < ny_; ++i){
            for(std::size_t j = 0; j < nx_; ++j){
                if(mtx.count(std::make_pair(i, j))){
                    values_.push_back(mtx.at(std::make_pair(i, j)));
                    cols_.push_back(j);
                    non_zero++;
                }
            }
            rows_.push_back(non_zero);
        }
    };
    double operator()(std::size_t i, std::size_t j) const {
        std::size_t start = rows_[i];
        std::size_t end = rows_[i + 1];
        double res = 0;
        for(std::size_t k = start; k < end; ++k){
            if(cols_[k] == j){
                res = values_[k];
                return res;
            }
        }
        return res;
    }

    Vector operator*(const Vector& v) const {
        std::size_t v_size = v.data().size();
        if(nx_ != v_size) {throw std::invalid_argument("Impossible to mult on vector");}
        std::vector<double> res;
        res.reserve(ny_);
        for(std::size_t i = 0; i < ny_; ++i){
            double current = 0;
            for(std::size_t k = rows_[i]; k < rows_[i + 1]; ++k){
                std::size_t q = cols_[k];
                current += values_[k] * v[q];
            }
            res.push_back(current);
        }
        return {res};
    }
};
/*Solve with sweep method*/
[[nodiscard]] Vector Solve(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const Vector& d);

/*qr decomposition*/
[[nodiscard]] std::pair<Matrix, Matrix> Qr(const Matrix& mtx);

/*solve with QR Ax = b*/
[[nodiscard]] Vector Solve(const Matrix& mtx, const Vector& b);

/*solve with sim method, sparse_mtx @ x = b with absolute error*/
[[nodiscard]] std::pair<Vector, double> Solve(const SparseMatrix& mtx, const Vector& b,
                                     const Vector& xBegin, std::size_t iter, double tau, double tol);

/*solve with jacobi method*/
[[nodiscard]] std::pair<Vector, double> Jacobi(const SparseMatrix& mtx, const Vector& b, const Vector& vBegin, std::size_t iter);

/*gauss-zeidel method*/
[[nodiscard]] std::pair<Vector, double> GaussZeidel(const SparseMatrix& mtx, const Vector& b, const Vector& vBegin, std::size_t iter);

/*solve with Chebyshov acceleration*/
[[nodiscard]] std::pair<Vector, double> Tn(const SparseMatrix& mtx, const Vector& b,
                                        const Vector& xBegin, std::size_t iter, std::pair<double, double> lambdas, double tolerance);