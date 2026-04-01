#include "Slae.hxx"
#include <cmath>
Vector Tridiagonal(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const Vector& d){
    std::size_t n = d.dim;
    if(b.size() != n or a.size() != (n - 1) or c.size() != (n - 1)){
        throw std::invalid_argument("Vectors don't match tridiagonal matrix elements");
    }
    std::vector<double> alpha(n - 1);
    std::vector<double> beta(n);

    alpha[0] = -(c[0] / b[0]);
    beta[0] = d[0] / b[0];
    for(std::size_t i = 1; i < (n - 1); ++i){
        alpha[i]  = -(c[i]) / (a[i - 1] * alpha[i - 1] + b[i]);
        beta[i] = (d[i] - beta[i - 1] * a[i - 1]) / (a[i - 1] * alpha[i - 1] + b[i]);
    }
    double q = (a[n - 2] * alpha[n - 2] + b[n - 1]);
    beta[n - 1] = (d[n - 1] - a[n - 2] * beta[n - 2]) / q;

    std::vector<double> x(n);
    x[n - 1] = beta[n - 1];
    for(int i = (n - 2); i >= 0; i--){
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    return Vector(x);
}
[[nodiscard]] Vector Solve(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const Vector& d){
    return Tridiagonal(a, b, c, d);
}
Vector Hausholder(const Vector& v, const Vector& x){
    double factor = -2 * (v * x) / (v * v);
    return x + v * factor;
}
[[nodiscard]] std::pair<Matrix, Matrix> Qr(const Matrix& mtx){
    if(mtx.data().empty()){throw std::invalid_argument("matrix is empty");}
    Matrix r = mtx;
    std::size_t rows = mtx.ny_;
    std::size_t cols = mtx.nx_;

    std::vector<double> qVec(rows * rows);
    for(std::size_t i = 0; i < rows; ++i){
        qVec[i * rows + i] = 1;
    }
    Matrix q(rows, rows, qVec);

    for(std::size_t l = 0; l < std::min(cols, rows); ++l){
        std::vector<double> xVec;
        for(std::size_t f = l; f < rows; ++f){
            xVec.push_back(mtx(f, l));
        }
        Vector x(xVec);
        double norm = Norm(x);
        double sign = (x[0] > 0) ? -norm : norm;
        std::vector<double> vVec(x.dim);
        vVec[0] = x[0] + sign;
        for(std::size_t i = 1; i < x.dim; ++i) {
            vVec[i] = x[i];
        }
        Vector v(vVec);
        if(v * v == 0){continue;}
        for(std::size_t j = l; j < cols; ++j) {
            std::vector<double> colVec;
            for(std::size_t i = l; i < rows; ++i) {
                colVec.push_back(mtx(i, j));
            }
            Vector column_j(colVec);
            Vector newColumn_j = Hausholder(v, column_j);
            for(std::size_t i = 0; i < newColumn_j.dim; ++i) {
                r(l + i, j) = newColumn_j[i];
            }
        }
        for(std::size_t i = 0; i < q.ny_; ++i) {
            std::vector<double> rowVec;
            for(std::size_t j = l; j < q.nx_; ++j) {
                rowVec.push_back(q(i, j));
            }
            Vector row_j(rowVec);
            Vector newRow_j = Hausholder(v, row_j);
            for(std::size_t j = 0; j < newRow_j.dim; ++j) {
                q(i, j + l) = newRow_j[j];
            }
        }
    }
    return {q, r};
}
[[nodiscard]] Vector Solve(const Matrix& mtx, const Vector& b){
    auto [q, r] = Qr(mtx);
    /*creating Q^Tb*/
    std::vector<double> y;
    y.reserve(q.ny_);
    for(std::size_t j = 0; j < q.ny_; ++j){
        double qTb = 0;
        for(std::size_t i = 0; i < q.nx_; ++i){
            qTb += q(i, j) * b[j];
        }
        y.push_back(qTb);
    }
    std::size_t n = b.dim;
    std::vector<double> x(n);
    /*beginnig from the end*/
    x[n - 1] = b[n - 1] / r(n - 1, n - 1);
    for(std::size_t l = n - 1; l > 0; --l){
        double linearCombination = 0;
        for(std::size_t j = l; j < n; ++j){
            linearCombination += r(l, j) * x[j];
        }
        double solveX_i = (y[l] - linearCombination) / r(l, l);
        x[l] = solveX_i;
    }
    return Vector(x);
}