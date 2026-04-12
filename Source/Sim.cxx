#include "Slae.hxx"
/*sim for sparse matrix: v = v - tau (a @ v - b)*/
[[nodiscard]] std::pair<Vector, double> Solve(const SparseMatrix& mtx, const Vector& b, const Vector& xBegin, 
                                    std::size_t iter, double tau, double tolerance, const L& callback){
    if(mtx.nx_ != mtx.ny_) {throw std::invalid_argument("matrix is not square");}
    Vector v = xBegin;
    double delta = 0;
    for(std::size_t l = 0; l < iter; ++l){
        Vector vCurrent = v - (mtx * v - b) * tau;
        /*current r = Ax - b*/
        double error = (mtx * vCurrent - b).linalg_norm;
        /*callback*/
        if(callback) { callback(l, error); }
        
        if(error < tolerance){
            double absoluteDelta = error;
            return {vCurrent, absoluteDelta};
        }

        delta = error;
        v = vCurrent;
    }
    double absoluteDelta = delta;
    return {v, absoluteDelta};
}
/*resolve zeros on diagonal elements*/
constexpr double singular = 1e-7;

[[nodiscard]] std::pair<Vector, double> Jacobi(const SparseMatrix& mtx, const Vector& b, 
                                const Vector& vBegin, std::size_t iter, double tolerance, const L& callback){
    if (mtx.nx_ != mtx.ny_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;

    for(std::size_t i = 0; i < iter; ++i){
        /*vector on i-iter*/
        Vector vCurrent = v;
        /*current error*/
        double error = (mtx * vCurrent - b).linalg_norm;

        /*callback*/
        if(callback) { callback(i, error); }

        if(error < tolerance) {
            double absoluteDelta = error;
            return {vCurrent, absoluteDelta};
        }
        for(std::size_t k = 0; k < mtx.ny_; ++k){
            double diagonalMtx = mtx(k, k);
            if(diagonalMtx == 0) {diagonalMtx += singular;}
            /*(l + u) @ x*/
            double lux = 0;
            for(std::size_t j = 0; j < mtx.nx_; ++j) {
                if (j != k) {lux += mtx(k, j) * v[j];}
            }
            vCurrent[k] = (1 / diagonalMtx) * (b[k] - lux);
        }
        delta = error;
        v = vCurrent;
    }
    double absoluteDelta = delta;
    return {v, absoluteDelta};
}
[[nodiscard]] std::pair<Vector, double> GaussZeidel(const SparseMatrix& mtx, const Vector& b, 
                                const Vector& vBegin, std::size_t iter, double tolerance, const L& callback){
    if (mtx.ny_ != mtx.nx_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;
    
    for(std::size_t i = 0; i < iter; ++i){
        /*current error*/
        double error = (mtx * v - b).linalg_norm;

        /*callback*/
        if (callback) { callback(i, error); };

        if(error < tolerance) {
            double absoluteDelta = error;
            return {v, absoluteDelta};
        }

        for(std::size_t k = 0; k < mtx.ny_; ++k){
            double diagonalMtx = mtx(k, k);
            if(diagonalMtx == 0) {diagonalMtx += singular;}
            double p = 0;
            for (std::size_t j = 0; j < mtx.nx_; ++j) {
                if (j != k) p += mtx(k, j) * v[j];
            } 
            v[k] = (1 / diagonalMtx) * (b[k] - p);
        }
        delta = error;
    }
    double absoluteDelta = delta;
    return {v, absoluteDelta};
}
[[nodiscard]] std::pair<Vector, double> Chebyshov(const SparseMatrix& mtx, const Vector& b, const Vector& xBegin, 
                    std::size_t iter, std::pair<double, double> lambdas, double tolerance, const L& callback) {
    std::size_t w = iter;
    std::size_t r = 0; while (w >>= 1) {++r;};
    std::vector<std::size_t> idx = {0};
    for(std::size_t i = 0; i < r; ++i) {
        /*2 ^ i*/
        std::size_t m = 1 << i;
        
        std::vector<std::size_t> c;
        for(const auto& k : idx) {
            c.push_back(k);
            c.push_back(2 * m - k - 1);
        }
        idx = c;
    }
    
    /*create Chebyshov polynomial roots*/
    auto lambda = [&] (std::size_t n, std::size_t s) -> double {
        double arg = std::numbers::pi_v<double> * (2 * s + 1) / (2 * n);
        return std::cos(arg);
    };

    std::vector<double> tau(iter);
    for(std::size_t i = 0; i < iter; ++i) {
        double ti = lambda(iter, idx[i] );
        double t_tilda = 0.5 * (lambdas.first + lambdas.second) + 0.5 * (lambdas.second - lambdas.first) * ti;
        tau[i] = 1 / t_tilda;
    }

    Vector v = xBegin;
    double delta = 0;
    for(std::size_t l = 0; l < iter; ++l){
        Vector vCurrent = v - (mtx * v - b) * tau[l];
        /*current r = Ax - b*/
        double error = (mtx * vCurrent - b).linalg_norm;
        /*callback*/
        if(callback) { callback(l, error); }

        if(error < tolerance){
            double absoluteDelta = error;
            return {vCurrent, absoluteDelta};
        }
        delta = error;
        v = vCurrent;
    }
    double absoluteDelta = delta;
    return {v, absoluteDelta};
}