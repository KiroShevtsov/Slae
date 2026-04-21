#include "Slae.hxx"
#include <utility>
/*sim for sparse matrix: v = v - tau (a @ v - b)*/
[[nodiscard]] std::pair<Vector, double> SimpleIteration(const Sparse& mtx, const Vector& b, 
                    const Vector& xBegin, std::size_t iter, double tau, double tolerance, 
                            const std::function<void(std::size_t, double)>& c){
    if(mtx.nx_ != mtx.ny_) {throw std::invalid_argument("matrix is not square");}
    Vector v = xBegin;
    double delta = 0;
    for(std::size_t l = 0; l < iter; ++l){
        Vector vCurrent = v - (mtx * v - b) * tau;
        /*current r = Ax - b*/
        double error = (mtx * vCurrent - b).linalg_norm;
        /*callback*/
        if(c) { c(l, error); }
        
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

[[nodiscard]] std::pair<Vector, double> Jacobi(const Sparse& mtx, const Vector& b, const Vector& vBegin, 
                std::size_t iter, double tolerance, const std::function<void(std::size_t, double)>& c){
    if (mtx.nx_ != mtx.ny_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;

    for(std::size_t i = 0; i < iter; ++i){
        /*vector on i-iter*/
        Vector vCurrent = v;
        /*current error*/
        double error = (mtx * vCurrent - b).linalg_norm;

        /*callback*/
        if(c) { c(i, error); }

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
[[nodiscard]] std::pair<Vector, double> GaussZeidel(const Sparse& mtx, const Vector& b, const Vector& vBegin,
                    std::size_t iter, double tolerance, const std::function<void(std::size_t, double)>& c){
    if (mtx.ny_ != mtx.nx_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;
    
    for(std::size_t i = 0; i < iter; ++i){
        /*current error*/
        double error = (mtx * v - b).linalg_norm;

        /*callback*/
        if (c) { c(i, error); };

        if(error < tolerance) {
            double absoluteDelta = error;
            return {v, absoluteDelta};
        }
        
        for(std::size_t k = 0; k < mtx.ny_; ++k){
            double diagonalMtx = mtx(k, k);
            if(diagonalMtx == 0) {diagonalMtx += singular;}
            double p = 0;
            for (std::size_t j = 0; j < mtx.nx_; ++j) {
                if (j != k) {p += mtx(k, j) * v[j];}
            } 
            v[k] = (1 / diagonalMtx) * (b[k] - p);
        }
        delta = error;
    }
    double absoluteDelta = delta;
    return {v, absoluteDelta};
}
[[nodiscard]] std::pair<Vector, double> Chebyshov(const Sparse& mtx, const Vector& b, const Vector& xBegin, 
                    std::size_t iter, const std::pair<double, double>& lambdas, double tolerance, 
                            const std::function<void(std::size_t, double)>& c) {
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
    auto [lMin, lMax] = lambdas;
    std::vector<double> tau(iter);
    for(std::size_t i = 0; i < iter; ++i) {
    /*create Chebyshov polynomial roots*/
        double arg = std::numbers::pi_v<double> * (2 * idx[i] + 1) / (2 * iter);
        double ti = std::cos(arg);
        double t_tilda = 0.5 * (lMin + lMax) + 0.5 * (lMax - lMin) * ti;
        tau[i] = 1 / t_tilda;
    }

    Vector v = xBegin;
    double delta = 0;
    for(std::size_t l = 0; l < iter; ++l){
        Vector vCurrent = v - (mtx * v - b) * tau[l];
        /*current r = Ax - b*/
        double error = (mtx * vCurrent - b).linalg_norm;
        /*callback*/
        if(c) { c(l, error); }

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
std::pair<Vector, double> Symmetric::GaussZeidel_S(const Sparse& mtx, const Vector& b, const Vector& vBegin, 
                std::size_t iter, double rho, const double& tolerance,
                    const std::function<void(std::size_t, double)>& c) {
    Vector v = vBegin;
    double delta = (mtx * vBegin - b).linalg_norm;

    auto fx = [&mtx, &b, &vBegin](const Vector& x) -> Vector {
        Vector u = x;
        std::size_t s = mtx.nx_;
        for(std::size_t i = 0; i < s; ++i) {
            double current = 0;
            for(std::size_t l = 0; l < i; ++l) {
                current += mtx(i, l) * u[l];
            }
            for(std::size_t l = i + 1; l < s; ++l) {
                current += mtx(i, l) * x[l];
            }
            double diagonalMtx = mtx(i, i);
            if (diagonalMtx == 0) {diagonalMtx += singular;}
            u[i] = (b[i] - current) / (diagonalMtx);
        }
        Vector res = u;
        for(long long i = static_cast<long long>(s) - 1; i >= 0; --i) {
            double current = 0;
            for(std::size_t j = 0; j < static_cast<std::size_t>(i); ++j) {
                current += mtx(i, j) * u[j];
            }
            for(std::size_t j = static_cast<std::size_t>(i) + 1; j < s; ++j) {
                current += mtx(i, j) * res[j];
            }
            double diagonalMtx = mtx(i, i);
            if(diagonalMtx == 0) {diagonalMtx += singular;}
            res[i] = (b[i] - current) / (diagonalMtx);
        }
        return res;
    };
    
    Vector yBegin = vBegin;
    /*first iteration symmetric Gauss-Zeidel*/
    Vector res = fx(yBegin);
    double w = 1;
    double wNext = 2 / (2 - rho * rho);

    for(std::size_t i = 2; i < iter + 1; ++i) {
        Vector z = fx(res);
        Vector current = (z - yBegin) * wNext + yBegin;
        yBegin = std::exchange(res, current);
        /*w[i]*/
        w = wNext;
        /*w[i + 1]*/
        wNext = 1 / (1 - 0.25 * (w * rho * rho));
        double error = (mtx * res - b).linalg_norm;

        /*callback*/
        if (c) {c(i, error);}

        if (error < tolerance) {
            double absoluteDelta = error;
            return {res, absoluteDelta};
        }
        delta = error;
    }
    double absoluteDelta = delta;
    return {res, absoluteDelta};
}