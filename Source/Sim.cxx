#include "Slae.hxx"
struct IterationArgs{
    SparseMatrix a_;
    Vector b_;
    Vector vBegin_;
    std::size_t iter_;
    double tolerance_;
    double tau_;
    IterationArgs(const SparseMatrix& m, const Vector& b, const Vector& vBegin, 
            std::size_t iter, double tolerance, double tau) : a_(m), b_(b), vBegin_(vBegin), iter_(iter), tolerance_(tolerance), tau_(tau){}
};
/*sim for sparse matrix: v = v - tau (a @ v - b)*/
std::pair<Vector, double> SimpleIteration(const IterationArgs& args){
    if(args.a_.nx_ != args.a_.ny_) {throw std::invalid_argument("matrix is not square");}
    Vector v = args.vBegin_;
    double delta = EuclidNorm(args.a_ * v - args.b_);
    for(std::size_t l = 0; l < args.iter_; ++l){
        Vector vCurrent = v - (args.a_ * v - args.b_) * args.tau_;
        /*current r = Ax - b*/
        double error = EuclidNorm(args.a_ * vCurrent - args.b_);
        if(error < args.tolerance_){
            double absoluteDelta = error;
            return {vCurrent, absoluteDelta};
        }
        delta = error;
        v = vCurrent;
    }
    double absoluteDelta = delta;
    return {v, absoluteDelta};
}
[[nodiscard]] std::pair<Vector, double> Solve(const SparseMatrix& mtx, const Vector& b, 
                                    const Vector& xBegin, std::size_t iter, double tau, double tol){
    IterationArgs args(mtx, b, xBegin, iter, tol, tau);
    auto solve = SimpleIteration(args);
    return solve;
}
[[nodiscard]] std::pair<Vector, double> Jacobi(const SparseMatrix& mtx, const Vector& b, const Vector& vBegin, std::size_t iter){
    if (mtx.nx_ != mtx.ny_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;
    /*anti-singular stabilisation*/
    const double lambda = 1e-3;
    for(std::size_t i = 0; i < iter; ++i){
        /*vector on i-iter*/
        Vector vCurrent = v;
        for(std::size_t k = 0; k < mtx.ny_; ++k){
            double diagonalMtx = mtx(k, k);
            if(diagonalMtx == 0) {diagonalMtx += lambda;}
            /*(l + u) @ x*/
            double lux = 0;
            for(std::size_t j = 0; j < mtx.nx_; ++j) {
                if (j != k) {lux += mtx(k, j) * vCurrent[j];}
            }
            vCurrent[k] = (1 / diagonalMtx) * (b[k] - lux);
        }
        v = vCurrent;
    }
    delta = EuclidNorm(mtx * v - b);
    return {v, delta};
}
[[nodiscard]] std::pair<Vector, double> GaussZeidel(const SparseMatrix& mtx, const Vector& b, const Vector& vBegin, std::size_t iter){
    if (mtx.ny_ != mtx.nx_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;
    /*anti-singular stabilisation*/
    const double lambda = 1e-3;
    for(std::size_t i = 0; i < iter; ++i){
        Vector vCurrent = v;
        for(std::size_t k = 0; k < mtx.ny_; ++k){
            double diagonalMtx = mtx(k, k);
            if(diagonalMtx == 0) {diagonalMtx += lambda;}
            double p = 0, q = 0;
            /*with x[i]*/
            for (std::size_t j = k + 1; j < mtx.nx_; ++j) {
                q += mtx(k, j) * vCurrent[j];
            }
            for(std::size_t j = 0; j < k; ++j) {
                p += mtx(k, j) * v[j];
            } 
            v[k] = (1 / diagonalMtx) * (b[k] - p - q);
        }
    }
    delta = EuclidNorm(mtx * v - b);
    return {v, delta};
}
[[nodiscard]] std::pair<Vector, double> Tn(const SparseMatrix& mtx, const Vector& b,
                                const Vector& xBegin, std::size_t iter, std::pair<double, double> lambdas ,double tolerance) {
    std::size_t r = std::log(iter) / std::log(2.0);
    std::vector<std::size_t> idx = {0};
    for(std::size_t i = 0; i < r; ++i) {
        std::size_t m = std::pow(2, i);
        
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
        // tau[i] = 1  / ti;
        double t_tilda = 0.5 * (lambdas.first + lambdas.second) + 0.5 * (lambdas.second - lambdas.first) * ti;
        tau[i] = 1 / t_tilda;
    }

    Vector v = xBegin;
    double delta = EuclidNorm(mtx * v - b);
    for(std::size_t l = 0; l < iter; ++l){
        Vector vCurrent = v - (mtx * v - b) * tau[l];
        /*current r = Ax - b*/
        double error = EuclidNorm(mtx * vCurrent - b);
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