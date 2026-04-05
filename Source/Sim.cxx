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
    double delta = Norm(args.a_ * v - args.b_);
    for(std::size_t l = 0; l < args.iter_; ++l){
        Vector vCurrent = v - (args.a_ * v - args.b_) * args.tau_;
        /*current r = Ax - b*/
        double error = Norm(args.a_ * vCurrent - args.b_);
        if(error < args.tolerance_){
            double absoluteDelta = std::abs(args.tolerance_ - error);
            return {v, absoluteDelta};
        }
        delta = error;
        v = vCurrent;
    }
    double absoluteDelta = std::abs(args.tolerance_ - delta);
    return {v, absoluteDelta};
}
[[nodiscard]] std::pair<Vector, double> Solve(const SparseMatrix& mtx, const Vector& b, 
                                    const Vector& xBegin, std::size_t iter, double tau, double tol){
    IterationArgs args(mtx, b, xBegin, iter, tol, tau);
    auto solve = SimpleIteration(args);
    return solve;
}
[[nodiscard]] std::pair<Vector, double> Jacobi(const SparseMatrix& mtx, const Vector b, const Vector& vBegin, std::size_t iter){
    if (mtx.nx_ != mtx.ny_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;
    /*anti-singular stabilisation*/
    const double lambda = 1e-3;
    for(std::size_t i = 0; i < iter; ++i){
        for(std::size_t k = 0; k < mtx.ny_; ++k){
            double diagonalMtx = mtx(k, k);
            if(diagonalMtx == 0) {diagonalMtx += lambda;}
            /*(l + u) @ x*/
            double lux = 0;
            for(std::size_t j = 0; j < mtx.ny_; ++j) {
                if (i != j) {lux += mtx(i, j) * v[j];}
            }
            v[k] = (1 / diagonalMtx) * (b[k] - lux);
        }
        delta = Norm(mtx * v - b);
    }
    return {v, delta};
}
[[nodiiscard]] std::pair<Vector, double> GaussZeidel(const SparseMatrix& mtx, const Vector& b, const Vector& vBegin, std::size_t iter){
    if (mtx.ny_ != mtx.nx_) {throw std::invalid_argument("matrix is not square");}
    Vector v = vBegin;
    double delta = 0;
    /*anti-singular stabilisation*/
    const double lambda = 1e-3;
    for(std::size_t i = 0; i < iter; ++i){
        for(std::size_t k = 0; k < mtx.ny_; ++k){
            double diagonalMtx = mtx(k, k);
            if(diagonalMtx == 0) {diagonalMtx += lambda;}
            double p = 0, q = 0;
            /*with x[i + 1]*/
            for(std::size_t j = 1; j < k; ++j) {p += mtx(k, j) * v[j];} 
            /*with x[i]*/
            for (std::size_t j = k + 1; j < mtx.ny_; ++j) {q += mtx(k, j) * v[j];}
            v[k] = (1 / diagonalMtx) * (b[k] - p - q);
        }
        delta = Norm(mtx * v - b);
    }
    return {v, delta};
}