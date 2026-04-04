#include "Slae.hxx"
/*sim for sparse matrix: v[i + 1] = p @ v[i]  + b*/
std::pair<Vector, double> SimpleIteration(const SparseMatrix& p, const Vector& b, const Vector& vBegin, std::size_t iter, double tol){
    Vector v = vBegin;
    double delta = Norm(v);
    for(std::size_t l = 0; l < iter; ++l){
        Vector current = p * v + b;
        double e = Norm(current);
        if(e < tol){
            return {Vector(v), e};
        }
        delta = e;
        v = current;
    }
    return {v, delta};
}
[[nodiscard]] Vector Solve(const SparseMatrix& mtx, const Vector& b, const Vector& xBegin, std::size_t iter, double tol){}