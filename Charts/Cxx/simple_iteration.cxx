#include <Slae.hxx>
#include <chrono>
#include <fstream>
#include <random>
#include <functional>
using key = std::pair<std::size_t, std::size_t>;
/*alpha -  sparsity of diagonal sparse matrix, if alpha == 1, mtx is sparse, and dense in another case*/
inline SparseMatrix CreateSparse(std::size_t nx, double alpha, double min = 0, double max = 10){
    std::map<key, double> sparse;
    
    std::random_device dev;
    std::mt19937 r(dev());
    std::uniform_real_distribution<double> values(min, max);
    /*the probability that a random element is non-zero*/
    std::uniform_real_distribution<double> probability(0, 1);
    
    for(std::size_t i = 0; i < nx; ++i){
        double val = (probability(r) <= alpha) ? values(r) : 0;
        sparse[std::make_pair(i, i)] = val;
    }
    return SparseMatrix(nx, nx, sparse);
}
/*create xBegin - may be choosen any*/
inline Vector CreateBeginX(std::size_t size){
    std::vector<double> res(size);
    for(std::size_t i = 0; i < size; ++i){
        res[i] = std::log(i + 1) * std::exp(std::sin(i * 9)) + std::exp(-0.9 * i);
    }
    return Vector(res);
}
/*there are three main methods*/
std::size_t countMethods = 3;
const double tolerance = 1e-4;
const double tau = 0.9;
const std::size_t size = 100;
const std::size_t iter = 100;

/*measure time of solve slae*/
using Solver = std::pair<Vector, double>;
inline double IterationTime(const std::function<Solver ()>& function){
    auto start = std::chrono::steady_clock::now();
    auto solve = function();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> tau = end - start;
    return tau.count();
}
int main(){
    std::ofstream times("times.txt");
    for(std::size_t i = 1; i <= size; ++i){
        std::vector<double> b(i, static_cast<double>(0));
        SparseMatrix diag = CreateSparse(i, 0.9);
        Vector xBegin = CreateBeginX(i);
        auto sim = [&]() {return Solve(diag, Vector(b), xBegin, iter, tau, tolerance);};
        double simTau = IterationTime(sim);

        auto jacobi = [&]() {return Jacobi(diag, b, xBegin, iter);};
        double jacobiTau = IterationTime(jacobi);

        auto gz = [&]() {return GaussZeidel(diag, b, xBegin, iter);};
        double gzTau = IterationTime(gz);
        times << i << " " << simTau << " " << jacobiTau << " " << gzTau << std::endl;
    }
    return 0;
}