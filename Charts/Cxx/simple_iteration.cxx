#include <Slae.hxx>
#include <chrono>
#include <fstream>
#include <random>
#include <functional>
#include <iostream>
using key = std::pair<std::size_t, std::size_t>;
inline SparseMatrix CreateSparse(std::size_t nx){
    std::map<key, double> m;
    for(std::size_t i = 0; i < nx; ++i) {
        m[std::make_pair(i,i)] = 10.0;
        if(i > 0) m[std::make_pair(i, i - 1)] = -5.0;
        if(i < nx - 1) m[std::make_pair(i, i + 1)] = -5.0;
    }
    SparseMatrix mtx(nx, nx, m);
    return mtx;
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
const double tolerance = 1e-10;
const double tau = 0.8;
const std::size_t size = 100;
const std::size_t iter = 128;

/*measure time of solve slae*/
inline double IterationTime(const std::function<std::pair<Vector, double> ()>& function){
    auto start = std::chrono::steady_clock::now();
    auto solve = function();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> tau = end - start;
    return tau.count();
}
int main(){
    try{
        std::ofstream times("times.txt");
        for(std::size_t i = 1; i <= size; ++i){
            std::vector<double> b(i, static_cast<double>(std::sin(i * 0.3)));
            SparseMatrix mtx = CreateSparse(i);
            Vector xBegin = CreateBeginX(i);

            /*measuring times*/
            auto sim = [&]() {return Solve(mtx, Vector(b), xBegin, iter, tau, tolerance);};
            double simTau = IterationTime(sim);

            auto jacobi = [&]() {return Jacobi(mtx, b, xBegin, iter);};
            double jacobiTau = IterationTime(jacobi);

            auto gz = [&]() {return GaussZeidel(mtx, b, xBegin, iter);};
            double gzTau = IterationTime(gz);

            auto cheb = [&]() {return Tn(mtx, b, xBegin, iter, std::make_pair(4, 6), tolerance);};
            double chebTau = IterationTime(cheb);

            times << i << " " << simTau << " " << jacobiTau << " " << gzTau << std::endl;
        }
    }
    catch(std::exception& e){
        std::cout << e.what() << std::endl;
        }
    return 0;
}