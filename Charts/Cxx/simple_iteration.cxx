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