#include <Slae.hxx>
#include <fstream>
#include <random>
#include <iostream>
/*alpha - is density of matrix*/
std::vector<double> Mtx(std::size_t nx, std::size_t ny, double alpha, double min = 0, double max = 10){
    std::size_t size = nx * ny;
    std::vector<double> res(size);
    
    std::random_device dev;
    std::mt19937 r(dev());
    std::uniform_real_distribution<double> values(min, max);
    std::uniform_real_distribution<double> prob_non_zero(0, 1);
    
    for(std::size_t i = 0; i < size; ++i){
        double l = prob_non_zero(r);
        double v = (l <= alpha) ? values(r) : 0;
        res[i] = v;
    }
    return res;
}
int main(){
    std::size_t n = 4;
    std::vector<double> mtxVec = Mtx(n, n, 1, 0, 1);
    Matrix mtx(n, n, mtxVec);
    auto [q, r] = Qr(mtx);
    std::ofstream file("data.txt");
    file << "X" << "\n" << mtx;
    file << "Q" << "\n" << q;
    file << "R" << "\n" << r;
    return 0;
}