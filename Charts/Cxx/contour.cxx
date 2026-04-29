#include <iostream>
#include <Slae.hxx>
#include <Utility.hxx>
using key = std::pair<std::size_t, std::size_t>;
constexpr std::size_t iter = 50;
constexpr double tau = 0.110;
constexpr double tolerance = 1e-10;
int main() {
    Vector b({1, 1, 1});
    std::map<key, double> mtxMap = {{{0,0}, 10}, {{0, 1}, 3}, {{0, 2}, 6}, {{1, 0}, 3}, {{1, 1}, 5}, 
                                {{1, 2}, 1}, {{2, 0}, 6}, {{2, 1}, 1}, {{2, 2}, 8}};
    Sparse mtx(3, 3, mtxMap);
    Vector vBegin({0, 0, 0});
    auto [solve, error] = SimpleIteration(mtx, b * (-1), vBegin, iter, tau, tolerance);
    for(const auto& item : solve.data()) {std::cout << item << std::endl;}

    std::cout << solve << std::endl;
    Vector n({-0.14, 0.19, -0.53});
    double l = n.linalg_norm;

    /*normalized eigen vectors*/
    Vector vMin({0.64, -0.48, -0.58});
    Vector vMax({0.75, 0.26, 0.60});

    SlaeIo::Output file("proj.txt");
    std::vector<Vector> result;
    for(std::size_t i = 1; i <= iter; ++i) {
        auto [current, errorCurrent] = SimpleIteration(mtx, b * (-1), vBegin, i, tau, tolerance);
        double factor = (current * n) / (l * l);
        Vector currentProjection = n * factor;
        Vector stepProj = current - currentProjection;

        /*x = alpha @ vMin + beta @ vMax in this basis*/
        double alpha = stepProj * (vMin);
        double beta  = stepProj * (vMax);
        file.stream_ << alpha << " " << beta << std::endl;
    }
    return 0;
}