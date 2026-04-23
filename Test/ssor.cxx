#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
#include <cmath>
using key = std::pair<std::size_t, std::size_t>;
constexpr std::size_t s = 10;
constexpr std::size_t iter = 128;
constexpr double tolerance = 1e-10;
constexpr double rho = 0.92;
TEST(Ssor, Symmetric_Gauss_Zeidel) {
    std::map<key, double> m;
    for(std::size_t i = 0; i < s; ++i) {
        m[std::make_pair(i,i)] = 10;
        if(i > 0) m[std::pair(i, i - 1)] = -5.0;
        if(i < s - 1) m[std::pair(i, i + 1)] = -5.0;
    }
    Sparse mtx(s, s, m);
    Vector solve({3, 6, 8.8, 11.2, 13, 14, 14, 12.8, 10.2, 6});
    std::vector<double> xBeginStl(s);
    for(std::size_t i = 0; i < s; ++i) {
        xBeginStl[i] = 0;
    }
    std::vector<double> bVec(s);
    for(std::size_t i = 0; i < s; ++i) {
        bVec[i] = i;
    }
    auto [ssorSolve, error] = Symmetric::GaussZeidel_S(mtx, Vector(bVec), Vector(xBeginStl), iter, rho, tolerance);

    for(std::size_t i =  0; i < s; ++i) {EXPECT_NEAR(solve[i], ssorSolve[i], tolerance);}
    EXPECT_LE(error, tolerance);
}