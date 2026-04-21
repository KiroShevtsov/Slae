#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
#include <format>
constexpr double tau = 0.2;
/*required tolerance*/
constexpr double tolerance = 1e-10;
constexpr std::size_t iter = 100;
using key = std::pair<std::size_t, std::size_t>;
[[nodiscard]] inline Sparse Unit(std::size_t size){
    std::map<key, double> e;
    for(std::size_t i = 0; i < size; ++i){
        e[std::make_pair(i, i)] = static_cast<double>(1);
    }
    return Sparse(size, size, e);
}
TEST(sim, b_zero){
    std::size_t n = 5;
    Sparse unit = Unit(n);
    std::vector<double> bVec(n);
    Vector b(bVec);

    /*xBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = std::log(i + 1) * std::sin(i) * std::exp(i);
    }
    Vector xBegin(xBeginVec);
    auto [simSolve, error] = SimpleIteration(unit, b, xBegin, iter, tau, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", error) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(simSolve[i], static_cast<double>(0), tolerance);
    }
}
TEST(sim, sim_with_fixed_mtx){
    std::size_t n = 2;
    double tau_sim = 0.2;
    std::map<key, double> mtxMap = {{{0,0}, 4}, {{0, 1}, 1}, {{1, 0}, 1}, {{1, 1}, 3}};
    Sparse mtx(n, n, mtxMap);
    std::vector<double> bVec = {5, 4};;
    Vector b(bVec);

    /*vBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = std::log(i + 1) * std::sin(i);
    }
    Vector xBegin(xBeginVec);
    auto [simSolve, error] = SimpleIteration(mtx, b, xBegin, iter, tau_sim, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", error) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(simSolve[i], static_cast<double>(1), tolerance);
    }
}
TEST(sim, jacobi_for_zero_b){
    std::size_t n = 5;
    Sparse unit = Unit(n);
    std::vector<double> bVec(n);
    /*zero vector*/
    Vector b(bVec);

    /*xBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i) * std::exp(i);
    }
    Vector xBegin(xBeginVec);
    auto [jacobiSolve, error] = Jacobi(unit, b, xBegin, iter, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", error) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(jacobiSolve[i], static_cast<double>(0), tolerance);
    }
}
TEST(sim, jacobi_for_fixed_mtx){
    std::size_t n = 3;
    std::map<key, double> mtxMap;
    for(std::size_t i = 0; i < n; ++i){
        mtxMap[std::make_pair(i, i)] = static_cast<double>(4);
        if(i > 0) {
            mtxMap[std::make_pair(i, i - 1)] = static_cast<double>(1);
        }
        if (i < n - 1) {
            mtxMap[std::make_pair(i, i + 1)] = static_cast<double>(1);
        }
    }
    Sparse mtx(n, n, mtxMap);
    std::vector<double> bVec = {5, 6, 5};
    Vector b(bVec);

    /*vBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i);
    }
    Vector xBegin(xBeginVec);
    auto [jacobiSolve, error] = Jacobi(mtx, b, xBegin, iter, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", error) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(jacobiSolve[i], static_cast<double>(1), tolerance);
    }
}
TEST(sim, gauss_zeidel_for_zero_b){
    std::size_t n = 5;
    Sparse unit = Unit(n);
    std::vector<double> bVec(n);
    /*zero vector*/
    Vector b(bVec);

    /*xBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i) * std::exp(i);
    }
    Vector xBegin(xBeginVec);
    auto [solve, error] = GaussZeidel(unit, b, xBegin, iter, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", error) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(solve[i], static_cast<double>(0), tolerance);
    }
}
TEST(sim, gauss_zeidel_for_mtx){
    std::size_t n = 2;
    std::map<key, double> mtxMap = {{{0,0}, 4}, {{0, 1}, 1}, {{1, 0}, 1}, {{1, 1}, 3}};
    Sparse mtx(n, n, mtxMap);
    std::vector<double> bVec = {5, 4};
    Vector b(bVec);

    /*vBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i);
    }
    Vector xBegin(xBeginVec);
    auto [solve, error] = GaussZeidel(mtx, b, xBegin, iter, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", error) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(solve[i], static_cast<double>(1), tolerance);
    }
}