#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
#include <format>
constexpr double tolerance = 1e-4;
const std::size_t iter = 128;

using key = std::pair<std::size_t, std::size_t>;
[[nodiscard]] inline SparseMatrix Unit(std::size_t size){
    std::map<key, double> e;
    for(std::size_t i = 0; i < size; ++i){
        e[std::make_pair(i, i)] = static_cast<double>(1);
    }
    return SparseMatrix(size, size, e);
}
TEST(chebyshov, b_zero){
    std::size_t n = 5;
    SparseMatrix unit = Unit(n);
    std::vector<double> bVec(n);
    Vector b(bVec);

    /*xBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = std::log(i + 1) * std::sin(i) * std::exp(i);
    }
    Vector xBegin(xBeginVec);
    auto pseudo_solve = Chebyshov(unit, b, xBegin, iter, std::make_pair(1, 1), tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(0), tolerance);
    }
}
TEST(chebyshov, with_fixed_mtx){
    std::size_t n = 2;
    /*fixed tau for this matrix*/
    std::map<key, double> mtxMap = {{{0,0}, 4}, {{0, 1}, 1}, {{1, 0}, 1}, {{1, 1}, 3}};
    SparseMatrix mtx(n, n, mtxMap);
    std::vector<double> bVec = {5, 4};;
    Vector b(bVec);

    /*vBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = std::log(i + 1) * std::sin(i);
    }
    Vector xBegin(xBeginVec);

    /*spectrum*/
    double l_max = 0.5 * (7 + std::sqrt(5));
    double l_min = 0.5 * (7 - std::sqrt(5));

    auto pseudo_solve = Chebyshov(mtx, b, xBegin, iter, std::make_pair(l_min, l_max), tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(1), tolerance);
    }
}
TEST(chebyshov, for_fixed_mtx){
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
    SparseMatrix mtx(n, n, mtxMap);
    std::vector<double> bVec = {5, 6, 5};
    Vector b(bVec);

    /*vBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i);
    }
    Vector xBegin(xBeginVec);
    double l_min = 4 - std::sqrt(2);
    double l_max = 4 + std::sqrt(2);

    auto pseudo_solve = Chebyshov(mtx, b, xBegin, iter, std::make_pair(l_min, l_max), tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(1), tolerance);
    }
}