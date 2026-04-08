#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
#include <format>
const double tau = 0.2;
/*required tolerance*/
const double tolerance = 1e-4;
std::size_t iter = 100;
using key = std::pair<std::size_t, std::size_t>;
[[nodiscard]] SparseMatrix Unit(std::size_t size){
    std::map<key, double> e;
    for(std::size_t i = 0; i < size; ++i){
        e[std::make_pair(i, i)] = static_cast<double>(1);
    }
    return SparseMatrix(size, size, e);
}
TEST(sim, b_zero){
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
    auto pseudo_solve = Solve(unit, b, xBegin, iter, tau, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(0), tolerance);
    }
}
TEST(sim, sim_with_fixed_mtx){
    std::size_t n = 2;
    /*interesting :/*/
    double tau_sim = 0.2;
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
    auto pseudo_solve = Solve(mtx, b, xBegin, iter, tau_sim, tolerance);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(1), tolerance);
    }
}
TEST(sim, jacobi_for_zero_b){
    std::size_t n = 5;
    SparseMatrix unit = Unit(n);
    std::vector<double> bVec(n);
    /*zero vector*/
    Vector b(bVec);

    /*xBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i) * std::exp(i);
    }
    Vector xBegin(xBeginVec);
    auto pseudo_solve = Jacobi(unit, b, xBegin, iter);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(0), tolerance);
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
    SparseMatrix mtx(n, n, mtxMap);
    std::vector<double> bVec = {5, 6, 5};
    Vector b(bVec);

    /*vBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i);
    }
    Vector xBegin(xBeginVec);
    auto pseudo_solve = Jacobi(mtx, b, xBegin, iter);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(1), tolerance);
    }
}
TEST(sim, gauss_zeidel_for_zero_b){
    std::size_t n = 5;
    SparseMatrix unit = Unit(n);
    std::vector<double> bVec(n);
    /*zero vector*/
    Vector b(bVec);

    /*xBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i) * std::exp(i);
    }
    Vector xBegin(xBeginVec);
    auto pseudo_solve = GaussZeidel(unit, b, xBegin, iter);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(0), tolerance);
    }
}
TEST(sim, gauss_zeidel_for_mtx){
    std::size_t n = 2;
    std::map<key, double> mtxMap = {{{0,0}, 4}, {{0, 1}, 1}, {{1, 0}, 1}, {{1, 1}, 3}};
    SparseMatrix mtx(n, n, mtxMap);
    std::vector<double> bVec = {5, 4};
    Vector b(bVec);

    /*vBegin - may be choosen any*/
    std::vector<double> xBeginVec(n);
    for(std::size_t i = 0; i < n; ++i){
        xBeginVec[i] = i * std::sin(i);
    }
    Vector xBegin(xBeginVec);
    auto pseudo_solve = GaussZeidel(mtx, b, xBegin, iter);
    std::cout << "absolute error: " << std::format("{:.10f}", pseudo_solve.second) << std::endl;
    std::cout << "required error: " << tolerance << std::endl;
    for(std::size_t i = 0; i < n; ++i){
        EXPECT_NEAR(pseudo_solve.first[i], static_cast<double>(1), tolerance);
    }
}
TEST(sim, degenerate_matrix_no_solution) {
    std::size_t n = 2;
    std::map<key, double> mtxMap = {
        {{0, 0}, 0.0}, {{0, 1}, 0.0},
        {{1, 0}, 0.0}, {{1, 1}, 1.0}
    };
    SparseMatrix A(n, n, mtxMap);
    std::vector<double> bVec = {1.0, 0.0};
    Vector b(bVec);
    
    Vector xBegin(std::vector<double>(n, 0.0));
    
    auto [solution, error] = Solve(A, b, xBegin, iter, tau, tolerance);
    EXPECT_GE(error, tolerance);
}
TEST(sim, degenerate_matrix_inf_solutions) {
    std::size_t n = 2;
    std::map<key, double> mtxMap = {
        {{0, 0}, 0.0}, {{0, 1}, 0.0},
        {{1, 0}, 0.0}, {{1, 1}, 1.0}
    };
    SparseMatrix A(n, n, mtxMap);
    std::vector<double> bVec = {0.0, 0.0};
    Vector b(bVec);

    Vector xBegin(std::vector<double>(n, 0.0));
    auto [solution, error] = Solve(A, b, xBegin, iter, tau, tolerance);
    
    EXPECT_LE(error, tolerance);
    EXPECT_NEAR(solution[1], 0.0, tolerance);
}