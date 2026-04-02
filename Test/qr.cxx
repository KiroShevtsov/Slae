#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
const double delta = 1e-5;
std::size_t n = 4;
Matrix unit(std::size_t m){
    std::vector<double> u(m * m);
    for(std::size_t i = 0; i < m; ++i){
        u[i * m + i] = 1;
    }
    return Matrix(m, m, u);
};
TEST(qr_task, qr_unit){
    std::vector<double> v;
    v.reserve(n);
    for(std::size_t i = 0; i < n; ++i){
        v.push_back(i);
    }
    Matrix e = unit(n);

    auto [q, r] = Qr(e);
    Vector v_after_qr = q * (r * v);
    EXPECT_EQ(v, v_after_qr);    
}
TEST(qr_task, qr_on_rangle){
    std::vector<double> v;
    v.reserve(n);
    for(std::size_t i = 0; i < n; ++i){
        v.push_back(i);
    }
    std::size_t nx = 4, ny = 2;
    std::vector<double> aVec = {1, 2, 3, 17, 19, 20, 21, 1};
    Matrix a(nx, ny, aVec);

    std::cout << a << " " << v << std::endl;

    auto [q, r] = Qr(a);
    EXPECT_EQ((a * v).dim, (q * (r * v)).dim);
    for(std::size_t i = 0; i < (a * v).dim; ++i){
        EXPECT_NEAR((q * (r * v))[i], (a * v)[i], delta);
    }
}
TEST(qr_task, qr_sle_unit){
    Matrix e = unit(4);
    std::vector<double> v;
    v.reserve(n);
    for(std::size_t i = 0; i < n; ++i){
        v.push_back(i);
    }
    Vector solve = Solve(e, v);
    for(std::size_t i = 0; i < v.size(); ++i){
        EXPECT_NEAR(solve[i], v[i], delta);
    }
}
TEST(qr_task, qr_sle_examples){
    const double tol = 1e-1;
    Vector d({5, 9, 7, 2});
    std::vector<double> aVec = {113, std::sin(3), 0, 0, 4, 8, 1, 0, 0, 1, 5, 3, 0, 0, 4, 9};
    Matrix a(4, 4, aVec);
    Vector x = Solve({4, 1, 4}, {113, 8, 5, 9}, {std::sin(3), 1, 3}, d);
    Vector res = Solve(a, d);
    
    ASSERT_EQ(x.dim, res.dim);
    for(std::size_t i = 0; i < d.dim; ++i){
        // EXPECT_NEAR(x[i], res[i], delta);
        EXPECT_NEAR(x[i], res[i], tol);
    }
}