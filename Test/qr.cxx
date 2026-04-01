#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
TEST(qr_task, qr_decomposition){
    const double delta = 1e-5;
    std::size_t n = 4;
    std::vector<double> v;
    v.reserve(n);
    auto unit = [&](std::size_t m)->Matrix{
        std::vector<double> u(m * m);
        for(std::size_t i = 0; i < m; ++i){
            u[i * m + i] = 1;
        }
        return Matrix(m, m, u);
    };
    for(std::size_t i = 0; i < n; ++i){
        v.push_back(i);
    }
    Matrix e = unit(n);

    auto [q, r] = Qr(e);
    Vector v_after_qr = q * (r * v);
    EXPECT_EQ(v, v_after_qr);

    std::size_t nxa = 4, nya = 2;
    std::vector<double> aVec = {1.0, 2.0, 3.0, 17.0, 19.0, 20.0, 23.0, 1.0};
    Matrix a(nxa, nya, aVec);
    auto [qa, ra] = Qr(a);
    Vector pseudo = qa * (ra * v);
    Vector res = a * v;
    EXPECT_EQ(res.dim, pseudo.dim);
    for(std::size_t i = 0; i < res.dim; ++i){
        EXPECT_NEAR(res[i], pseudo[i], delta);
    }

    std::size_t nxb = 4, nyb = 1;
    std::vector<double> bVec = {0, 1, 1, -1};
    Matrix b(nxb, nyb, bVec);
    auto [qb, rb] = Qr(b);
    /*(qr) * v = 0*/
    EXPECT_NEAR(0, (qb * (rb * v))[0], delta);
}