#include <Slae.hxx>
#include <chrono>
#include <fstream>
#include <iostream>
using key = std::pair<std::size_t, std::size_t>;
inline SparseMatrix CreateSparse(std::size_t nx){
    std::map<key, double> m;
    for(std::size_t i = 0; i < nx; ++i) {
        m[std::make_pair(i,i)] = 10.1;
        if(i > 0) m[std::pair(i, i - 1)] = -5.0;
        if(i < nx - 1) m[std::pair(i, i + 1)] = -5.0;
    }
    SparseMatrix mtx(nx, nx, m);
    return mtx;
}
constexpr std::size_t iter = 128;
constexpr std::size_t size = 100;
constexpr double tolerance = 1e-11;
constexpr double tau = 0.09;

int main() {
    std::ofstream plot("errors.txt");
    std::ofstream times("times.txt");
    
    SparseMatrix mtx = CreateSparse(size);
    Vector b = Vector(std::vector<double>(size, 10.0));
    Vector x0(std::vector<double>(size, 100));

    auto logger = [&plot](const std::string& name) {
        return [&plot, name](std::size_t it, double e) {
            plot << name << " " << it << " " << e << "\n";
        };
    };
    
    auto l1 = Solve(mtx, b, x0, iter, tau, tolerance, logger("sim"));
    auto l2 = Chebyshov(mtx, b, x0, iter, {0.10, 20.09}, tolerance, logger("chebyshov"));
    auto l3 = GaussZeidel(mtx, b, x0, iter, tolerance, logger("gauss-zeidel"));
    auto l4 = Jacobi(mtx, b, x0, iter, tolerance, logger("jacobi"));

    double t_sim = 0, t_gz = 0, t_jacobi = 0, t_ch = 0;

    Vector s = x0;
    for (int i = 0; i <= iter; i++) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = Solve(mtx, b, s, iter, tau, tolerance);
        auto end = std::chrono::steady_clock::now();
        s = v;
        t_sim += std::chrono::duration<double>(end - start).count();
        times << "sim" << " " << t_sim << " " << e << "\n";
    }
    
    Vector g_z = x0;
    for (int i = 0; i <= iter; i++) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = GaussZeidel(mtx, b, g_z, iter, tolerance);
        auto end = std::chrono::steady_clock::now();
        g_z = v;
        t_gz += std::chrono::duration<double>(end - start).count();
        times << "gauss-zeidel" << " " << t_gz << " " << e << "\n";
    }
    
    Vector j = x0;
    for (int i = 0; i < iter; i++) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = Jacobi(mtx, b, j, iter, tolerance);
        auto end = std::chrono::steady_clock::now();
        j = v;
        t_jacobi += std::chrono::duration<double>(end - start).count();
        times << "jacobi" << " " << t_jacobi << " " << e << "\n";
    }

    Vector ch = x0;
    for(int i = 0; i < iter; ++i) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = Chebyshov(mtx, b, ch, iter, {0.10, 20.09}, tolerance);
        auto end = std::chrono::steady_clock::now();
        ch = v;
        t_ch += std::chrono::duration<double>(end - start).count();
        times << "chebyshov" << " " << t_ch << " " << e << "\n";
    }
    return 0;
}