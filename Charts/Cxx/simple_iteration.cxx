#include <Slae.hxx>
#include <chrono>
#include <Utility.hxx>
#include <iostream>
using key = std::pair<std::size_t, std::size_t>;
inline Sparse CreateSparse(std::size_t nx){
    std::map<key, double> m;
    for(std::size_t i = 0; i < nx; ++i) {
        m[std::make_pair(i,i)] = 10.1;
        if(i > 0) m[std::pair(i, i - 1)] = -5.0;
        if(i < nx - 1) m[std::pair(i, i + 1)] = -5.0;
    }
    Sparse mtx(nx, nx, m);
    return mtx;
}
constexpr std::size_t iter = 128;
constexpr std::size_t size = 100;
constexpr double tolerance = 1e-11;
constexpr double tau = 0.09;
/*for symmetric qz*/
constexpr double rho = 0.9;
constexpr std::pair<double, double> lambdas = {0.105, 20.09};
int main() {
    auto [lMin, lMax] = lambdas;
    // std::ofstream plot("errors.txt");
    SlaeIo::Output plot("errors.txt");
    SlaeIo::Output times("times.txt");
    
    Sparse mtx = CreateSparse(size);
    Vector b = Vector(std::vector<double>(size, 10.0));
    Vector x0(std::vector<double>(size, 100));

    auto logger = [&plot](const std::string& name) {
        return [&plot, name](std::size_t it, double e) {
            plot.stream_ << name << " " << it << " " << e << "\n";
        };
    };
    
    auto l1 = SimpleIteration(mtx, b, x0, iter, tau, tolerance, logger("sim"));
    auto l2 = Chebyshov(mtx, b, x0, iter, {0.10, 20.09}, tolerance, logger("chebyshov"));
    auto l3 = GaussZeidel(mtx, b, x0, iter, tolerance, logger("gauss-zeidel"));
    auto l4 = Jacobi(mtx, b, x0, iter, tolerance, logger("jacobi"));
    auto l5 = Symmetric::GaussZeidel_S(mtx, b, x0, iter, rho, tolerance, logger("symmetric_gz"));

    double t_sim = 0, t_gz = 0, t_jacobi = 0, t_ch = 0, t_symmetric = 0;

    Vector s = x0;
    for (int i = 1; i <= iter; i++) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = SimpleIteration(mtx, b, s, i, tau, tolerance);
        auto end = std::chrono::steady_clock::now();
        s = v;
        t_sim += std::chrono::duration<double>(end - start).count();
        times.stream_ << "sim" << " " << t_sim << " " << e << "\n";
    }
    
    Vector g_z = x0;
    for (int i = 1; i <= iter; i++) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = GaussZeidel(mtx, b, g_z, i, tolerance);
        auto end = std::chrono::steady_clock::now();
        g_z = v;
        t_gz += std::chrono::duration<double>(end - start).count();
        times.stream_ << "gauss-zeidel" << " " << t_gz << " " << e << "\n";
    }
    
    Vector j = x0;
    for (int i = 1; i <= iter; i++) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = Jacobi(mtx, b, j, i, tolerance);
        auto end = std::chrono::steady_clock::now();
        j = v;
        t_jacobi += std::chrono::duration<double>(end - start).count();
        times.stream_ << "jacobi" << " " << t_jacobi << " " << e << "\n";
    }

    Vector ch = x0;
    for(int i = 1; i <= iter; ++i) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = Chebyshov(mtx, b, ch, i, lambdas, tolerance);
        auto end = std::chrono::steady_clock::now();
        ch = v;
        t_ch += std::chrono::duration<double>(end - start).count();
        times.stream_ << "chebyshov" << " " << t_ch << " " << e << "\n";
    }

    Vector symmetric_gz = x0;
    for(int i = 1; i <= iter; ++i) {
        auto start = std::chrono::steady_clock::now();
        auto [v, e] = Symmetric::GaussZeidel_S(mtx, b, symmetric_gz, i, rho, tolerance);
        auto end = std::chrono::steady_clock::now();
        symmetric_gz = v;
        t_symmetric += std::chrono::duration<double>(end - start).count();
        times.stream_ << "sym_gz" << " " << t_symmetric << " " << e << "\n";
    }
    return 0;
}