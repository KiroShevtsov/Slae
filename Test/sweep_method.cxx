#include <gtest/gtest.h>
#include "../Source/Slae.hxx"
TEST(Sweep, SweepTest){
    std::vector<double> solve1 = {1, 4, 5};
    Vector  ksi = Solve({0, 0}, {1, 1, 1}, {0, 0}, solve1);
    ASSERT_EQ(ksi , Vector({1, 4, 5}));
    std::vector<double> solve = {
        115.0 / 297,     
        20.0 / 27,    
        4077.0 / 2673,      
        -1218.0 / 2673
        
    };
    Vector d({5, 9, 7, 2});
    Vector x = Solve({4, 1, 4}, {11, 8, 5, 9}, {1, 1, 3}, d);
    const double tolerance = 1e-3;
    for(std::size_t i = 0; i < d.dim; ++i){
        EXPECT_NEAR(x[i], solve[i], tolerance);
    }
}
