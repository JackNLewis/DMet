//
// Created by jack lewis on 06/02/2022.
//
#include "gtest/gtest.h"
#include <fstream>
#include <vector>
#include <mpfr.h>
#include <DMet/PointDistances.h>


using std::cout;
using std::endl;
using std::vector;

////================== Chebyshev TESTS ===================//
TEST(ChebyshevTests, WorkingPoints) {
    //test with there being a difference of epsilon and producing a value that isn't zero
    vector<double> v1 {1.0,50.0,1.0};
    vector<double> v2 {1.0,0.0,51.0};
    vector<double> v3 {0.0,0.0,1.0};
    vector<double> v4 {-7.0,0.0,1.0};

    mpfr_t res;
    mpfr_init(res);
    DMet::PointDistances::getChebyshev(res,v1,v2);
    double res1 = mpfr_get_d(res,GMP_RNDN);

    DMet::PointDistances::getChebyshev(res,v3,v4);
    double res2 = mpfr_get_d(res,GMP_RNDN);

    mpfr_clear(res);

    EXPECT_DOUBLE_EQ(50.0, res1);
    EXPECT_DOUBLE_EQ(7.0, res2);
}
