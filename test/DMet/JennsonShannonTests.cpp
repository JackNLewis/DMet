//
// Created by jack lewis on 14/02/2022.
//

#include "gtest/gtest.h"
#include <fstream>
#include <vector>
#include <mpfr.h>
#include <DMet/DistribDistance.h>


using std::cout;
using std::endl;
using std::vector;


TEST(JensenShannon, WorkingPoints) {
    //test with there being a difference of epsilon and producing a value that isn't zero
    vector<double> v1 {0.9,0.05,0.05};
    vector<double> v2 {0.1,0.8,0.1};

    mpfr_t res;
    mpfr_init(res);

    //measure time compared to scipy
    //measure over 1000 examples
    //standard deviation
    DMet::Distrib::JensenShannon(res, v1, v2);

    double res1 = mpfr_get_d(res,GMP_RNDN);
    mpfr_clear(res);
//
    EXPECT_DOUBLE_EQ(0.3877869489141309, res1);
}