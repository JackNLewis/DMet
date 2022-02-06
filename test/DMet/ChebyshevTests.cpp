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
    mpfr_t res;
    mpfr_init(res);
    DMet::PointDistances::getChebyshev(res,v1,v2);

    mpfr_clear(res);
//    EXPECT_DOUBLE_EQ(1.0, DMet::PointDistances::getChebyshev(arr1, arr2, size1, size2));

//    double arr3[] = {8.0,0.0,1.0};
//    double arr4[] = {0.0,0.0,1.0};
//    EXPECT_DOUBLE_EQ(8.0, DMet::PointDistances::getChebyshev(arr3, arr4, size1, size2));
}
