//
// Created by jack lewis on 24/03/2022.
//
#include "gtest/gtest.h"
#include <fstream>
#include <vector>
#include <mpfr.h>
#include <DMet/PointDistances.h>
#include <DMet/EqWidthBin.h>


using std::cout;
using std::endl;
using std::vector;
using DMet::PointDistances::getMinkowski;
using DMet::PointDistances::getEuclidean;

TEST(Demo, InfintePval) {
    vector<double> v1 {1.0, 1.0,50.0};
    vector<double> v2 {1.0,0,0};

    mpfr_t res;
    mpfr_init(res);

    getMinkowski(res,v1,v2,std::numeric_limits<double>::infinity());
    double resCheck1 = mpfr_get_d(res,GMP_RNDN);

    EXPECT_FLOAT_EQ(50,resCheck1);
    mpfr_clear(res);
}

TEST(Demo, Overflow) {
    double max = std::numeric_limits<double>::max();
    vector<double> v1{1.0e+308, 1.0e+308, 1.0e+308, 1.0e+308};
    vector<double> v2{0, 0, 0, 0};
    mpfr_t res;
    mpfr_inits(res, NULL);
    getEuclidean(res, v1, v2);
//    mpfr_printf("Overflow Result %.5Re\n", res); //result "2.0e308"
}


/**
[-inf,3][-inf,2]
[-inf,3][2,4]
[-inf,3][4,inf]
[3,6][-inf,2]
[3,6][2,4]
[3,6][4,inf]
[6,inf][-inf,2]
[6,inf][2,4]
[6,inf][4,inf]
 */
TEST(Demo, Binning){
    vector<vector<double>> vect
            {
                    {1, 5},
                    {4, 1},
                    {9, 6},
                    {0,0}
            };

    DMet::EqWidthBin bin = DMet::EqWidthBin();
    bin.setRanges(vect);
    bin.generateBins(3);
    bin.assignBins(vect);
//    bin.printBins(bin);
    bin.getPDF();
}