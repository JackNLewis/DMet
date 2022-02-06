//
// Created by jack lewis on 03/11/2021.
//

#include "gtest/gtest.h"
#include "DMet/PointDistances.h"
#include "DMet/EqWidthBin.h"
#include "DMet/DistribDistance.h"
#include <fstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>


using std::cout;
using std::endl;
using std::vector;
//================== KL DIV TESTS ===================//

//
//TEST(KLTest, BasicTest) {
//    vector<double> v1 {0.1,0.9};
//    vector<double> v2 {0.2,0.8};
//
//    DMet::Distrib::KLDiv(v1,v2);
//
//}
//
//TEST(Binning, BinningTests) {
//    DMet::EqWidthBin bin = DMet::EqWidthBin();
//    vector<vector<double>> vect
//            {
//                    {1, 5, 3},
//                    {4, 1, 6},
//                    {9, 6, 12},
//                    {0,0,0}
//            };
//    bin.setRanges(vect);
//    bin.generateBins(3);
//    bin.assignBins(vect);
//
//    cout << "[" << std::flush;
//    for(double d: bin.getPDF()){
//        cout << d << "," <<std::flush;
//    }
//    cout << "]"<< endl;
//
//    bin.clearBins();
//    vector<vector<double>> vect2
//            {
//                    {1, 5, 3},
//                    {4, 1, 6},
//                    {9, 6, 12},
//                    {0,4,0}
//            };
//    bin.assignBins(vect2);
//    cout << "[" << std::flush;
//    for(double d: bin.getPDF()){
//        cout << d << "," <<std::flush;
//    }
//    cout << "]" << endl;
//
//}

//
//TEST(GMPTest, BasicTest) {
//    double v1[] = {0.0,0.0};
//    double v2[] = {1.7976931348623157e+308, 1.7976931348623157e+308};
//
//    int size1 = sizeof(v1)/sizeof(v1[0]);
//    int size2 = sizeof(v2)/sizeof(v2[0]);
//
//    DMet::PointDistances::getMinkowski(NULL,v1,v2,size1,size2,2,8);
//    DMet::PointDistances::getMinkowski(NULL,v2,v1,size1,size2,2,8);
//}
//
////================== MINKOWSKI TESTS ===================//


////================== EUCLIDEAN TESTS ===================//
//TEST(EuclideanTests, WorkingPoints) {
//    double arr1[] = {1.0f, 0.0f, 0.0f};
//    double arr2[] = {0.0f,1.0f,0.0f};
//    int size1 = sizeof(arr1) / sizeof(arr1[0]);
//    int size2 = sizeof(arr2) / sizeof(arr2[0]);
//    EXPECT_FLOAT_EQ(1.4142135623730951, DMet::PointDistances::getEuclidean(arr1, arr2, size1, size2));
//}
//
//TEST(EuclideanTests, SamePoints) {
//    //test on same point
//    double arr1[] = {1.0f, 2.0f, 3.0f, 4.0f};
//    double arr2[] = {1.0f, 2.0f, 3.0f, 4.0f};
//    int size1 = sizeof(arr1) / sizeof(arr1[0]);
//    int size2 = sizeof(arr2) / sizeof(arr2[0]);
//    EXPECT_FLOAT_EQ(0, DMet::PointDistances::getEuclidean(arr1, arr2, size1, size2));
//}
//
//TEST(EuclideanTests, EpsilonDistance) {
//    //test with there being a difference of epsilon and producing a value that isn't zero
//    double arr1[] = {1.0,1.0,1.0};
//    double arr2[] = {1.0,1.0,1.0+std::numeric_limits<float>::epsilon()};
//    int size1 = sizeof(arr1) / sizeof(arr1[0]);
//    int size2 = sizeof(arr2) / sizeof(arr2[0]);
//    EXPECT_DOUBLE_EQ(std::numeric_limits<float>::epsilon(), DMet::PointDistances::getEuclidean(arr1, arr2, size1, size2));
//}
//
////================== Chebyshev TESTS ===================//
//TEST(ChebyshevTests, WorkingPoints) {
//    //test with there being a difference of epsilon and producing a value that isn't zero
//    double arr1[] = {1.0,1.0,1.0};
//    double arr2[] = {1.0,0.0,1.0};
//    int size1 = sizeof(arr1) / sizeof(arr1[0]);
//    int size2 = sizeof(arr2) / sizeof(arr2[0]);
//    EXPECT_DOUBLE_EQ(1.0, DMet::PointDistances::getChebyshev(arr1, arr2, size1, size2));
//
//    double arr3[] = {8.0,0.0,1.0};
//    double arr4[] = {0.0,0.0,1.0};
//    EXPECT_DOUBLE_EQ(8.0, DMet::PointDistances::getChebyshev(arr3, arr4, size1, size2));
//}
