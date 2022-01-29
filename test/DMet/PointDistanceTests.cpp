//
// Created by jack lewis on 03/11/2021.
//

#include "gtest/gtest.h"
#include "DMet/PointDistances.h"
#include <fstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>

//using namespace DMet::Utils;
using std::cout;
using std::endl;


//================== KL DIV TESTS ===================//
//
//TEST(Binning, BinningTests) {
//    EqWidthBin bin = EqWidthBin();
//    vector<vector<double>> vect
//            {
//                    {1, 5, 3,4},
//                    {4, 1, 6,5},
//                    {9, 6, 12,8},
//                    {0,0,0,0}
//            };
//    bin.setRanges(vect);
//    bin.generateBins(3);
//}


TEST(GMPTest, BasicTest) {
    double v1[] = {0.0,0.0};
    double v2[] = {1.7976931348623157e+308, 1.7976931348623157e+308};

    int size1 = sizeof(v1)/sizeof(v1[0]);
    int size2 = sizeof(v2)/sizeof(v2[0]);
    DMet::getMinkowski(NULL,v1,v2,size1,size2,2,8);
    DMet::getMinkowski(NULL,v2,v1,size1,size2,2,8);
}

//================== MINKOWSKI TESTS ===================//
TEST(MinkowskiTests, GeneralTests) {
    std::ifstream data("/Users/jacklewis/CLionProjects/DistanceMetrics/DistMets_tests/tests.csv");
    std::string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }
    // get next line in file
    while ( getline(data,line) ){
        std::stringstream ss(line);
        //split line
        double arr1[4];
        double arr2[4];
        int index = 0;
        double pval5;
        double pval10;
        while(std::getline(ss,field,',')){
            if(index <= 3){ // place in first array
                arr1[index] = std::stod(field);
            }else if(index <= 7){ // place in second array
                arr2[index-4] = std::stod(field);
            }else if(index == 8){ // result for p val 5
                pval5 = std::stod(field);
            }else{ // reslut for pval 10
                pval10 = std::stod(field);
            }
            index++;
        }
        int size1 = sizeof(arr1)/sizeof(arr1[0]);
        int size2 = sizeof(arr2)/sizeof(arr2[0]);


        EXPECT_FLOAT_EQ(pval5,DMet::getMinkowski(arr1,arr2,size1,size2,5));
        EXPECT_FLOAT_EQ(pval10,DMet::getMinkowski(arr1,arr2,size1,size2,10));
    }

}

TEST(MinkowskiTests, WorkingPoints) {
    double pVals[] = {1.0f,3.0f};
    double arr1[] = {1.0f, 0.0f, 0.0f};
    double arr2[] = {0.0f,1.0f,0.0f};
    double res[] = {2.0f,1.2599210498948732};
    int size1 = sizeof(arr1) / sizeof(arr1[0]);
    int size2 = sizeof(arr2) / sizeof(arr2[0]);
    for(int i=0;i<sizeof(res) / sizeof(res[0]);i++){
        EXPECT_FLOAT_EQ(res[i], DMet::getMinkowski(arr1, arr2, size1, size2,pVals[i]));
    }
}

TEST(MinkowskiTests, Oveflow) {
    double maxVal = std::numeric_limits<double>::max();
    double arr1[] = {maxVal, maxVal, maxVal};
    double arr2[] = {1.0f,0.0f,0.0f};
    int size1 = sizeof(arr1) / sizeof(arr1[0]);
    int size2 = sizeof(arr2) / sizeof(arr2[0]);
    EXPECT_FLOAT_EQ(HUGE_VAL,DMet::getMinkowski(arr1,arr2,size1,size2,2));

}

TEST(MinkowskiTests, SamePoints) {
    double pVals[] = {1.0f,10.0f,100.0f};
    double arr1[] = {1.0f, 0.0f, 0.0f};
    double arr2[] = {1.0f,0.0f,0.0f};
    int size1 = sizeof(arr1) / sizeof(arr1[0]);
    int size2 = sizeof(arr2) / sizeof(arr2[0]);
    for(int i=0;i<1;i++){
        EXPECT_FLOAT_EQ(0,DMet::getMinkowski(arr1,arr2,size1,size2,pVals[i]));
    }

    double arr3[] = {1.50625103e+308,1.34702226e+308,1.38046575e+307};
    for(int i=0;i<1;i++){
        EXPECT_FLOAT_EQ(0,DMet::getMinkowski(arr3,arr3,size1,size2,pVals[i]));
    }
}

//================== EUCLIDEAN TESTS ===================//
TEST(EuclideanTests, WorkingPoints) {
    double arr1[] = {1.0f, 0.0f, 0.0f};
    double arr2[] = {0.0f,1.0f,0.0f};
    int size1 = sizeof(arr1) / sizeof(arr1[0]);
    int size2 = sizeof(arr2) / sizeof(arr2[0]);
    EXPECT_FLOAT_EQ(1.4142135623730951, DMet::getEuclidean(arr1, arr2, size1, size2));
}

TEST(EuclideanTests, SamePoints) {
    //test on same point
    double arr1[] = {1.0f, 2.0f, 3.0f, 4.0f};
    double arr2[] = {1.0f, 2.0f, 3.0f, 4.0f};
    int size1 = sizeof(arr1) / sizeof(arr1[0]);
    int size2 = sizeof(arr2) / sizeof(arr2[0]);
    EXPECT_FLOAT_EQ(0, DMet::getEuclidean(arr1, arr2, size1, size2));
}

TEST(EuclideanTests, EpsilonDistance) {
    //test with there being a difference of epsilon and producing a value that isn't zero
    double arr1[] = {1.0,1.0,1.0};
    double arr2[] = {1.0,1.0,1.0+std::numeric_limits<float>::epsilon()};
    int size1 = sizeof(arr1) / sizeof(arr1[0]);
    int size2 = sizeof(arr2) / sizeof(arr2[0]);
    EXPECT_DOUBLE_EQ(std::numeric_limits<float>::epsilon(), DMet::getEuclidean(arr1, arr2, size1, size2));
}

//================== Chebyshev TESTS ===================//
TEST(ChebyshevTests, WorkingPoints) {
    //test with there being a difference of epsilon and producing a value that isn't zero
    double arr1[] = {1.0,1.0,1.0};
    double arr2[] = {1.0,0.0,1.0};
    int size1 = sizeof(arr1) / sizeof(arr1[0]);
    int size2 = sizeof(arr2) / sizeof(arr2[0]);
    EXPECT_DOUBLE_EQ(1.0, DMet::getChebyshev(arr1, arr2, size1, size2));

    double arr3[] = {8.0,0.0,1.0};
    double arr4[] = {0.0,0.0,1.0};
    EXPECT_DOUBLE_EQ(8.0, DMet::getChebyshev(arr3, arr4, size1, size2));
}
