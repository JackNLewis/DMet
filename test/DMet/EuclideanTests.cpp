//
// Created by jack lewis on 04/03/2022.
//

#include "gtest/gtest.h"
#include <fstream>
#include <vector>
#include <mpfr.h>
#include <DMet/PointDistances.h>
#include <iostream>
#include <chrono>

using std::cout;
using std::endl;
using std::vector;
using DMet::PointDistances::getEuclidean;
using std::ifstream;
using std::string;
using namespace std::chrono;

TEST(EuclideanTests, General) {
    ifstream data("/Users/jacklewis/Documents/work/year3/DMet/test/scripts/eucl_general");
    string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }


    while(std::getline(data,field)){
        // get next line in file
        while ( getline(data,line) ){
            std::stringstream ss(line);
            vector<double> v1;
            vector<double> v2;
            double ans;
            //split line
            int index = 0;
            while(std::getline(ss,field,',')){
                if(index < 4){ // place in first array
                    v1.push_back(std::stod(field));
                }else if(index < 8){ // place in second array
                    v2.push_back(std::stod(field));
                }
                else{ // result for pval
                    ans = std::stod(field);
                }
                index++;
            }
            mpfr_t res;
            mpfr_init(res);

            auto start = high_resolution_clock::now(); // start clock
            getEuclidean(res,v1,v2);
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);

            double resDouble = mpfr_get_d(res,GMP_RNDN);
            mpfr_clear(res);
            EXPECT_FLOAT_EQ(resDouble, ans);
        }
    }
}

/**
 * Test to library produces consistent results with the same points are input
 */
TEST(EuclideanTests, SamePoints){
    ifstream data("/Users/jacklewis/Documents/work/year3/DMet/test/scripts/eucl_same");
    string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }

    while(std::getline(data,field)) {
        // get next line in file
        while (getline(data, line)) {
            std::stringstream ss(line);
            vector<double> v1;
            vector<double> v2;
            double ans = 0.0;
            //split line
            int index = 0;
            while (std::getline(ss, field, ',')) {
                if (index < 4) { // place in first array
                    v1.push_back(std::stod(field));
                } else if (index < 8) { // place in second array
                    v2.push_back(std::stod(field));
                }
                index++;
            }
            mpfr_t res;
            mpfr_init(res);

            auto start = high_resolution_clock::now(); // start clock
            getEuclidean(res, v1, v2);
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);

            double resDouble = mpfr_get_d(res, GMP_RNDN);
            mpfr_clear(res);
            EXPECT_FLOAT_EQ(resDouble, ans);
        }
    }
}

/**
 * Test to see what happens when multiple infinite value are input. Tested on both input vectors
 * Tests for an input of infinite and negative infinity
 */
TEST(EuclideanTests, MachinePrecision){
    double eps = std::numeric_limits<double>::epsilon();
    vector<double> v1 {1.0,1.0,1.0};
    vector<double> v2 {1.0+eps,1.0,1.0};

    mpfr_t res;
    mpfr_init(res);
    getEuclidean(res, v1, v2);
//    mpfr_printf("Res %.5Re\n",res);
    double res_d = mpfr_get_d(res,GMP_RNDN);
    EXPECT_FLOAT_EQ(res_d,eps);
    mpfr_clear(res);
}

/**
 * Test to see if the method throws an error when incompatable sized arrays are input
 */
TEST(EuclideanTests, IncompatableSizes){
    vector<double> v1 {0.0,1.0,2.0};
    vector<double> v2 {0.0,1.0};
    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(getEuclidean(res, v1, v2),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test for whether v1 is empty or v2 is empty
 */
TEST(EuclideanTests, EmptyArguments){
    vector<double> v1;
    vector<double> v2;
    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(getEuclidean(res, v1, v2),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test to see what happens when a single infinite value is input. Tested on both input vectors
 */
TEST(EuclideanTests, SingleInfinite){
    double inf = std::numeric_limits<double>::infinity();

    vector<double> v1;
    vector<double> v2;
    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(getEuclidean(res, v1, v2),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test to see what happens when multiple infinite value are input. Tested on both input vectors
 * Tests for an input of infinite and negative infinity
 */
TEST(EuclideanTests, MultipleInfinite){

}
