//
// Created by jack lewis on 07/03/2022.
//

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
#include <math.h>

using std::cout;
using std::endl;
using std::vector;
using DMet::PointDistances::getManhattan;
using std::ifstream;
using std::string;
using namespace std::chrono;

TEST(ManhattanTests, General) {
    ifstream data("/Users/jacklewis/Documents/work/year3/DMet/test/scripts/man_general");
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
            getManhattan(res,v1,v2);
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
TEST(ManhattanTests, SamePoints){
    ifstream data("/Users/jacklewis/Documents/work/year3/DMet/test/scripts/same_points");
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
            getManhattan(res, v1, v2);
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);

            double resDouble = mpfr_get_d(res, GMP_RNDN);
            mpfr_clear(res);
            EXPECT_FLOAT_EQ(resDouble, ans);
        }
    }
}

/**
 * Test for whether v1 is empty or v2 is empty
 */
TEST(ManhattanTests, EmptyArguments){
    vector<double> v1;
    vector<double> v2;
    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(getManhattan(res, v1, v2),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test to see if the method throws an error when incompatable sized arrays are input
 */
TEST(ManhattanTests, IncompatableSizes){
    vector<double> v1 {0.0,1.0,2.0};
    vector<double> v2 {0.0,1.0};
    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(getManhattan(res, v1, v2),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test to see what happens when a single infinite value is input. Tested on both input vectors
 */
TEST(ManhattanTests, SingleInfinite){
    double inf = std::numeric_limits<double>::infinity();
    vector<double> v1 {inf,1.0,1.0};
    vector<double> v2 {2.0,1.0,3.0};
    mpfr_t res;
    mpfr_init(res);
    getManhattan(res, v1, v2);
    EXPECT_EQ(mpfr_get_d(res,GMP_RNDN), inf);
    mpfr_clear(res);


    vector<double> v3 {2.0,1.0,1.0};
    vector<double> v4 {inf,1.0,3.0};
    mpfr_t res2;
    mpfr_init(res2);
    getManhattan(res2, v3, v4);
    EXPECT_EQ(mpfr_get_d(res2,GMP_RNDN), inf);
    mpfr_clear(res2);
}

/**
 * Test to see what happens when multiple infinite value are input. Tested on both input vectors
 * Tests for an input of infinite and negative infinity
 */
TEST(ManhattanTests, MultipleInfinite){
    //need to see how to test for nan
    double inf = std::numeric_limits<double>::infinity();
    vector<double> v1 {inf,1.0,1.0};
    vector<double> v2 {inf,1.0,3.0};
    mpfr_t res;
    mpfr_init(res);
    getManhattan(res, v1, v2);
    double res_d = mpfr_get_d(res,GMP_RNDN);
    EXPECT_TRUE(isnan(res_d));
    mpfr_clear(res);

    vector<double> v3 {2.0,1.0,inf};
    vector<double> v4 {inf,1.0,3.0};
    mpfr_t res2;
    mpfr_init(res2);
    getManhattan(res2, v3, v4);
    EXPECT_EQ(mpfr_get_d(res2,GMP_RNDN), inf);
    mpfr_clear(res2);
}

/**
 * Tests for values larger than max value of 64 bit floating point
 */
TEST(ManhattanTests, Overflow){
    double max = std::numeric_limits<double>::max();
    vector<double> v1 {1.7e308,1.7e308};
    vector<double> v2 {0.0,0.0};
    string ans_string  = "3.4e308";
    mpfr_t res,ans;
    mpfr_inits(res,ans,NULL);
    mpfr_set_str(ans,ans_string.c_str(),10,GMP_RNDN);
    getManhattan(res, v1, v2);

    mpfr_printf("Result: %.5Re\n",res);
    EXPECT_TRUE(mpfr_cmp(res,ans) == 0);
    mpfr_clears(res,ans,NULL);
}

