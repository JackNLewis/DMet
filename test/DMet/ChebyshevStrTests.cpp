//
// Created by jack lewis on 07/03/2022.
//

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
using DMet::PointDistances::getChebyshev;
using std::ifstream;
using std::string;
using namespace std::chrono;

/**
 * General tests for manhattan distance. Should all compute a value
 */
TEST(ChebyshevStrTests, General) {

    namespace fs = std::__fs::filesystem;
    fs::path p = std::__fs::filesystem::current_path();
    ifstream data(p.parent_path().parent_path().string() + "/test/scripts/cheb_general");

    string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }

    while(std::getline(data,field)){
        // get next line in file
        while ( getline(data,line) ){
            std::stringstream ss(line);
            vector<string> v1;
            vector<string> v2;
            double ans;
            //split line
            int index = 0;
            while(std::getline(ss,field,',')){
                if(index < 4){ // place in first array
                    v1.push_back(field);
                }else if(index < 8){ // place in second array
                    v2.push_back(field);
                }
                else{ // result for ans
                    ans = std::stod(field);
                }
                index++;
            }
            mpfr_t res;
            mpfr_init(res);

            auto start = high_resolution_clock::now(); // start clock
            getChebyshev(res,v1,v2);
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
TEST(ChebyshevStrTests, SamePoints){
    namespace fs = std::__fs::filesystem;
    fs::path p = std::__fs::filesystem::current_path();
    ifstream data(p.parent_path().parent_path().string() + "/test/scripts/same_points");
    string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }

    while(std::getline(data,field)) {
        // get next line in file
        while (getline(data, line)) {
            std::stringstream ss(line);
            vector<string> v1;
            vector<string> v2;
            double ans = 0.0;
            //split line
            int index = 0;
            while (std::getline(ss, field, ',')) {
                if (index < 4) { // place in first array
                    v1.push_back(field);
                } else if (index < 8) { // place in second array
                    v2.push_back(field);
                }
                index++;
            }
            mpfr_t res;
            mpfr_init(res);

            auto start = high_resolution_clock::now(); // start clock
            getChebyshev(res, v1, v2);
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
TEST(ChebyshevStrTests, EmptyArguments){
    vector<string> v1;
    vector<string> v2;
    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(getChebyshev(res, v1, v2),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test to see if the method throws an error when incompatable sized arrays are input
 */
TEST(ChebyshevStrTests, IncompatableSizes){
    vector<string> v1 {"1.0","2.0","3.0"};
    vector<string> v2 {"0.0","1.0"};
    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(getChebyshev(res, v1, v2),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test to see what happens when a single infinite value is input. Tested on both input vectors
 */
TEST(ChebyshevStrTests, SingleInfinite){
    double inf = std::numeric_limits<double>::infinity();
    vector<string> v1 {"inf","1.0","1.0"};
    vector<string> v2 {"2.0","1.0","3.0"};
    mpfr_t res;
    mpfr_init(res);
    getChebyshev(res, v1, v2);
    EXPECT_EQ(mpfr_get_d(res,GMP_RNDN), inf);
    mpfr_clear(res);


    vector<string> v3 {"2.0","1.0","1.0"};
    vector<string> v4 {"inf","1.0","3.0"};
    mpfr_t res2;
    mpfr_init(res2);
    getChebyshev(res2, v3, v4);
    EXPECT_EQ(mpfr_get_d(res2,GMP_RNDN), inf);
    mpfr_clear(res2);
}

/**
 * Test to see what happens when multiple infinite value are input. Tested on both input vectors
 * Tests for an input of infinite and negative infinity
 */
TEST(ChebyshevStrTests, MultipleInfinite){
    //need to see how to test for nan
    double inf = std::numeric_limits<double>::infinity();
    vector<string> v1 {"inf","1.0","1.0"};
    vector<string> v2 {"inf","1.0","3.0"};
    mpfr_t res;
    mpfr_init(res);
    getChebyshev(res, v1, v2);
    double res_d = mpfr_get_d(res,GMP_RNDN);
    EXPECT_TRUE(isnan(res_d));
    mpfr_clear(res);

    vector<string> v3 {"2.0","1.0","inf"};
    vector<string> v4 {"inf","1.0","3.0"};
    mpfr_t res2;
    mpfr_init(res2);
    getChebyshev(res2, v3, v4);
    EXPECT_EQ(mpfr_get_d(res2,GMP_RNDN), inf);
    mpfr_clear(res2);
}

/**
 * Tests for values smaller than precision of 64 bit floating point
 */
TEST(ChebyshevStrTests, Underflow){
    vector<string> v1 {"2.0e-400","1.0e-350"};
    vector<string> v2 {"0","0"};
    string ans_string  = "2.0e-350";
    mpfr_t res,ans;
    mpfr_inits(res,ans,NULL);
    mpfr_set_str(ans,ans_string.c_str(),10,GMP_RNDN);
    getChebyshev(res, v1, v2);

    mpfr_printf("Result: %.5Re\n",res);
    EXPECT_TRUE(mpfr_cmp(res,ans));
    mpfr_clears(res,ans,NULL);
}

TEST(ChebyshevStrTests, Overflow){
    double max = std::numeric_limits<double>::max();
    vector<string> v1 {"2.0e+400","1.0e+308"};
    vector<string> v2 {"0","0"};
    string ans_string  = "2.0e+400";
    mpfr_t res,ans;
    mpfr_inits(res,ans,NULL);
    mpfr_set_str(ans,ans_string.c_str(),10,GMP_RNDN);
    getChebyshev(res, v1, v2);

    mpfr_printf("Result: %.5Re\n",res);
    EXPECT_TRUE(mpfr_cmp(res,ans)==0);
    mpfr_clears(res,ans,NULL);
}

