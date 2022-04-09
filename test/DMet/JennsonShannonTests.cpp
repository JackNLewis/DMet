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
using std::ifstream;
using std::string;

TEST(JSTests, SimpleTest) {
    //test with there being a difference of epsilon and producing a value that isn't zero
    vector<double> v1 {0.5,0,0,0.25,0,0,0,0.25,0};
    vector<double> v2 {0.25,0,0,0.25,0,0.25,0,0.25,0};

    mpfr_t res;
    mpfr_init(res);

    DMet::Distrib::JensenShannon(res, v1, v2);

    double res1 = mpfr_get_d(res,GMP_RNDN);
    mpfr_clear(res);
    EXPECT_DOUBLE_EQ(0.10788077716941784, res1);
}

//
TEST(JSTests, WorkingPoints) {
    namespace fs = std::__fs::filesystem;
    fs::path p = std::__fs::filesystem::current_path();
    ifstream data(p.parent_path().parent_path().string() + "/test/scripts/js_general");

    string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }

    mpfr_t res;
    mpfr_init(res);
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
                if(index == 0){
                    ans = std::stod(field);
                }
                else if(index < 5+1){ // place in first array
                    v1.push_back(std::stod(field));
                }else if(index < 10+1){ // place in second array
                    v2.push_back(std::stod(field));
                }
                index++;
            }
            DMet::Distrib::JensenShannon(res,v1,v2);
            double resDouble = mpfr_get_d(res,GMP_RNDN);

            EXPECT_FLOAT_EQ(resDouble, ans);
        }
    }
    mpfr_clear(res);
}

//
TEST(JSTests, Symmetric) {
    namespace fs = std::__fs::filesystem;
    fs::path p = std::__fs::filesystem::current_path();
    ifstream data(p.parent_path().parent_path().string() + "/test/scripts/js_general");
    string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }

    mpfr_t res;
    mpfr_init(res);
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
                if(index == 0){
                    ans = std::stod(field);
                }
                else if(index < 5+1){ // place in first array
                    v1.push_back(std::stod(field));
                }else if(index < 10+1){ // place in second array
                    v2.push_back(std::stod(field));
                }
                index++;
            }

            DMet::Distrib::JensenShannon(res,v1,v2);
            double res1 = mpfr_get_d(res,GMP_RNDN);
            DMet::Distrib::JensenShannon(res,v2,v1);
            double res2 = mpfr_get_d(res,GMP_RNDN);

            EXPECT_FLOAT_EQ(res1,res2);
        }
    }
    mpfr_clear(res);

}

/**
 * Tests for when inputs in x are zero
 * Should have no effect on the sum since xlog(x/y) = 0 if x = 0
 * So if extra values of 0 are added it shouldn't change the distance
 */
TEST(JSTests, ZeroDiff){
    vector<double> v1 {0.9,0.1};
    vector<double> v2 {0.5,0.5};

    vector<double> v3 {0.9,0.1,0.0};
    vector<double> v4 {0.5,0.5,0.0};
    mpfr_t res;
    mpfr_init(res);

    DMet::Distrib::JensenShannon(res,v1,v2);
    double res1 = mpfr_get_d(res,GMP_RNDN);

    DMet::Distrib::JensenShannon(res,v3,v4);
    double res2 = mpfr_get_d(res,GMP_RNDN);

    EXPECT_FLOAT_EQ(res1,res2);
    EXPECT_FLOAT_EQ(res1,0.10174922507919676);

    mpfr_clear(res);
}

/**
 * Test for when x>0 but y=0 this will lead to xlog(x/0) = inf
 */
TEST(JSTests, Infinity){
    vector<double> v1 {0.8,0.1,0.1};
    vector<double> v2 {0.5,0.0,0.5};

    mpfr_t res;
    mpfr_init(res);

    DMet::Distrib::KLDiv(res,v1,v2);
    double res1 = mpfr_get_d(res,GMP_RNDN);

    EXPECT_FLOAT_EQ(res1,std::numeric_limits<double>::infinity());
    mpfr_clear(res);
}

/**
 * Test for when x>0 but y=0 this will lead to xlog(x/0) = inf
 */
TEST(JSTests, InvalidPdf){
    vector<double> v1 {0.1,0.1,0.1};
    vector<double> v2 {0.5,0.0,0.5};

    vector<double> v3 {2.0,3.4,7,2};
    vector<double> v4 {0.5,0.0,0.5};

    vector<double> v5 {0.1,0.8,0.1};
    vector<double> v6 {2.0,3.4,7,2};

    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(DMet::Distrib::JensenShannon(res,v1,v2),std::invalid_argument);
    EXPECT_THROW(DMet::Distrib::JensenShannon(res,v3,v4),std::invalid_argument);
    EXPECT_THROW(DMet::Distrib::JensenShannon(res,v5,v6),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test for when arrays are different lengths
 */
TEST(JSTests, IncompatableSizes){
    vector<double> v1 {0.1,0.9};
    vector<double> v2 {0.5,0.0,0.5};

    vector<double> v3 {0.1,0.8,0.1};
    vector<double> v4 {0.5,0.5};

    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(DMet::Distrib::JensenShannon(res,v1,v2),std::invalid_argument);
    EXPECT_THROW(DMet::Distrib::JensenShannon(res,v3,v4),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * A simple test to make sure the KL works with the binning class to produce correct output
 */
TEST(JSTests, Binning){
    vector<vector<double>> v1
            {
                    {1, 1},
                    {4, 1},
                    {9, 6},
                    {0,0}
            };

    vector<vector<double>> v2
            {
                    {0, 2},
                    {5, 3},
                    {3, 18},
                    {9,9}
            };

    /*
     * Equal width binning across both dimensions
     * Dimension 1 - 0, 9
     * Dimension 2 - 0, 18
     *
     * Using 3 bins
     * Dim 1
     * (0,3)(3,6)(6,9)
     * Dim 2
     * (0,6)(6,12)(12,18)
     * Cartesian product
     * ((0,3)(0,6)) - 2, 1
     * ((0,3)(6,12)) - 0, 0
     * ((0,3)(12,18)) - 0, 0
     * ((3,6)(0,6)) - 1, 0
     * ((3,6)(6,12)) - 0, 0
     * ((3,6)(12,18)) - 0, 1
     * ((6,9)(0,6)) - 0, 0
     * ((6,9)(6,12)) - 1, 1
     * ((6,9)(12,18)) - 0, 1
     *
     * 2,0,0,1,0,0,0,1,0
     * 1,0,0,0,0,1,0,1,1
     */
    mpfr_t res;
    mpfr_init(res);
    DMet::Distrib::JensenShannon(res,v1,v2,3);
    EXPECT_FLOAT_EQ(0.10788077716941784,mpfr_get_d(res,GMP_RNDN));
//    mpfr_printf("Result: %.5Re\n",res);
    mpfr_clear(res);
}
