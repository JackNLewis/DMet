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
using DMet::Distrib::KLDiv;
using std::ifstream;
using std::string;

TEST(KLTests, WorkingPoints) {
    ifstream data("/Users/jacklewis/Documents/work/year3/DMet/test/scripts/kl_general");
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

            DMet::Distrib::KLDiv(res,v1,v2);
            double resDouble = mpfr_get_d(res,GMP_RNDN);
            EXPECT_FLOAT_EQ(resDouble, ans);
        }

    }
    mpfr_clear(res);
}

TEST(KLTests, NonSymetric) {
   vector<double> v1 {0.1,0.2,0.7};
   vector<double> v2 {0.3,0.4,0.3};

   mpfr_t res;
   mpfr_init(res);

   KLDiv(res,v1,v2);
   double res1 = mpfr_get_d(res,GMP_RNDN);
   KLDiv(res,v2,v1);
   double res2 = mpfr_get_d(res,GMP_RNDN);
   //distance x to y
   EXPECT_FLOAT_EQ(0.3446178372922425,res1);
   //distacen y to x
   EXPECT_FLOAT_EQ(0.3526532007082499,res2);
   mpfr_clear(res);
}

/**
 * Tests for when inputs in x are zero
 * Should have no effect on the sum since xlog(x/y) = 0 if x = 0
 * So if extra values of 0 are added it shouldn't change the distance
 */
TEST(KLTests, ZeroDiff){
    vector<double> v1 {0.9,0.1};
    vector<double> v2 {0.5,0.5};

    vector<double> v3 {0.9,0.1,0.0};
    vector<double> v4 {0.5,0.5,0.0};
    mpfr_t res;
    mpfr_init(res);

    DMet::Distrib::KLDiv(res,v1,v2);
    double res1 = mpfr_get_d(res,GMP_RNDN);

    DMet::Distrib::KLDiv(res,v3,v4);
    double res2 = mpfr_get_d(res,GMP_RNDN);

    cout << "res 1: " << res1 << endl;
    cout << "res 2: " << res2 << endl;
    EXPECT_FLOAT_EQ(res1,res2);
    EXPECT_FLOAT_EQ(res1,0.3680642071684971);
    mpfr_clear(res);
}

/**
 * Test for when x>0 but y=0 this will lead to xlog(x/0) = inf
 */
TEST(KLTests, Infinity){
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
TEST(KLTests, InvalidPdf){
    vector<double> v1 {0.1,0.1,0.1};
    vector<double> v2 {0.5,0.0,0.5};

    vector<double> v3 {2.0,3.4,7,2};
    vector<double> v4 {0.5,0.0,0.5};

    vector<double> v5 {0.1,0.8,0.1};
    vector<double> v6 {2.0,3.4,7,2};

    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(DMet::Distrib::KLDiv(res,v1,v2),std::invalid_argument);
    EXPECT_THROW(DMet::Distrib::KLDiv(res,v3,v4),std::invalid_argument);
    EXPECT_THROW(DMet::Distrib::KLDiv(res,v5,v6),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * Test for when x>0 but y=0 this will lead to xlog(x/0) = inf
 */
TEST(KLTests, IncompatableSizes){
    vector<double> v1 {0.1,0.9};
    vector<double> v2 {0.5,0.0,0.5};

    vector<double> v3 {0.1,0.8,0.1};
    vector<double> v4 {0.5,0.5};

    mpfr_t res;
    mpfr_init(res);
    EXPECT_THROW(DMet::Distrib::KLDiv(res,v1,v2),std::invalid_argument);
    EXPECT_THROW(DMet::Distrib::KLDiv(res,v3,v4),std::invalid_argument);
    mpfr_clear(res);
}

/**
 * A simple test to make sure the KL works with the binning class to produce correct output
 */
TEST(KLTests, Binning){
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
    DMet::Distrib::KLDiv(res,v1,v2,3);
    mpfr_printf("Result: %.5Re\n",res);
    mpfr_clear(res);
}
