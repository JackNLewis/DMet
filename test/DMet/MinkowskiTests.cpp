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
using DMet::PointDistances::getMinkowski;


TEST(MinkowskiTests, WorkingInRange) {
    namespace fs = std::__fs::filesystem;
    fs::path p = std::__fs::filesystem::current_path();
    std::ifstream data(p.parent_path().parent_path().string() + "/test/scripts/minkowskiTests.csv");

    std::string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }
    // get next line in file
    while ( getline(data,line) ){
        std::stringstream ss(line);
        //split line
        vector<double> v1;
        vector<double> v2;

        int index = 0;
        double pval5;
        double pval10;
        while(std::getline(ss,field,',')){
            if(index <= 3){ // place in first array
                v1.push_back(std::stod(field));
            }else if(index <= 7){ // place in second array
                v2.push_back(std::stod(field));
            }else if(index == 8){ // result for p val 5
                pval5 = std::stod(field);
            }else{ // reslut for pval 10
                pval10 = std::stod(field);
            }
            index++;
        }
        mpfr_t res;
        mpfr_init(res);
        getMinkowski(res,v1,v2,5);
        double p5 = mpfr_get_d(res,GMP_RNDN);
        getMinkowski(res,v1,v2,10);
        double p10 = mpfr_get_d(res,GMP_RNDN);
        mpfr_clear(res);
        EXPECT_FLOAT_EQ(pval5, p5);
        EXPECT_FLOAT_EQ(pval10, p10);
    }

}

TEST(MinkowskiTests, SamePoints) {
    double pVals[] = {1.0,10.0,100.0};
    vector<double> v1 {1.0, 0.0, 0.0};
    vector<double> v2 {1.0,0.0,0.0};
    for(int i=0;i<1;i++){
        mpfr_t res;
        mpfr_init(res);
        getMinkowski(res,v1,v2,pVals[i]);
        double resCheck = mpfr_get_d(res,GMP_RNDN);

        EXPECT_FLOAT_EQ(0,resCheck);
        mpfr_clear(res);
    }

    vector<double> v3 {1.50625103e+308,1.34702226e+308,1.38046575e+307};
    for(int i=0;i<1;i++){
        mpfr_t res;
        mpfr_init(res);
        getMinkowski(res,v3,v3,pVals[i]);
        double resCheck = mpfr_get_d(res,GMP_RNDN);
        EXPECT_FLOAT_EQ(0,resCheck);
        mpfr_clear(res);
    }
}

TEST(MinkowskiTests, IncompatableSizes) {

    vector<double> v1 {1.0, 0.0, 0.0};
    vector<double> v2 {1.0,0.0};

    vector<double> v3 {1.50625103e+308};
    mpfr_t res;
    mpfr_init(res);

    EXPECT_THROW(getMinkowski(res,v1,v2,2),std::invalid_argument);
    EXPECT_THROW(getMinkowski(res,v3,v1,2),std::invalid_argument);
    mpfr_clear(res);
}

TEST(MinkowskiTests, SingleInfinite) {
    vector<double> v1 {1.0, std::numeric_limits<double>::infinity(),3.0};
    vector<double> v2 {1.0,0,0};

    mpfr_t res;
    mpfr_init(res);

    getMinkowski(res,v1,v2,2);
    double resCheck1 = mpfr_get_d(res,GMP_RNDN);

    getMinkowski(res,v2,v1,2);
    double resCheck2 = mpfr_get_d(res,GMP_RNDN);
    EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(),resCheck1);
    EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(),resCheck2);
    mpfr_clear(res);
}

TEST(MinkowskiTests, InfintePval) {
    vector<double> v1 {1.0, 1.0,50.0};
    vector<double> v2 {1.0,0,0};

    mpfr_t res;
    mpfr_init(res);

    getMinkowski(res,v1,v2,std::numeric_limits<double>::infinity());
    double resCheck1 = mpfr_get_d(res,GMP_RNDN);

    EXPECT_FLOAT_EQ(50,resCheck1);
    mpfr_clear(res);
}

TEST(MinkowskiTests, NegInfintePval) {
    vector<double> v1 {1.0, 1.0,50.0};
    vector<double> v2 {1.0,0,0};
    double inf = std::numeric_limits<double>::infinity();
    mpfr_t res;
    mpfr_init(res);

    getMinkowski(res,v1,v2,-inf);
    double resCheck1 = mpfr_get_d(res,GMP_RNDN);

    EXPECT_FLOAT_EQ(0,resCheck1);
    mpfr_clear(res);
}

TEST(MinkowskiTests, VariableLength){
    namespace fs = std::__fs::filesystem;
    fs::path p = std::__fs::filesystem::current_path();
    std::ifstream data(p.parent_path().parent_path().string() + "/test/scripts/mink_var_length");

    std::string line,field;

    if (!data.is_open()) {
        cout << "file not open" << endl;
    }
    // get next line in file
    mpfr_t res;
    mpfr_init(res);
    int inputSizes[] {10,100,1000};
    int index = 0;
    while (getline(data,line)){
        std::stringstream ss(line);
        //split line
        vector<double> v1;
        vector<double> v2;
        double ans;
        int i=0;

        while(getline(ss,field,',') && i<inputSizes[index]*2+1){
            if(i==0){
                ans = std::stod(field);
            }else if(i < inputSizes[index]+1){
                v1.push_back(std::stod(field));
            }else if(i < inputSizes[index]*2+1){
                v2.push_back(std::stod(field));
            }
            i++;
        }
        getMinkowski(res,v1,v2,10);
        double res_d = mpfr_get_d(res,GMP_RNDN);
        EXPECT_FLOAT_EQ(res_d,ans);
        index++;
    }
    mpfr_clear(res);
}

