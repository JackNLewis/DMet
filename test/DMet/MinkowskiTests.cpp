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



TEST(MinkowskiTests, WorkingInRange) {
    std::ifstream data("/Users/jacklewis/Documents/work/year3/DMet/test/minkowskiTests.csv");
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
        DMet::PointDistances::getMinkowski(res,v1,v2,5);
        double p5 = mpfr_get_d(res,GMP_RNDN);
        DMet::PointDistances::getMinkowski(res,v1,v2,10);
        double p10 = mpfr_get_d(res,GMP_RNDN);
        mpfr_clear(res);
        EXPECT_FLOAT_EQ(pval5, p5);
        EXPECT_FLOAT_EQ(pval10, p10);
    }

}

TEST(MinkowskiTests, SamePoints) {
    double pVals[] = {1.0f,10.0f,100.0f};
    vector<double> v1 {1.0f, 0.0f, 0.0f};
    vector<double> v2 {1.0f,0.0f,0.0f};
    for(int i=0;i<1;i++){
        mpfr_t res;
        mpfr_init(res);
        DMet::PointDistances::getMinkowski(res,v1,v2,pVals[i]);
        double resCheck = mpfr_get_d(res,GMP_RNDN);

        EXPECT_FLOAT_EQ(0,resCheck);
        mpfr_clear(res);
    }

    vector<double> v3 {1.50625103e+308,1.34702226e+308,1.38046575e+307};
    for(int i=0;i<1;i++){
        mpfr_t res;
        mpfr_init(res);
        DMet::PointDistances::getMinkowski(res,v3,v3,pVals[i]);
        double resCheck = mpfr_get_d(res,GMP_RNDN);
        EXPECT_FLOAT_EQ(0,resCheck);
        mpfr_clear(res);
    }
}

//TEST(MinkowskiTests, WorkingPoints) {
//    double pVals[] = {1.0f,2.0f};
//    vector<double> v1 {1.0f, 0.0f, 0.0f};
//    vector<double> v2 {0.0f,2.0f,0.0f};
//    for(int i=0;i<sizeof(pVals) / sizeof(pVals[0]);i++){
//        DMet::PointDistances::getMinkowski(NULL,v1, v2, pVals[i],200);
//    }
//}


//TEST(MinkowskiTests, Oveflow) {
//    double maxVal = std::numeric_limits<double>::max();
//    double arr1[] = {maxVal, maxVal, maxVal};
//    double arr2[] = {1.0f,0.0f,0.0f};
//    int size1 = sizeof(arr1) / sizeof(arr1[0]);
//    int size2 = sizeof(arr2) / sizeof(arr2[0]);
//    EXPECT_FLOAT_EQ(HUGE_VAL,DMet::PointDistances::getMinkowski(arr1,arr2,size1,size2,2));
//
//}
//

