//
// Created by jack lewis on 07/03/2022.
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

TEST(EuclideanStrTests, General) {
    ifstream data("/Users/jacklewis/Documents/work/year3/DMet/test/scripts/eucl_general");
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