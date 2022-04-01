//
// Created by jack lewis on 13/02/2022.
//
#include "gtest/gtest.h"
#include "DMet/EqWidthBin.h"
#include <vector>
#include <random>
#include <chrono>

using std::cout;
using std::endl;
using std::flush;
using std::vector;

using namespace std::chrono;

/**
 * Test that the assign bins function works
 */
TEST(Binning, BinningGeneral){
    DMet::EqWidthBin bin = DMet::EqWidthBin();
    vector<vector<double>> vect
            {
                    {1, 5},
                    {4, 1},
                    {9, 6},
                    {0,0}
            };
    bin.setRanges(vect);
    bin.generateBins(3);
    bin.assignBins(vect);
    //should create values in order [1,0,1,1,0,0,0,0,1]
    vector<double> values {1,0,1,1,0,0,0,0,1};

    for(int i=0;i<bin.bins.size();i++){
        int values_size = bin.bins[i].values.size();
        EXPECT_FLOAT_EQ(values_size,values[i]);
    }
    vector<double> t{-std::numeric_limits<double>::infinity(),4};
    bin.assignPoint(t);
}


/**
 * Test that correct bins are created
 */
TEST(Binning, BinCreation) {
    double inf = std::numeric_limits<double>::infinity();
    DMet::EqWidthBin bin = DMet::EqWidthBin();
    vector<vector<double>> vect
            {
                    {1, 5},
                    {4, 1},
                    {9, 6},
                    {0,0}
            };
    bin.setRanges(vect);
    bin.generateBins(3);

    vector<vector<vector<double>>> res{
            {{-inf,3},{-inf,2}},
            {{-inf,3},{2,4}},
            {{-inf,3},{4,inf}},
            {{3,6},{-inf,2}},
            {{3,6},{2,4}},
            {{3,6},{4,inf}},
            {{6,inf},{-inf,2}},
            {{6,inf},{2,4}},
            {{6,inf},{4,inf}},
    };

    //compare the ranges
    for(int i=0;i<bin.bins.size();i++){
        vector<vector<double>> ranges = bin.bins[i].range;
        for(int j=0;j<ranges.size();j++){ //go through each dimension
            for(int k=0;k<ranges.size();k++){//for each bound of dimension
                EXPECT_FLOAT_EQ(ranges[j][k],res[i][j][k]);
            }
        }
    }
}


/**
 * Test the bins work with infinite values
[-inf,3][-inf,2]
[-inf,3][2,4]
[-inf,3][4,inf]
[3,6][-inf,2]
[3,6][2,4]
[3,6][4,inf]
[6,inf][-inf,2]
[6,inf][2,4]
[6,inf][4,inf]
 Adding point (-inf,4) will be added to bin 3
 */
TEST(Binning, InfiniteVals) {
    DMet::EqWidthBin bin = DMet::EqWidthBin();
    vector<vector<double>> vect
            {
                    {1, 5},
                    {4, 1},
                    {9, 6},
                    {0,0}
            };
    vector<double> before{0.25,0,0.25,0.25,0,0,0,0,0.25};
    vector<double> after{0.2,0,0.4,0.2,0,0,0,0,0.2};
    bin.setRanges(vect);
    bin.generateBins(3);
    bin.assignBins(vect);
    bin.printBins(bin);
    vector<double> pdf_before (bin.getPDF());
    vector<double> t{-std::numeric_limits<double>::infinity(),4};// pdf should increase in bin {-inf-3,2-4}
    bin.assignPoint(t);
    cout << endl;vector<double> pdf_after (bin.getPDF());
    for(int i=0;i<pdf_before.size();i++){
        EXPECT_FLOAT_EQ(pdf_before[i],before[i]);
        EXPECT_FLOAT_EQ(pdf_after[i],after[i]);
    }
}

// /**
// * Set a time limit of 5 seconds,
// * see how many dimensions it gets up to
// */
TEST(Binning, MaxDims) {
    for(int z=0;z<11;z++){
        double lower_bound = 0;
        double upper_bound = 10000;
        std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
        std::default_random_engine re;

        int size = 8;
        vector<vector<double>> vect;
        //generate 100 data points with a single dimension
        for(int i=0;i<100;i++){
            vector<double> v;
            for(int j=0;j<size;j++){
                double random_d = unif(re);
                v.push_back(random_d);
            }
            vect.push_back(v);
        }

        //increase the dimension and calculate time of generating bins
        DMet::EqWidthBin bin = DMet::EqWidthBin();

        auto start = high_resolution_clock::now();
        bin.setRanges(vect);
        bin.generateBins(3);
        bin.assignBins(vect);
        auto end = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(end - start);
        cout << duration.count() <<", " <<flush;
        if(duration.count() > 2000000){
            cout << "Max Dim: "<< size;
            break;
        }
    }

}