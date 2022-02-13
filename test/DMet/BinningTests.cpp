//
// Created by jack lewis on 13/02/2022.
//
#include "gtest/gtest.h"
#include "DMet/PointDistances.h"
#include "DMet/EqWidthBin.h"
#include "DMet/DistribDistance.h"
#include <fstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>

using std::cout;
using std::endl;
using std::flush;
using std::vector;


TEST(Binning, BinningRanges){

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
    vector<double> t{-std::numeric_limits<double>::infinity(),4};
    bin.assignPoint(t);

    for(auto a : bin.bins){
        cout << a << endl;
    }
}



TEST(Binning, BinningTests) {
    DMet::EqWidthBin bin = DMet::EqWidthBin();
    vector<vector<double>> vect
            {
                    {1, 5, 3},
                    {4, 1, 6},
                    {9, 6, 12},
                    {0,0,0}
            };
    bin.setRanges(vect);
    bin.generateBins(3);
    bin.assignBins(vect);
//
//    cout << "[" << std::flush;
//    for(double d: bin.getPDF()){
//        cout << d << "," <<std::flush;
//    }
//    cout << "]"<< endl;

    bin.clearBins();
    vector<vector<double>> vect2
            {
                    {1, 5, 3},
                    {4, 1, 6},
                    {9, 6, 12},
                    {0,4,0}
            };
    bin.assignBins(vect2);
//    cout << "[" << std::flush;
//    for(double d: bin.getPDF()){
//        cout << d << "," <<std::flush;
//    }
//    cout << "]" << endl;

}