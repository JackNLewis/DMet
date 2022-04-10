//
// Created by jack lewis on 21/01/2022.
//

#include "DMet/EqWidthBin.h"
#include <iostream>
#include "DMet/utils.h"

using std::cout;
using std::endl;
using std::flush;

vector<vector<double>> DMet::EqWidthBin::setRanges(vector<vector<double>> &data) {
    int columns = data[0].size();
    // insert bounds -ininity and +infinity for all dimensions
    for(int i=0; i< columns;i++){
        vector<double> r{std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
        ranges.push_back(r);
    }

    // go through each row and update ranges of each dimensions
    for(vector<double> row : data){
        for(int i=0; i<columns;i++){
            //if smaller than min bound of that dimension
            if(row[i] == std::numeric_limits<double>::infinity() || row[i] == -std::numeric_limits<double>::infinity()){
                continue;
            }
            else if(row[i] < ranges[i][0]){
                ranges[i][0] = row[i];
            }else if(row[i] > ranges[i][1]){
                ranges[i][1] = row[i];
            }
        }
    }
    return data;
}

void DMet::EqWidthBin::generateBins(int arity) {
    // vectors where each dimension is binned
    vector<vector<vector<double>>> binnedDims;

    for(int i=0;i<ranges.size();i++){
        //get range
        vector<vector<double>> singleBinnedDim;
        double diff = ranges[i][1] - ranges[i][0];
        double width = diff/arity;
        for(int j=0;j<arity;j++){
            vector<double> t{ranges[i][0] + j*width, ranges[i][0] + (j+1)*width};
            if(j==0){
                t[0] = -std::numeric_limits<double>::infinity();
            }
            if(j == arity-1){
                t[1] = std::numeric_limits<double>::infinity();
            }
            singleBinnedDim.push_back(t);
        }
        binnedDims.push_back(singleBinnedDim);
    }

    // vectors which holds the cartesian product of all the binned dimensions
    vector<vector<vector<double>>> cartesianRanges;

    //starting from last dimension do cartesion onto a binnedDims
    for(int i= binnedDims.size() - 1; i >= 0; i--){
        DMet::Utils::cartesian(cartesianRanges, binnedDims[i]);
    }

    //assign the ranges the bin struct
    for(auto a : cartesianRanges){
        struct Bin bin;
        bin.range = a;
        bins.push_back(bin);
    }
}

/**
 * assigns data points to correct bins if it falls in the bins range.
 * It falls in the range if the point x is
 *      lower bound <= x < upper bound
 * if upper bound = inf then point is in range to make sure points with infinity are assigned
 *      lower bound <= x <= upper bound

 * @param data
 */
void DMet::EqWidthBin::assignBins(vector<vector<double>> &data){
    for(auto point: data){
        assignPoint(point);
    }
}

void DMet::EqWidthBin::assignPoint(vector<double> &point) {
    for(Bin &bin : bins){
        int dimensions = bin.range.size();
        bool addToBin = true;
        for(int i=0; i<dimensions;i++){
            vector<double> dimRange = bin.range[i];
            //if bin has upper bound infinity allow to be equal
            bool inRange;
            if(dimRange[1] == std::numeric_limits<double>::infinity()){
                inRange = (point[i] >= dimRange[0] && point[i] <= dimRange[1]);
            }else{
                inRange = (point[i] >= dimRange[0] && point[i] < dimRange[1]);
            }


            addToBin = addToBin && inRange;
        }
        if(addToBin){
            bin.values.push_back(point);
        }
    }
}


//std::ostream &DMet::operator<<(std::ostream &out, const DMet::EqWidthBin::Bin &bin) {
//    //print all the ranges of the bin
//    for(auto r : bin.range){
//        out << "[" << r[0] << "," << r[1] << "]" << std::flush;
//    }
//    out << ": " << bin.values.size();
//    return out;
//}

void DMet::EqWidthBin::clearBins() {
    for(auto &bin : bins){
        bin.values.clear();
    }
}

vector<double> DMet::EqWidthBin::getPDF() {
    vector<double> pdf;
    int count = 0;
    for(auto &bin : bins){
        pdf.push_back(bin.values.size());
        count+=bin.values.size();
    }

    for(double &d : pdf){
        d /= count;
    }

    return pdf;
}

//void DMet::EqWidthBin::printBins(const DMet::EqWidthBin &bin) {
//    for(auto b : bin.bins){
//        cout << b << endl;
//    }
//}


