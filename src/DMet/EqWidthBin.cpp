//
// Created by jack lewis on 21/01/2022.
//

#include "DMet/EqWidthBin.h"
#include <iostream>

using std::cout;
using std::endl;
using std::flush;



vector<vector<double>> DMet::EqWidthBin::setRanges(vector<vector<double>> &data) {
    int columns = data[0].size();
    int rows = data.size();
    // insert bounds -ininity and +infinity for all dimensions
    for(int i=0; i< columns;i++){
        vector<double> r{std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
        ranges.push_back(r);
    }

    // go through each row and update ranges of each dimensions
    for(vector<double> row : data){
        for(int i=0; i<columns;i++){
            //if smaller than min bound of that dimension
            if(row[i] < ranges[i][0]){
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
            singleBinnedDim.push_back(t);
        }
        binnedDims.push_back(singleBinnedDim);
    }


    // vectors which holds the cartesian product of all the binned dimensions
    vector<vector<vector<double>>> cartesianRanges;

    //starting from last dimension do cartesion onto a binnedDims
    for(int i= binnedDims.size() - 1; i >= 0; i--){
        cartesian(cartesianRanges, binnedDims[i]);
    }

    //assign the ranges the bin struct
    for(auto a : cartesianRanges){
        struct Bin bin;
        bin.range = a;
        bins.push_back(bin);
    }
}

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
            bool inRange = (point[i] >= dimRange[0] && point[i] < dimRange[1]);
            addToBin = addToBin && inRange;
        }
        if(addToBin){
            bin.values.push_back(point);
        }
    }
}

vector<vector<vector<double>>> DMet::EqWidthBin::cartesian(vector<vector<vector<double>>>& vec1, vector<vector<double>>& vec2){
    //no dimesnions yet
    if(vec1.size() == 0){
        for(auto v2 : vec2){
            vector<vector<double>> t{v2};
            vec1.push_back(t);
        }
    }
    //else multiply vec1 vec2.size times to get every combination
    else{
        vector<vector<vector<double>>> copy(vec1);
        int n = copy.size();
        for(int i=0;i<vec2.size()-1;i++){
            vec1.insert(vec1.end(),copy.begin(),copy.end());
        }

        int startPos = 0;
        for(auto d : vec2){
            //copy it n times into newly created vector res
            for(int i=0;i<n;i++){
                vec1[startPos+i].insert(vec1[startPos+i].begin(), d);
            }
            startPos+=n;
        }
    }
    return vec1;
}

std::ostream &DMet::operator<<(std::ostream &out, const DMet::EqWidthBin::Bin &bin) {
    //print all the ranges of the bin

    for(auto r : bin.range){
        out << "[" << r[0] << "," << r[1] << "]" << std::flush;
    }
    out << ": " << bin.values.size() << endl;
    return out;
}

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
