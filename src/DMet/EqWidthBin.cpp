//
// Created by jack lewis on 21/01/2022.
//

#include "DMet/EqWidthBin.h"
#include <iostream>

using std::cout;
using std::endl;
using std::flush;

DMet::EqWidthBin::EqWidthBin() {}



vector<vector<double>> DMet::EqWidthBin::setRanges(vector<vector<double>>& data) {
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

//    for(vector<double> row : ranges){
//        cout << "Min: " << row[0] << " Max: " << row[1] <<  endl;
//    }
    return data;
}

void DMet::EqWidthBin::generateBins(int arity) {
    vector<vector<vector<double>>> res;
    //start from last columns create and array and do cartesian of that with res;
    for(int i=0;i<ranges.size();i++){
        //get range
        vector<vector<double>> binRanges;
        double diff = ranges[i][1] - ranges[i][0];
        double width = diff/arity;
        for(int j=0;j<arity;j++){
            vector<double> t{ranges[i][0] + j*width, ranges[i][0] + (j+1)*width};
            binRanges.push_back(t);
        }
        res.push_back(binRanges);
    }

    for(int i=res.size() ; i>0;i--){
        cout << "Dim: " << i << ": "<< flush;
        for(auto dimR: res[i-1]){
            cout <<  "[" << dimR[0] <<"," << dimR[1] << "] "<< flush;
        }
        cout << endl;
    }

    vector<vector<vector<double>>> binRanges;
    //starting from last dimension do cartesion onto a res
    for(int i=res.size()-1;i>=0;i--){
        cartesian(binRanges,res[i]);
    }

    for(auto a : binRanges){
        //for each vector of ranges
        for(auto dims : a){
            cout << "[" << dims[0] <<"," << dims[1] << "]" << std::flush;
        }
        cout<<endl;

        //assign the bins the range
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
