//
// Created by jack lewis on 21/01/2022.
//

#pragma once
#include <fstream>
#include <vector>

using std::vector;

class EqWidthBin {
    struct Bin{
        int index;
        vector<vector<double>> range;
        vector<vector<double>> values;
    };

public:
    vector<vector<double>> ranges;
    vector<vector<vector<double>>> bins;

    EqWidthBin();
    vector<vector<double>> setRanges(vector<vector<double>>& data);
    void generateBins(int arity);
    vector<vector<vector<double>>> cartesian(vector<vector<vector<double>>>& vec1, vector<vector<double>>& vec2);

};


