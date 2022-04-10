//
// Created by jack lewis on 22/11/2021.
//

#pragma once
#include <vector>
//

using std::vector;

namespace DMet{ namespace Utils{

    bool isPDF(std::vector<double> &v);


    /*!
        * Generates the cartesian product of two lists of vectors
        *
        * @param vec1
        * @param vec2
        * @return
        */
    vector<vector<vector<double>>> cartesian(vector<vector<vector<double>>> &vec1, vector<vector<double>> &vec2);

    }}
