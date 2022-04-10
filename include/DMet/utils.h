//
// Created by jack lewis on 22/11/2021.
//

#pragma once
#include <vector>
//

using std::vector;

namespace DMet{ namespace Utils{

    /***
     * Function to check is a vector is a valid pdf by summing to 1.0
     * @param v the input probability density function
     * @return true for a valid probability density function
     */
    bool isPDF(std::vector<double> &v);


    /*!
    * Generates the cartesian product of vec2 and every vector in vec1 uses vec1 as a working result.
    *
    * @param vec1 vector 1 used to as a working result
    * @param vec2 vector 2
    * @return the cartesian product of vector 1 and vector 2
    */
    vector<vector<vector<double>>> cartesian(vector<vector<vector<double>>> &vec1, vector<vector<double>> &vec2);

    }}
