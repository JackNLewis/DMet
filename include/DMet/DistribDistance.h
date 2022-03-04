//
// Created by jack lewis on 22/11/2021.
//

#pragma once
#include <vector>

namespace DMet { namespace Distrib{
    using std::vector;

    void KLDiv(mpfr_t &res, vector<double> &v1, vector<double> &v2, bool pdfCheck);
    void KLDiv(mpfr_t &res, vector<double> &v1, vector<double> &v2);
    void KLDiv(mpfr_t &res, vector<vector<double>> &v1, vector<vector<double>> &v2); // allows binning


    void JensenShannon(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);
    void JensenShannon(mpfr_t &res, vector<vector<double>> &v1, vector<vector<double>> &v2);
}}

