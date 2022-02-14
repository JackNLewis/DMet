//
// Created by jack lewis on 22/11/2021.
//

#pragma once
#include <vector>

namespace DMet { namespace Distrib{
    using std::vector;

    void KLDiv(mpfr_t &res, std::vector<double> &v1, std::vector<double> &v2);
    void JennsonShannon(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);

}}

