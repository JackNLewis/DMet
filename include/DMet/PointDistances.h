//
// Created by jack lewis on 03/11/2021.
//

#pragma once
#include <gmp.h>
#include <mpfr.h>

using std::vector;

/*!
 * Outer namespace of library
 */
namespace DMet{
    /*!
     * namespace for functions to compute point distances
     */
    namespace PointDistances{

        void getMinkowski(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, double pvalue);
        void getMinkowskiInfP(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, bool pos);
        void getManhattan(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);
        void getEuclidean(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);
        void getChebyshev(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);
}}



