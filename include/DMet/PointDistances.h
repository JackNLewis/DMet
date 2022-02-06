//
// Created by jack lewis on 03/11/2021.
//

#pragma once
#include <gmp.h>
#include <mpfr.h>


/*!
 * Outer namespace of library
 */
namespace DMet{
    /*!
     * namespace for functions to compute point distances
     */
    namespace PointDistances{

        void getMinkowski(mpfr_t &res, std::vector<double> &vector1, std::vector<double> &vector2, double pvalue);
        void getManhattan(mpfr_t &res, std::vector<double> &vector1, std::vector<double> &vector2);
        void getEuclidean(mpfr_t &res, std::vector<double> &vector1, std::vector<double> &vector2);
        void getChebyshev(mpfr_t &res, std::vector<double> &vector1, std::vector<double> &vector2);
}}



