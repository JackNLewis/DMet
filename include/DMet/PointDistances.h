//
// Created by jack lewis on 03/11/2021.
//

#pragma once
#include <gmp.h>
#include <mpfr.h>

using std::vector;
using std::string;
/*!
 * Outer namespace of library
 */
namespace DMet{
    /*!
     * namespace for functions to compute point distances
     */
    namespace PointDistances{

        // Distances accepting double as inputs
        /*!
        * Function that computes the distance between vector 1 and vector 2
        * Uses Minkowski distance
        * Uses precision 200
        * (∑|Xi−Yi|^p)^1/p
        *
        * @param vector1
        * @param vector2
        * @param size1 size of vector1
        * @param size2 size of vector2
        * @param pvalue pvalue
        * @return floating point result
        */
        void getMinkowski(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, double pvalue);

        /**
         * Calculates minkowski distance when the p value is either +inf or -inf. This is equivalent to finding the
         * largest distance when p is +inf or finding the smallest distance when p is -inf.
         *
         * @param res
         * @param vector1
         * @param vector2
         * @param pos
         */
        void getMinkowskiInfP(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, bool pos);

        /*!
         * Function that computes the Manhattan distance.
         * Uses the implementation of Minkowski distance with a p value set as 1
         * @param vector1
         * @param vector2
         * @param size1 size of vector1
         * @param size2 size of vector2
         * @return floating point result
        */
        void getManhattan(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);

        /*!
         * Function that computes the Euclidean distance.
         * Uses the implementation of Minkowski distance with a p value set as 2
         *
         * @param vector1 input vector 1
         * @param vector2 input  vector 2
         * @param precision the precision of the result
         * @return floating point result
         */
        void getEuclidean(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);

        /*!
         * Function that computes the Chebeshev distance.
         * Uses a different implementation than minkowski
         *
         * @param vector1
         * @param vector2
         * @return floating point result
         */
        void getChebyshev(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);

        //distances accepting strings for higher precision
        void getMinkowski(mpfr_t &res, vector<string> &vector1, vector<string> &vector2, double pvalue);
        void getMinkowskiInfP(mpfr_t &res, vector<string> &vector1, vector<string> &vector2, bool pos);
        void getManhattan(mpfr_t &res, vector<string> &vector1, vector<string> &vector2);
        void getEuclidean(mpfr_t &res, vector<string> &vector1, vector<string> &vector2);
        void getChebyshev(mpfr_t &res, vector<string> &vector1, vector<string> &vector2);

}}



