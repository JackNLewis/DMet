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

        /*!
        * Function that computes the distance between two points represented by their cartesian coordinated
        * Uses Minkowski distance
        * Uses precision 200
        * (∑|Xi−Yi|^p)^1/p
        *
        * @param res mpfr variable to store result
        * @param vector1 point 1
        * @param vector2 point 2
        * @param pvalue pvalue
        */
        void getMinkowski(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, double pvalue);

        /**
         * Calculates the minkowski distance between two points represented by their cartesian coordinates
         * when the p value is either +inf or -inf. This is equivalent to finding the
         * largest distance when p is +inf or finding the smallest distance when p is -inf.
         *
         * @param res mpfr variable to store result
         * @param vector1 point 1
         * @param vector2 point 2
         * @param pos flag to determine if Infinity is positive or negative. For positive set to true.
         */
        void getMinkowskiInfP(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, bool pos);

        /*!
         * Function that computes the Manhattan distance between two points .
         * Uses the implementation of Minkowski distance with a p value set as 1
         *
         * @param res mpfr variable to store result
         * @param vector1 point 1
         * @param vector2 point 2
        */
        void getManhattan(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);

        /*!
         * Function that computes the Euclidean distance between 2 points represented as cartesian coordinates.
         * Uses the implementation of Minkowski distance with a p value set as 2
         *
         * @param res mpfr variable to store result
         * @param vector1 point 1
         * @param vector2 input  vector 2
         */
        void getEuclidean(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);

        /*!
         * Function that computes the Chebyshev distance between 2 points represented as cartesian coordinates.
         *
         * @param res mpfr variable to store result
         * @param vector1 point 1
         * @param vector2 point 2
         */
        void getChebyshev(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);

        //distances accepting strings for higher precision
        /***
        * Function that computes the distance between two points represented by their cartesian coordinated
        * Uses Minkowski distance
        * Uses precision 200
        * (∑|Xi−Yi|^p)^1/p
        *
        * @param res mpfr variable to store result
        * @param vector1 point 1
        * @param vector2 point 2
        * @param pvalue pvalue
        */
        void getMinkowski(mpfr_t &res, vector<string> &vector1, vector<string> &vector2, double pvalue);

        /**
         * Calculates the minkowski distance between two points represented by their cartesian coordinates
         * when the p value is either +inf or -inf. This is equivalent to finding the
         * largest distance when p is +inf or finding the smallest distance when p is -inf.
         *
         * @param res mpfr variable to store result
         * @param vector1 point 1
         * @param vector2 point 2
         * @param pos flag to determine if Infinity is positive or negative. For positive set to true.
         */
        void getMinkowskiInfP(mpfr_t &res, vector<string> &vector1, vector<string> &vector2, bool pos);

        /*!
         * Function that computes the Manhattan distance between two points .
         * Uses the implementation of Minkowski distance with a p value set as 1
         *
         * @param res mpfr variable to store result
         * @param vector1 point 1
         * @param vector2 point 2
        */
        void getManhattan(mpfr_t &res, vector<string> &vector1, vector<string> &vector2);

        /*!
        * Function that computes the Euclidean distance between 2 points represented as cartesian coordinates.
        * Uses the implementation of Minkowski distance with a p value set as 2
        *
        * @param res mpfr variable to store result
        * @param vector1 point 1
        * @param vector2 input  vector 2
        */
        void getEuclidean(mpfr_t &res, vector<string> &vector1, vector<string> &vector2);

        /*!
        * Function that computes the Chebyshev distance between 2 points represented as cartesian coordinates.
        *
        * @param res mpfr variable to store result
        * @param vector1 point 1
        * @param vector2 point 2
        */
        void getChebyshev(mpfr_t &res, vector<string> &vector1, vector<string> &vector2);

}}



