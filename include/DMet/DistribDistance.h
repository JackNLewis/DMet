//
// Created by jack lewis on 22/11/2021.
//

#pragma once
#include <vector>
using std::vector;

namespace DMet { namespace Distrib{

    /***
     * Computes the Kullback Leibler Divergence between two probability distributions.
     *
     * ∑ Xi log(Xi/Yi)
     * if Xi = 0 && Yi >= 0 then Yi
     * if Xi > 0 && Yi = 0 then inf
     *
     * @param res mpfr variable to store result
     * @param v1 vector 1 representing probability distribution
     * @param v2 vector 2 representing probability distribution
     * @param pdfCheck Flag to turn of pdf check used for Jensen Shannon
     */
    void KLDiv(mpfr_t &res, vector<double> &v1, vector<double> &v2, bool pdfCheck);

    /***
     * Computes the Kullback Leibler Divergence between two probability distributions.
     *
     *  KL_Div = ∑ Xi log(Xi/Yi)
     * if Xi = 0 && Yi >= 0 then Yi
     * if Xi > 0 && Yi = 0 then inf
     *
     * @param res mpfr variable to store result
     * @param v1 vector 1 representing probability distribution
     * @param v2 vector 2 representing probability distribution
     */
    void KLDiv(mpfr_t &res, vector<double> &v1, vector<double> &v2);

    /***
     * Computes the Kullback Leibler Divergence between two sets of n dimensional points. Data is first binned to generate
     * probability distributions then used to calculate KL divergence. Data is binned using equal width binning where each dimension
     * is split into k amount of bins with k representing the arity.
     *
     * @param res mpfr variable to store result
     * @param v1 vector 1
     * @param v2 vector 2
     * @param arity Number of bins per dimension
     */
    void KLDiv(mpfr_t &res, vector<vector<double>> &v1, vector<vector<double>> &v2, int arity);


    /***
     * Computes the Jensen Shannon Divergence between 2 probability distributions
     *
     * JS(P||Q) = 1/2KL_Div(P||M) + 1/2KL_Div(Q||M)
     * M = 1/2(P+Q)
     * where KL_Div is the Kullback Leibler Divergence
     *
     * @param res the mpfr variable to store the result
     * @param vector1 probability distribution 1
     * @param vector2 probability distribution 2
     */
    void JensenShannon(mpfr_t &res, vector<double> &vector1, vector<double> &vector2);


    /***
     * Computes the Jensen Shannon Distance from 2 vectors of point data.
     * The data is first binned using equal width binning into k bins per dimension where k is the arity.
     * The probability distribution is then retrieved the Jensen Shannon divergence is calculated.
     *
     * @param res mpfr variable to store result
     * @param v1 vector 1
     * @param v2 vector 2
     * @param arity number of bins per dimension
     */
    void JensenShannon(mpfr_t &res, vector<vector<double>> &v1, vector<vector<double>> &v2, int arity);
}}

