//
// Created by jack lewis on 03/11/2021.
//

#include <iostream>
#include <cmath>
#include "DMet/PointDistances.h"
#include <limits>
#include <mpfr.h>
#include <gmp.h>
#include <vector>


namespace DMet { namespace PointDistances{
        using std::vector;

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
        void getMinkowski(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, double pvalue){
            double inf = std::numeric_limits<double>::infinity();
            //check arrays are compatible sizes
            if(vector1.empty() || vector2.empty()){
                throw std::invalid_argument("vector is empty");
            }
            if(vector1.size() != vector2.size()){
                throw std::invalid_argument("vector1 and vector2 are incompatible sizes");
            }

            //check pvalue is +inf or -infinity
            if(pvalue == inf){
                getMinkowskiInfP(res,vector1,vector2,true);
                return;
            }else if(pvalue == -inf){
                getMinkowskiInfP(res,vector1,vector2,false);
                return;
            }

            mpfr_t sum,x,y,p;
            mpfr_inits(sum,x,y,p,NULL);
            mpfr_set_d(sum,0,GMP_RNDN);

            for(int i=0; i<vector1.size(); i++){
                // inf - inf is undefined so set sum to nan
                if((vector1[i] == inf)
                && (vector2[i] == inf)){
                    mpfr_set_nan(res);
                    mpfr_clears(x,y,p,sum,NULL);
                    return;
                }
                if((vector1[i] == inf)
                   || (vector2[i] == inf)){
                    mpfr_set_inf(res,1);
                    mpfr_clears(x,y,p,sum,NULL);
                    return;
                }
                mpfr_set_d(x,vector1[i],GMP_RNDN);
                mpfr_set_d(y,vector2[i],GMP_RNDN);
                mpfr_sub(y,x,y,GMP_RNDN);
                mpfr_init_set_d(p,pvalue,GMP_RNDD);
                mpfr_abs(y,y,GMP_RNDD);
                mpfr_pow(y,y,p,GMP_RNDN);
                mpfr_add(sum,sum,y,GMP_RNDN);
            }
            mpfr_rootn_ui(sum,sum,pvalue,GMP_RNDN);
            mpfr_set(res,sum,GMP_RNDN);
            mpfr_clears(x,y,p,sum,NULL);
        }

        /**
         * Calculates minkowski distance when the p value is either +inf or -inf. This is equivalent to finding the
         * largest distance when p is +inf or finding the smallest distance when p is -inf.
         *
         * @param res
         * @param vector1
         * @param vector2
         * @param pos
         */
        void getMinkowskiInfP(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, bool pos){
            mpfr_t max_min,x ,y;
            mpfr_init_set_d(max_min, 0, GMP_RNDN);
            mpfr_inits(x,y,NULL);

            for(int i=0;i<vector1.size();i++){
                mpfr_set_d(x,vector1[i],GMP_RNDN);
                mpfr_set_d(y,vector2[i],GMP_RNDN);
                mpfr_sub(y,x,y,GMP_RNDN);
                mpfr_abs(y,y,GMP_RNDD);
                //if positive p and is greater than max
                //or negative p and less than min
                if((pos && mpfr_greater_p(y, max_min)) || (!pos && mpfr_less_p(y, max_min))){
                    mpfr_set(max_min, y, GMP_RNDN);
                }
            }
            mpfr_set(res, max_min, GMP_RNDN);
            mpfr_clears(x, y, max_min, NULL);
        }

        /*!
         * Function that computes the Euclidean distance.
         * Uses the implementation of Minkowski distance with a p value set as 2
         *
         * @param vector1 input vector 1
         * @param vector2 input  vector 2
         * @param precision the precision of the result
         * @return floating point result
         */
        void getEuclidean(mpfr_t &res, std::vector<double> &vector1, std::vector<double> &vector2){
            getMinkowski(res, vector1, vector2,2);
        }

        /*!
         * Function that computes the Manhattan distance.
         * Uses the implementation of Minkowski distance with a p value set as 1
         * @param vector1
         * @param vector2
         * @param size1 size of vector1
         * @param size2 size of vector2
         * @return floating point result
        */
        void getManhattan(mpfr_t &res, vector<double> &vector1, vector<double> &vector2){
            getMinkowski(res, vector1, vector2 ,1);
        }

        /*!
         * Function that computes the Chebeshev distance.
         * Uses a different implimentation than minkowski
         *
         * @param vector1
         * @param vector2
         * @return floating point result
         */
         //Test for infinite values
        void getChebyshev(mpfr_t &res, vector<double> &vector1, vector<double> &vector2){
             if(vector1.empty() || vector2.empty()){
                 throw std::invalid_argument("vector is empty");
             }
             if(vector1.size() != vector2.size()){
                 throw std::invalid_argument("vector1 and vector2 are incompatible sizes");
             }
             getMinkowskiInfP(res,vector1,vector2,true);
        }


        //==========================================================//
        //===========OVERLOADED METHODS WITH STRING INPUT===========//
        //==========================================================//

        /**
         * Retrieves Minkowski distance. Accepts vector of strings as inputs if higher precision than double is needed.
         *
         * @param res
         * @param vector1
         * @param vector2
         * @param pvalue
         */
        void getMinkowski(mpfr_t &res, vector<string> &vector1, vector<string> &vector2, double pvalue){
            double inf = std::numeric_limits<double>::infinity();
            //check arrays are compatible sizes
            if(vector1.empty() || vector2.empty()){
                throw std::invalid_argument("vector is empty");
            }
            if(vector1.size() != vector2.size()){
                throw std::invalid_argument("vector1 and vector2 are incompatible sizes");
            }

            //check pvalue is +inf or -infinity
            if(pvalue == inf){
                getMinkowskiInfP(res,vector1,vector2,true);
                return;
            }else if(pvalue == -inf){
                getMinkowskiInfP(res,vector1,vector2,false);
                return;
            }

            mpfr_t sum,x,y,p;
            mpfr_inits(sum,x,y,p,NULL);
            mpfr_set_d(sum,0,GMP_RNDN);

            for(int i=0; i<vector1.size(); i++){

                // inf - inf is undefined so set sum to nan
                if(mpfr_set_str(x,vector1[i].c_str(),10,GMP_RNDN) == -1 ||
                        mpfr_set_str(y,vector2[i].c_str(),10,GMP_RNDN) == -1){
                    throw std::invalid_argument("String input not a valid float");
                }

                if(mpfr_inf_p(x) && mpfr_inf_p(y)){
                    mpfr_set_nan (res);
                    mpfr_clears(x,y,p,sum,NULL);
                    return;
                }
                if(mpfr_inf_p(x) && mpfr_inf_p(y)){
                    mpfr_set_nan (res);
                    mpfr_clears(x,y,p,sum,NULL);
                    return;
                }

                mpfr_sub(y,x,y,GMP_RNDN);
                mpfr_init_set_d(p,pvalue,GMP_RNDD);
                mpfr_abs(y,y,GMP_RNDD);
                mpfr_pow(y,y,p,GMP_RNDN);
                mpfr_add(sum,sum,y,GMP_RNDN);
            }
            mpfr_rootn_ui(sum,sum,pvalue,GMP_RNDN);
            mpfr_set(res,sum,GMP_RNDN);
            mpfr_clears(x,y,p,sum,NULL);
        }

        void getMinkowskiInfP(mpfr_t &res, vector<string> &vector1, vector<string> &vector2, bool pos){
            mpfr_t max_min,x ,y;
            mpfr_init_set_d(max_min, 0, GMP_RNDN);
            mpfr_inits(x,y,NULL);

            for(int i=0;i<vector1.size();i++){
                if(mpfr_set_str(x,vector1[i].c_str(),10,GMP_RNDN) == -1 ||
                mpfr_set_str(y,vector2[i].c_str(),10,GMP_RNDN) == -1){
                    throw std::invalid_argument("String is not a valid float");
                }

                mpfr_sub(y,x,y,GMP_RNDN);
                mpfr_abs(y,y,GMP_RNDD);
                //if positive p and is greater than max
                //or negative p and less than min
                if((pos && mpfr_greater_p(y, max_min)) || (!pos && mpfr_less_p(y, max_min))){
                    mpfr_set(max_min, y, GMP_RNDN);
                }
            }
            mpfr_set(res, max_min, GMP_RNDN);
            mpfr_clears(x, y, max_min, NULL);
        }

        void getManhattan(mpfr_t &res, vector<string> &vector1, vector<string> &vector2){
            getMinkowski(res,vector1,vector2,1);
        }

        void getEuclidean(mpfr_t &res, vector<string> &vector1, vector<string> &vector2){
            getMinkowski(res,vector1,vector2,2);
        }

        void getChebyshev(mpfr_t &res, vector<string> &vector1, vector<string> &vector2){
            getMinkowskiInfP(res,vector1,vector2,true);
        }
} }
