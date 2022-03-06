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
            long precision = 200;
            //check arrays are compatible sizes
            if(vector1.empty() || vector2.empty()){
                throw std::invalid_argument("vector is empty");
            }
            if(vector1.size() != vector2.size()){
                throw std::invalid_argument("vector1 and vector2 are incompatible sizes");
            }
            //check pvalue is +inf or -infinity
            if(pvalue == inf){
                getChebyshev(res,vector1,vector2);
                return;
            }else if(pvalue == -inf){
                //oposite of chebshev
            }

            mpfr_t sum,x,y,p;
            mpfr_inits2(precision,sum,x,y,NULL);
            mpfr_set_d(sum,0,GMP_RNDN);

            for(int i=0; i<vector1.size(); i++){
                // inf - inf is undefined so set sum to nan
                if((vector1[i] == inf)
                && (vector2[i] == inf)){
                    mpfr_set_d(res,NAN,GMP_RNDN);
                    mpfr_clears(x,y,p,sum,NULL);
                    break;
                }
                else if(mpfr_get_d(sum,GMP_RNDN) == inf){
                    continue;
                }
                else if(vector1[i] == inf || vector2[i] == inf){
                    mpfr_set_d(res,inf,GMP_RNDN);
                    continue;
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
            if(vector1.size() != vector2.size()){
                std::cout << "Array sizes not compatable" << std:: endl;
                return;
            }

            mpfr_t max,x ,y;
            mpfr_init_set_d(max,0,GMP_RNDN);
            mpfr_inits(x,y,NULL);

            for(int i=0;i<vector1.size();i++){
                mpfr_set_d(x,vector1[i],GMP_RNDN);
                mpfr_set_d(y,vector2[i],GMP_RNDN);
                mpfr_sub(y,x,y,GMP_RNDN);
                mpfr_abs(y,y,GMP_RNDD);
                if(mpfr_cmp(y,max)>0){
                    mpfr_set(max,y,GMP_RNDN);
                }
            }
//            mpfr_printf("Result %.64Re \n", max);
            mpfr_set(res,max,GMP_RNDN);
            mpfr_clears(x,y,max,NULL);
        }
} }
