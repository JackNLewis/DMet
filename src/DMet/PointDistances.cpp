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
         * uses precision 200
         * (∑(Xi−Yi)^p)^1/p
         *
         * @param vector1
         * @param vector2
         * @param size1 size of vector1
         * @param size2 size of vector2
         * @param pvalue pvalue
         * @return floating point result
         */
        void getMinkowski(mpfr_t &res, vector<double> &vector1, vector<double> &vector2, double pvalue){
            long precision = 200;
            mpfr_t sum;
            mpfr_init2(sum,precision);
            mpfr_set_d(sum,0,GMP_RNDN);
            mpfr_t x,y,p;
            mpfr_inits2(precision,x,y,NULL);

            for(int i=0; i<vector1.size(); i++){
                mpfr_set_d(x,vector1[i],GMP_RNDN);
                mpfr_set_d(y,vector2[i],GMP_RNDN);
                mpfr_sub(y,x,y,GMP_RNDN);
                mpfr_init_set_d(p,pvalue,GMP_RNDD);
                mpfr_abs(y,y,GMP_RNDD);
                mpfr_pow(y,y,p,GMP_RNDN);
                mpfr_add(sum,sum,y,GMP_RNDN);
//                mpfr_printf("sum %.64Re \n", sum);
            }

            mpfr_clears(x,y,p,NULL);
            mpfr_rootn_ui(sum,sum,pvalue,GMP_RNDN);
            mpfr_printf("Result %.64Re \n", sum);
            mpfr_set(res,sum,GMP_RNDN);
            mpfr_clear(sum);
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
        void getEuclidean(mpfr_t &res, vector<double> &vector1, vector<double> &vector2){
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
         * @param vector1
         * @param vector2
         * @param size1 size of vector1
         * @param size2 size of vector2
         * @return floating point result
         */
         //Test for infinite values
        void getChebyshev(mpfr_t &res, vector<double> &vector1, vector<double> &vector2){
            if(vector1.size() != vector2.size()){
                std::cout << "Array sizes not compatable" << std:: endl;
                return;
            }

            double max = 0;
            for(int i=0;i<vector1.size();i++){
                double temp = std::abs(vector1[i] - vector2[i]);
                max = (temp > max) ? temp : max;

            }
            std::cout << "Chebyshev distance: " << max << std::endl;
        }
} }




//double DMet::PointDistances::getMinkowski(double vector1[], double vector2[], int size1, int size2, double pvalue){
//    if(vector1 == nullptr || vector2 == nullptr){
//        std::cout << "Array cannot be NULL" << std::endl;
//        return -1;
//    }
//    if(size1 != size2){
//        std::cout << "Array sizes not compatable" << std:: endl;
//        return -1;
//    }
//    if(pvalue == 0){
//        //throw exception
//    }
//    //if overflow return infinty
//    double sum = 0;
//    for(int i=0; i<size1;i++){
//        sum += pow(std::abs((vector1[i] - vector2[i])) ,pvalue);
//    }
//    double res = pow(sum, float(1.0f/pvalue));
//    return res;
//}
//


//
//
