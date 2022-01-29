//
// Created by jack lewis on 03/11/2021.
//

#include <iostream>
#include <cmath>
#include "DMet/PointDistances.h"
#include <limits>
#include <mpfr.h>
#include <gmp.h>

void DMet::testGMP() {
    mpf_t x, y, result;
    mpf_init_set_d(x,7123134);
    mpf_init_set_d(y,2353134);
    mpf_init(result);
    mpf_mul(result, x, y);
    gmp_printf("    %.Ff\n"
               "*\n"
               "    %.Ff\n"
               "--------------------\n"
               "%.Ff\n", x, y, result);

    /* free used memory */
    mpf_clear(x);
    mpf_clear(y);
    mpf_clear(result);
}

void DMet::getMinkowski(mpf_t *res, double vector1[], double vector2[], int size1, int size2, unsigned long pvalue,long precision){
    mpfr_t sum;
    mpfr_init2(sum,precision);
    mpfr_set_d(sum,0,GMP_RNDN);
    mpfr_t x,y;
    mpfr_inits2(precision,x,y,NULL);
    for(int i=0;i<size1;i++){
        mpfr_set_d(x,vector1[i],GMP_RNDN);
        mpfr_set_d(y,vector2[i],GMP_RNDN);
        mpfr_sub(y,y,x,GMP_RNDN);
        mpfr_abs(y,y,GMP_RNDN);
        mpfr_add(sum,sum,y,GMP_RNDN);
        mpfr_printf("Sum %.11Re \n", sum);
    }
    mpfr_clears(x,y,NULL);

    mpfr_rootn_ui(sum,sum,pvalue,GMP_RNDN);
    mpfr_printf("Result %.64Re \n", sum);
    mpfr_clear(sum);
}

/*!
 * Function that computes the distance between vector 1 and vector 2
 * Uses Minkowski distance
 * (∑(Xi−Yi)^p)^1/p
 *
 * @param vector1
 * @param vector2
 * @param size1 size of vector1
 * @param size2 size of vector2
 * @param pvalue pvalue
 * @return floating point result
 */
double DMet::getMinkowski(double vector1[], double vector2[], int size1, int size2, double pvalue){
    if(vector1 == nullptr || vector2 == nullptr){
        std::cout << "Array cannot be NULL" << std::endl;
        return -1;
    }
    if(size1 != size2){
        std::cout << "Array sizes not compatable" << std:: endl;
        return -1;
    }
    if(pvalue == 0){
        //throw exception
    }
    //if overflow return infinty
    double sum = 0;
    for(int i=0; i<size1;i++){
        sum += pow(std::abs((vector1[i] - vector2[i])) ,pvalue);
    }
    double res = pow(sum, float(1.0f/pvalue));
    return res;
}

/*!
 * Function that computes the Euclidean distance.
 * Uses the implementation of Minkowski distance with a p value set as 2
 *
 * @param vector1
 * @param vector2
 * @param size1 size of vector1
 * @param size2 size of vector2
 * @return floating point result
 */
double DMet::getEuclidean(double vector1[], double vector2[], int size1, int size2){
    return getMinkowski(vector1, vector2, size1,size2,2);
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
double DMet::getManhattan(double vector1[], double vector2[], int size1, int size2){
    return getMinkowski(vector1, vector2, size1,size2,1);
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
double DMet::getChebyshev(double vector1[], double vector2[], int size1, int size2) {
    if(vector1 == NULL || vector2 == NULL){
        std::cout << "Array cannot be NULL" << std::endl;
        return -1;
    }
    if(size1 != size2){
        std::cout << "Array sizes not compatable" << std:: endl;
        return -1;
    }

    int max = 0;
    for(int i=0;i<size1;i++){
        if(max == NULL){
            max = std::abs(vector1[i] - vector2[i]);
        }else{
            int temp = std::abs(vector1[i] - vector2[i]);
            max = (temp > max) ? temp : max;
        }
    }
    return max;
}