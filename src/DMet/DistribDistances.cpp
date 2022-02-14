//
// Created by jack lewis on 22/11/2021.
//
#include <iostream>
#include <math.h>
#include <mpfr.h>
#include "DMet/DistribDistance.h"
#include "DMet/utils.h"

using std::cout;
using std::endl;
using std::vector;
namespace DMet { namespace Distrib{

    /*!
     *
     * @param x - this is the probability distribution x
     * @param y  - this is the probability distribution y
     * @return the KL distance between the two distributions
     */
    void KLDiv(mpfr_t &res, vector<double> &v1, vector<double> &v2){
        //test both add to 1
        if(v1.size() != v2.size()){
            cout << "Input vectors are incompattable sizes" << endl;
            return;
        }
        if(!DMet::Utils::isPDF(v1) || !DMet::Utils::isPDF(v2)){
            cout << "Input arrays not a PDF" << endl;
            return;
        }
        if(DMet::Utils::containsZero(v2)){
            cout << "Array 2 contains zero value" << endl;
            return;
        }

        mpfr_t sum, x,y;
        mpfr_init_set_d(sum ,0 ,GMP_RNDN);
        mpfr_inits(x,y,NULL);

        for(int i=0; i<v1.size(); i++){
            if(v2[i] == 0){
                continue;
            }
            mpfr_set_d(x,v1[i],GMP_RNDN);
            mpfr_set_d(y,v2[i],GMP_RNDN);

            // res = v1[i] * log(v1[i] / v2[i]);
            mpfr_div(y,x,y,GMP_RNDN);
            mpfr_log(y,y,GMP_RNDN);
            mpfr_mul(x,x,y,GMP_RNDN);


            //sum += res;
            mpfr_add(sum,sum,x,GMP_RNDN);
        }
//        mpfr_printf("Result %.64Re \n", sum);
        mpfr_set(res,sum,GMP_RNDN);
        mpfr_clears(sum,x,y,NULL);
    }

//    void JennsonShannon(mpfr_t &res, vector<double> &vector1, vector<double> &vector2){
//        mpfr_t
//    }
}}