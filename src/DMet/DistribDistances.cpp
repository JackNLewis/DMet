//
// Created by jack lewis on 22/11/2021.
//
#include <iostream>
#include <mpfr.h>
#include <DMet/EqWidthBin.h>
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
    void KLDiv(mpfr_t &res, vector<double> &v1, vector<double> &v2, bool pdfCheck){
        if(v1.size() != v2.size()){
            cout << "Input vectors are incompattable sizes" << endl;
            throw std::invalid_argument("Input vector sizes are incompatible");
        }
        if(!DMet::Utils::isPDF(v1) || !DMet::Utils::isPDF(v2)){
            if(pdfCheck){
                cout << "Input arrays not a PDF" << endl;
                throw std::invalid_argument("Input array is not a pdf");
            }
        }

        mpfr_t sum, x,y;
        mpfr_init_set_d(sum ,0 ,GMP_RNDN);
        mpfr_inits(x,y,NULL);

        for(int i=0; i<v1.size(); i++){

            if(v1[i] == 0 && v2[i]>=0){ // if divides by zero continue
                continue;
            }
            else if(v1[i] >0 && v2[i] == 0){
                cout<< "set to inf" <<endl;
                mpfr_set_d(res,std::numeric_limits<double>::infinity(),GMP_RNDN);
                return;
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
        mpfr_set(res,sum,GMP_RNDN);
        mpfr_clears(sum,x,y,NULL);
    }

    void KLDiv(mpfr_t &res, vector<double> &v1, vector<double> &v2){
        KLDiv(res,v1,v2,true);
    }

    void KLDiv(mpfr_t &res, vector<vector<double>> &v1, vector<vector<double>> &v2, int arity){
        DMet::EqWidthBin bin = DMet::EqWidthBin();
        vector<vector<double>> copy (v1);
        copy.insert(copy.end(),v2.begin(),v2.end());
        bin.setRanges(copy);
        bin.generateBins(arity);
        bin.assignBins(v1);
        vector<double> pdf1 = bin.getPDF();
        for(double d: pdf1){
            cout << d << endl;
        }
        cout<<std::endl;
        bin.clearBins();
        bin.assignBins(v2);
        vector<double> pdf2 = bin.getPDF();
        for(double d: pdf2){
            cout << d << endl;
        }
        KLDiv(res,pdf1,pdf2,true);
        mpfr_printf("Result: %.5Re\n",res);
    }

    void JensenShannon(mpfr_t &res, vector<double> &v1, vector<double> &v2){
        if(v1.size() != v2.size()){
            cout << "Input vectors are incompattable sizes" << endl;
            throw std::invalid_argument("Input vector sizes are incompatible");
        }
        if(!DMet::Utils::isPDF(v1) || !DMet::Utils::isPDF(v2)){
            cout << "Input arrays not a PDF" << endl;
            throw std::invalid_argument("Input array is not a pdf");
        }

        vector<double> m;
        for(int i=0;i<v1.size();i++){
            double mval = 0.5*(v1[i] +v2[i]);
            m.push_back(mval);
        }
        mpfr_t kl1,kl2;
        mpfr_inits(kl1, kl2,NULL);

        KLDiv(kl1, v1, m);
        KLDiv(kl2, v2, m);

        mpfr_add(kl1,kl1,kl2,GMP_RNDN);
        mpfr_div_d(kl1,kl1,2,GMP_RNDN);

        mpfr_set(res,kl1,GMP_RNDN);

        mpfr_clears(kl1,kl2,NULL);
    }

    void JensenShannon(mpfr_t &res, vector<vector<double>> &v1, vector<vector<double>> &v2,int arity){
        DMet::EqWidthBin bin = DMet::EqWidthBin();
        vector<vector<double>> copy (v1);
        copy.insert(copy.end(),v2.begin(),v2.end());
        bin.setRanges(copy);
        bin.generateBins(arity);
        bin.assignBins(v1);
        vector<double> pdf1 (bin.getPDF());
        bin.clearBins();
        bin.assignBins(v2);
        vector<double> pdf2 (bin.getPDF());
        cout << "[" << std::flush;
        for(double d: pdf1){
            cout << d << "," << std::flush;
        }
        cout << "]" << endl;
        cout << "[" << std::flush;
        for(double d: pdf2){
            cout << d << "," << std::flush;
        }
        cout << "]" << endl;
        JensenShannon(res,pdf1,pdf2);
    }
}}