//
// Created by jack lewis on 22/11/2021.
//
#include <iostream>
#include <math.h>
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
    double KLDiv(vector<double> &v1, vector<double> &v2){
        //test both add to 1
        if(v1.size() != v2.size()){
            cout << "Input vectors are incompattable sizes" << endl;
            return -1;
        }
        if(!DMet::Utils::isPDF(v1) || !DMet::Utils::isPDF(v2)){
            cout << "Input arrays not a PDF" << endl;
            return -1;
        }
        if(DMet::Utils::containsZero(v2)){
            cout << "Array 2 contains zero value" << endl;
            return -1;
        }

        double sum = 0;

        for(int i=0; i<v1.size(); i++){
            if(v2[i] == 0){
                continue;
            }
            double res = v1[i] * log(v1[i] / v2[i]);
            cout << "pos " << i << ": " << res << endl;
            sum += res;

        }

        //for all in size
        cout << "Distance: " << sum <<  endl;
        return sum;
    }



}}