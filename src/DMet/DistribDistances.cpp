//
// Created by jack lewis on 22/11/2021.
//
#include <iostream>
#include <math.h>
#include "DMet/DistribDistance.h"
#include "DMet/utils.h"

namespace DMet { namespace Distrib{

    /*!
     *
     * @param x - this is the probability distribution x
     * @param y  - this is the probability distribution y
     * @return the KL distance between the two distributions
     */
    double KLDiv(double vector1[], double vector2[], int size1, int size2){
        //test both add to 1
        if(size1 != size2){
            std::cout << "Input vectors are incompattable sizes" << std::endl;
            return -1;
        }
        std::cout << size1 << std::endl;
        if(!DMet::Utils::isPDF(vector1,size1) || !DMet::Utils::isPDF(vector2,size2)){
            std::cout << "Input arrays not a PDF" << std::endl;
            return -1;
        }
        if(DMet::Utils::containsZero(vector2,size2)){
            std::cout << "Array 2 contains zero value" << std::endl;
            return -1;
        }

        double sum = 0;
        double prob1 = 0;
        double prob2 = 0;

        for(int i=0; i<size1; i++){
            double res = vector1[i] * log(vector1[i] / vector2[i]);
            std::cout << "loop " << i << " " << res << std::endl;
            sum += vector1[i] * log(vector1[i] / vector2[i]);
            prob1 += vector1[i];
            prob2 += vector2[i];
        }
        float epsilon = 0.001;
        if(prob1 < 1.0 - epsilon || prob1 > 1.0 + epsilon || prob2 < 1.0 - epsilon || prob2 > 1.0 + epsilon){
            std::cout << "Input is not a probability distribution" << std::endl;
            return -1;
        }

        //for all in size
        std::cout << "Distance: " << sum <<  std::endl;
        return sum;
    }



}}