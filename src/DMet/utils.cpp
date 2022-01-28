//
// Created by jack lewis on 22/11/2021.
//
#include <vector>
#include <iostream>
#include <optional>
#include <utility>
#include <cmath>

namespace DMet{ namespace Utils{

    bool isPDF(double arr[], int size1){
        double sum = 0;
        for(int i=0;i<size1;i++){
            if(arr[i] <0 ){
                return false;
            }
            sum += arr[i];
        }
        if(abs(1.0 - sum) < 0.01){
            return true;
        }else{
            return false;
        }
    }


    bool containsZero(double arr[], int size1){
        for(int i=0;i<size1;i++){
            if(arr[i] == 0){
                return false;
            }
        }
        return true;
    }

}}
