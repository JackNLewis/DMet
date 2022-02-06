//
// Created by jack lewis on 22/11/2021.
//
#include <vector>
#include <iostream>
#include <optional>
#include <utility>
#include <cmath>

using std::vector;
using std::cout;
using std::endl;

namespace DMet{ namespace Utils{

    bool isPDF(vector<double> &v){
        double sum = 0;
        for(int i=0;i<v.size();i++){
            if(v[i] < 0 ){
                return false;
            }
            sum += v[i];
        }
        if(abs(1.0 - sum) < 0.01){
            return true;
        }else{
            return false;
        }
    }


    bool containsZero(std::vector<double> &v){
        for(double i : v){
            cout << i <<endl;
            if(i == 0){
                return true;
            }
        }
        return false;
    }

}}
