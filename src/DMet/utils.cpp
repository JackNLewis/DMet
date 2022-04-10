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



    vector<vector<vector<double>>> cartesian(vector<vector<vector<double>>>& vec1, vector<vector<double>>& vec2){
        //no dimesnions yet
        if(vec1.size() == 0){
            for(auto v2 : vec2){
                vector<vector<double>> t{v2};
                vec1.push_back(t);
            }
        }

        //else multiply vec1 vec2.size times to get every combination
        else{
            vector<vector<vector<double>>> copy(vec1);
            int n = copy.size();
            for(int i=0;i<vec2.size()-1;i++){
                vec1.insert(vec1.end(),copy.begin(),copy.end());
            }

            int startPos = 0;
            for(auto d : vec2){
                //copy it n times into newly created vector res
                for(int i=0;i<n;i++){
                    vec1[startPos+i].insert(vec1[startPos+i].begin(), d);
                }
                startPos+=n;
            }
        }
        return vec1;
    }

}}
