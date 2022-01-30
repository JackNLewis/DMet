//
// Created by jack lewis on 21/01/2022.
//

#pragma once
#include <fstream>
#include <vector>

using std::vector;

namespace DMet {
    /*!
     * Class for equal width binning
     */
    class EqWidthBin {

        struct Bin {
            int index;
            vector<vector<double>> range;
            vector<vector<double>> values;
        };

    public:
        /**
         * Stores the ranges of each dimension
         */
        vector<vector<double>> ranges;

        /**
         * Stores the different vector of bins
         */
        vector<vector<vector<double>>> bins;

        /*!
         * Method to set the ranges for every dimension given data
         * @param data
         * @return
         */
        vector<vector<double>> setRanges(vector<vector<double>> &data);

        /*!
         * Generates the bins from the ranges for each dimension then creates the cartesian product of these bins
         *
         * @param arity
         */
        void generateBins(int arity);

        /*!
         * Generates the cartesian product of two lists of vectors
         *
         * @param vec1
         * @param vec2
         * @return
         */
        vector<vector<vector<double>>> cartesian(vector<vector<vector<double>>> &vec1, vector<vector<double>> &vec2);

    };
}


