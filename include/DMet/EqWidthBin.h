//
// Created by jack lewis on 21/01/2022.
//

#pragma once
#include <fstream>
#include <vector>
#include <ostream>

using std::vector;

namespace DMet {
    /*!
     * Class for equal width binning
     */
    class EqWidthBin {
        public:

        /*!
         * Structure for a bin
         */
        struct Bin {
            /** Contains the vector of ranges for each dimension **/
            vector<vector<double>> range;

            /** Contains a list of the raw values for each dimension **/
            vector<vector<double>> values;

            /** override operator << in order to print the bin on screen **/
            friend std::ostream& operator<< (std::ostream& out, const Bin &point);
        };

        /**
         * Stores the ranges of each dimension
         */
        vector<vector<double>> ranges;

        /**
         * Stores the different vector of bins
         */
        vector<Bin> bins;

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
         * Assign the values from parameter data to the bins that have been generated.
         * Throws an error if the bins have not been generated.
         *
         * @param data
         */
        void assignBins(vector<vector<double>>& data);

        /*!
         * Assingns a single point to its correct bin
         * @param point
         */
        void assignPoint(vector<double>& point);

        /*!
         * Clears all of the values in the bins so each bin is empty
         */
        void clearBins();

        /*!
         * Returns the probability density function of the binned values
         * @return
         */
        vector<double> getPDF();
        /*!
         * Generates the cartesian product of two lists of vectors
         *
         * @param vec1
         * @param vec2
         * @return
         */
         //Move to utility
        vector<vector<vector<double>>> cartesian(vector<vector<vector<double>>> &vec1, vector<vector<double>> &vec2);

    };
}


