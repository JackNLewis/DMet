//
// Created by jack lewis on 03/11/2021.
//

#pragma once
#include <gmp.h>
namespace DMet{
    int add(int a, int b);
    double getEuclidean(double vector1[], double vector2[], int size1, int size2);
    double getMinkowski(double vector1[], double vector2[], int size1, int size2, double pvalue);
    double getManhattan(double vector1[], double vector2[], int size1, int size2);
    double getChebyshev(double vector1[], double vector2[], int size1, int size2);
    void testGMP();

    void
    getMinkowski(mpf_t *res, double *vector1, double *vector2, int size1, int size2, unsigned long pvalue, long precision);
}



