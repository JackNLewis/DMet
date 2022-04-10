.. DMet documentation master file, created by
   sphinx-quickstart on Sat Jan 29 19:34:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DMet's documentation!
================================
DMet's purpose is to compute various distance metrics such difference between points and distributions.
The use of arbitrary precision floating point numbers means it can work out distances which would usually result in an overflow or underflow error.
The library also contains a binning class which enables you to input raw data up to n dimensions and produces the
probability density function.

There is a simple code snippet on how to use MPFR below.
Follow this link to find out more about MPFR: https://www.mpfr.org/mpfr-current/mpfr.html

.. code-block:: c

   //decleare and initialise variables
   mpfr_t res;
   mpfr_init(res);

   //use library to comupte distances e.g. euclidean distance
   vector<double> v1 {3.0,1.0,1.0};
   vector<double> v2 {2.0,1.0,3.0};
   DMet::getEuclidean(res,v1,v2);

   //Can print value using
   mpfr_printf("Result %.5Re\n", mpfr_variable);

   //if it fits in a double can convert it
   double d = mpfr_get_d(res, GMP_RNDN);

   //clear variables
   mpfr_clear(res);

View the full code on `Github`_.

.. _Github: https://github.com/JackNLewis/DMet
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pointdocs
   distribdocs
   binningdocs




