.. DMet documentation master file, created by
   sphinx-quickstart on Sat Jan 29 19:34:40 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DMet's documentation!
================================
DMet's purpose is to compute various distance metrics such difference between points and distributions.
The use of aribtary presicion floating point numbers means it can work out distances which would usually result in an over or underflow error if you are working with
spesifically large or small floating point numbers. The library also contains a binning class enabling you to input raw data up to n dimensions and produces the
probability density function.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pointdocs
   distribdocs
   binningdocs




