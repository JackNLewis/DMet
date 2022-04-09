# DMet

DMet is an open source C++ library made to compute distance functions with arbitary precision. Features of DMet include multidimesnional binning, compatability
with infinite values and also arbitary floating point precision. The distances functions contained in DMet are:
- Minkowski 
- Manhattan 
- Euclidean
- Chebyshev
- Kullback Leibler
- Jensen Shannon

## Installation
DMet requires the library GMP and MPFR for the arbitary floating point precision which need to be installed seperately. You can use homebrew to install
these easily:

        brew install gmp
        
        brew install mpfr
      
These should be installed in the directorys

/usr/local/Cellar/gmp/6.2.1_1/lib/libgmp.a  
/usr/local/Cellar/mpfr/4.1.0/lib/libmpfr.a

Alternatively you can download it form the websites:

https://gmplib.org/#DOWNLOAD 

https://www.mpfr.org/mpfr-current/#download

and follow the instuctions to install from there.

Clone the Github

        git clone https://git-teaching.cs.bham.ac.uk/mod-ug-proj-2021/jxl1358.git
        
## Running tests from terminal
Move into the DMet project 

        cd DMet
        
Then build the project using CMake and place files in a new build directory

        cmake -S . -B build

To run test files go into build test folder

        cd build/test

Then run

        make
        
Then call the google tests executable 

        ./Google_tests_run
   
