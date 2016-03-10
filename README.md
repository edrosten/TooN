# TooN

## Compiling and installing

To install on a unix system:

./configure && make && sudo make install

To verify that everything works, you can optioinally run

make test


If you use any LAPACK based features (SVD.h, LU.h, QR_Lapack.h,
SymEigen.h, Lapack_Cholesky.h) you will need to link against LAPACK,
probably using -llapack and perhaps -lblas.

## System compatibility

The code is ready to go and should work on any system (unix or non-unix)
without configuring or compiling.


## Documentation

[![Documentation Status](https://codedocs.xyz/edrosten/TooN.svg)](https://codedocs.xyz/edrosten/TooN/)

Latest documentation here: https://codedocs.xyz/edrosten/TooN/ or just run Doxygen.

Documentation for latest release is here: http://www.edwardrosten.com/cvd/toon/html-user/index.html


## Status of unit tests

[![Build Status](https://drone.io/edrosten/TooN/status.png)](https://drone.io/edrosten/TooN/latest)

