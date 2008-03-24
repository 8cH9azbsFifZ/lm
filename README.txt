Code Description:
=================
This code solves numerically the monodisperse and binary OZ / PY/HNC equations 
in two dimensions. The program was written in 2005 as a part of my diploma
thesis.

Status:
=======
This is a pure working code, which has to be adjusted in regard to potential
and units for your own needs. There are nearly no comments and some variable
descriptions are in german.

Dependencies:
=============
The code uses the free iterator lmfit[1]. This code has been put together
in marquardt.c. For the original version have a look a the reference.
Furthermore the GSL[2] can be used optionally for calculating the zeros of 
the bessel functions.

1: J. Wuttke: http://sourceforge.net/projects/lmfit
2: Gnu Scientific Library: http://www.gnu.org/software/gsl
