#ifndef GlobalsMod
#define GlobalsMod

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define true 1
#define false 0

// Optional: double precision
//#define DoublePrec

#ifdef DoublePrec
typedef double RealType;
typedef double complex ComplexType;
#else
typedef float RealType;
//typedef float complex ComplexType;
typedef float _Complex ComplexType;
#endif

#endif
