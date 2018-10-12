#ifndef LEGENDRE_POLYNOMIALS_H
#define LEGENDRE_POLYNOMIALS_H

#include "global_definitions.h"
#include "LargeFactorHelper.h"

class LegendrePolynom //uses recurrence relation, plus saves values for 2 maximum ls for faster evaluation of the next l (faster performance for the same x and rising l)
{
protected:
	//these are for caching (improves speed for P(x) for the same x and rising l)
	unsigned int l_last;
	unsigned int l_last1;
	long double P_last;
	long double P_last1;
	long double x_last;
public:
	LegendrePolynom();
	long double operator ()(long double x, unsigned int l);
};

long double Int_PlPl_1_0 (unsigned int n, unsigned int k);
long double Int_PlPl_0_1 (unsigned int n, unsigned int k);
long double Int_PlPl_transf_1_0 (unsigned int n, unsigned int k);
long double Int_PlPl_transf_0_1 (unsigned int n, unsigned int k);
long double Int_PlPl (unsigned int n, unsigned int k, long double from, long double to, long double dx);
long double Int_PlPl_transf (unsigned int n, unsigned int k, long double from, long double to, long double dx);

#endif
