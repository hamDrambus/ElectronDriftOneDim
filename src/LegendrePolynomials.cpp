#include "LegendrePolynomials.h"

LegendrePolynom::LegendrePolynom() {
	l_last = 0;
	l_last1 = 0;
	P_last = 1;
	P_last1 = 1;
	x_last = 1;
}

long double LegendrePolynom::operator ()(long double x, unsigned int l) {
	if (x==x_last) {
		if (l==l_last1)
			return P_last1;
		if (l==l_last)
			return P_last;
		if (l>l_last) {
			//iterate till l_last==l since l=l_last
			long double mem = 0;
			for (unsigned int i=(l_last+1);i<=l;++i) {
				mem = ((2*i-1)*x*P_last-(i-1)*P_last1)/i;
				P_last1 = P_last;
				P_last = mem;
				l_last1 = l_last; //==i-1
				l_last = i;
			}
			return P_last;
		}
		//in case l<l_last iterate since l=0:
	}
	l_last = 1;
	l_last1 = 0;
	P_last = x;
	P_last1 = 1;
	x_last = x;
	//iterate since l=0
	if (0==l)
		return 1;
	long double mem = 0;
	for (unsigned int i=2;i<=l;++i) {
		mem = ((2*i-1)*x*P_last-(i-1)*P_last1)/i;
		P_last1 = P_last;
		P_last = mem;
		l_last1 = l_last; //==i-1
		l_last = i;
	}
	return P_last;
}

//this formula is obtained by polynomial expansion of legendre polynomials
//Integral from -1 to 0
long double Int_PlPl_1_0 (unsigned int n, unsigned int k)
{
	if ((n%2)==(k%2))
		return (n==k ? 1/(long double)(2*n+1) : 0);
	LargeFactorHelper out;
	out.MultiplyBy(0);
	LargeFactorHelper Coef;
	for (unsigned int f = (0==n%2 ? n/2 : 1+ n/2); f<=n; ++f)
		for (unsigned int l = (0==k%2 ? k/2 : 1+ k/2); l<=k; ++l) {
			Coef.Clear();
			Coef.MultiplyByFactorial(2*f);
			Coef.MultiplyByFactorial(2*l);
			Coef.DivideByFactorial(2*f-n);
			Coef.DivideByFactorial(2*l-k);
			Coef.DivideByFactorial(f);
			Coef.DivideByFactorial(n-f);
			Coef.DivideByFactorial(l);
			Coef.DivideByFactorial(k-l);
			Coef.DivideBy(2*f-n+2*l-k+1);
			if (0==((2*f-n+2*l-k)%2))
				out += Coef;
			else
				out -= Coef;
		}
	out.DivideBy(2, n+k);
	return out.Output();
}

//Integral from 0 to 1
long double Int_PlPl_0_1 (unsigned int n, unsigned int k)
{
	if ((n%2)==(k%2))
		return (n==k ? 1/(long double)(2*n+1) : 0);
	LargeFactorHelper out;
	out.MultiplyBy(0);
	LargeFactorHelper Coef;
	for (unsigned int f = (0==n%2 ? n/2 : 1+ n/2); f<=n; ++f)
		for (unsigned int l = (0==k%2 ? k/2 : 1+ k/2); l<=k; ++l) {
			Coef.Clear();
			Coef.MultiplyByFactorial(2*f);
			Coef.MultiplyByFactorial(2*l);
			Coef.DivideByFactorial(2*f-n);
			Coef.DivideByFactorial(2*l-k);
			Coef.DivideByFactorial(f);
			Coef.DivideByFactorial(n-f);
			Coef.DivideByFactorial(l);
			Coef.DivideByFactorial(k-l);
			Coef.DivideBy(2*f-n+2*l-k+1);
			out += Coef;
		}
	out.DivideBy(2, n+k);
	return out.Output();
}

//Integral of Pn(x)Pk(x)(1-x)dx from -1 to 0
long double Int_PlPl_transf_1_0 (unsigned int n, unsigned int k)
{
	LargeFactorHelper out;
	out.MultiplyBy(0);
	LargeFactorHelper Coef;
	for (unsigned int f = (0==n%2 ? n/2 : 1+ n/2); f<=n; ++f)
		for (unsigned int l = (0==k%2 ? k/2 : 1+ k/2); l<=k; ++l) {
			Coef.Clear();
			Coef.MultiplyByFactorial(2*f);
			Coef.MultiplyByFactorial(2*l);
			Coef.DivideByFactorial(2*f-n);
			Coef.DivideByFactorial(2*l-k);
			Coef.DivideByFactorial(f);
			Coef.DivideByFactorial(n-f);
			Coef.DivideByFactorial(l);
			Coef.DivideByFactorial(k-l);
			Coef.DivideBy(2*f-n+2*l-k+1);
			Coef.DivideBy(2*f-n+2*l-k+2);
			Coef.MultiplyBy(4*f-2*n+4*l-2*k+3);
			if (0==((2*f-n+2*l-k)%2))
				out+=Coef;
			else
				out-=Coef;
		}
	out.DivideBy(2, n+k);
	return out.Output();
}
//Integral of Pn(x)Pk(x)(1-x)dx from 0 to 1
long double Int_PlPl_transf_0_1 (unsigned int n, unsigned int k)
{
	LargeFactorHelper out;
	out.MultiplyBy(0);
	LargeFactorHelper Coef;
	for (unsigned int f = (0==n%2 ? n/2 : 1+ n/2); f<=n; ++f)
		for (unsigned int l = (0==k%2 ? k/2 : 1+ k/2); l<=k; ++l) {
			Coef.Clear();
			Coef.MultiplyByFactorial(2*f);
			Coef.MultiplyByFactorial(2*l);
			Coef.DivideByFactorial(2*f-n);
			Coef.DivideByFactorial(2*l-k);
			Coef.DivideByFactorial(f);
			Coef.DivideByFactorial(n-f);
			Coef.DivideByFactorial(l);
			Coef.DivideByFactorial(k-l);
			Coef.DivideBy(2*f-n+2*l-k+1);
			Coef.DivideBy(2*f-n+2*l-k+2);
			out += Coef;
		}
	out.DivideBy(2, n+k);
	return out.Output();
}

long double Int_PlPl (unsigned int n, unsigned int k, long double from, long double to, long double dx)
{
	if (from>to) {
		long double mem = from;
		from = to;
		to = mem;
	}
	if (dx>(to-from))
		dx= (to-from)/dx;
	long double out = 0;
	long double x = from;
	LegendrePolynom P;
	for (unsigned int ix = 0, i_end = 1 + (unsigned int)((to-from)/dx); ix != i_end; ++ix) {
		out += ((x+dx)>to ? to - x : dx) * P(x,n)* P(x,k);
		x+=dx;
	}
	return out;
}

long double Int_PlPl_transf (unsigned int n, unsigned int k, long double from, long double to, long double dx)
{
	if (from>to) {
		long double mem = from;
		from = to;
		to = mem;
	}
	if (dx>(to-from))
		dx= (to-from)/dx;
	long double out = 0;
	long double x = from;
	LegendrePolynom P;
	for (unsigned int ix = 0, i_end = 1 + (unsigned int)((to-from)/dx); ix != i_end; ++ix) {
		out += ((x+dx)>to ? to - x : dx) * P(x,n)* P(x,k)*(1-x);
		x+=dx;
	}
	return out;
}
