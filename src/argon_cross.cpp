#include "argon_cross.h"

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
			for (int i=(l_last+1);i<=l;++i) {
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
	for (int i=2;i<=l;++i) {
		mem = ((2*i-1)*x*P_last-(i-1)*P_last1)/i;
		P_last1 = P_last;
		P_last = mem;
		l_last1 = l_last; //==i-1
		l_last = i;
	}
	return P_last;
}

long double argon_cross_elastic_diff (long double E, long double theta) {
	if (E <= 1) {

	}
	return 0;
}


