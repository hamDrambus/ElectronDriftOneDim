#include <math>
#include <iostream>
#include <fstream>

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

//E in eV, theta in radians, output is in cm-17
long double argon_cross_elastic_diff (long double E, long double theta);
