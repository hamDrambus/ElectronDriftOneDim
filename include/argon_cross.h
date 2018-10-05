#include <math>
#include <iostream>
#include <fstream>
#include <vector>

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

class LargeFactorHelper //assumed that only numbers less than N are in the input
{
protected:
    unsigned int n_max;
    void SimplifyFactorial(void);
	//void Factorize(void);
	void Simplify(void);
	std::vector<unsigned int> dividers;
	std::vector<unsigned int> denominators;
	std::vector<unsigned int> prime_dividers; //according to the prime_list
	std::vector<unsigned int> prime_denominators; //according to the prime_list
	std::vector<unsigned int> prime_list; //all prime numbers <= n_max, starting from 2
	std::vector<unsigned int> factor_dividers;
	std::vector<unsigned int> factor_denominators;
public:
	LargeFactorHelper(unsigned int N_max);
	void MultiplyByFactorial(unsigned int n);
	void DivideByFactorial(unsigned int n);
	void MultiplyBy(unsigned int n, unsigned int power = 1);
	void DivideBy(unsigned int n, unsigned int power = 1);
	long double Output(void);
	void Print(void);
	void Clear(void);
};

//E in eV, theta in radians, output is in cm-17
long double argon_cross_elastic_diff (long double E, long double theta);
