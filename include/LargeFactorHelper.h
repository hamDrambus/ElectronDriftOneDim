#ifndef LARGE_FACTOR_HELPER_H
#define LARGE_FACTOR_HELPER_H

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "global_definitions.h"

class LargeFactorHelper //assumed that only numbers less than N are in the input
{
protected:
    unsigned long int n_max;
    static unsigned long int global_n_max;
    void SimplifyFactorial(void);
	void Simplify(void);
	std::vector<unsigned long int> prime_dividers; //according to the prime_list
	std::vector<unsigned int> prime_denominators; //according to the prime_list
	static std::vector<unsigned long int> g_prime_list; //all prime numbers <= n_max, starting from 2
	std::vector<unsigned long int> factor_dividers;
	std::vector<unsigned long int> factor_denominators;
	void ExtendPrimeList(unsigned long int N_max);
	void MultiplyBy0 (void);
	enum sign {plus = 1, zero = 0, minus = -1} sign_;
public:
	LargeFactorHelper(void);
	void MultiplyByFactorial(unsigned long int n);
	void DivideByFactorial(unsigned long int n);
	void MultiplyBy(unsigned long int n, unsigned int power = 1);
	void DivideBy(unsigned long int n, unsigned int power = 1);
	long double Output(void);
	void Print(void);
	void Clear(void);

	LargeFactorHelper& operator+= (const LargeFactorHelper& added);
	LargeFactorHelper& operator-= (const LargeFactorHelper& added);
};

#endif
