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

void LargeFactorHelper::SimplifyFactorial(void)
{
	while (true) {
		unsigned int divider_max_factor = 1;
		unsigned int divider_max_index = 0;
		unsigned int denominator_max_factor = 1;
		unsigned int denominator_max_index = 0;
		for (unsigned int i = 0, i_end = factor_dividers.size(); i<i_end; ++i) {
			divider_max_factor = std::max(divider_max_factor, factor_dividers[i]);
			divider_max_index = i;
		}
		for (unsigned int i = 0, i_end = factor_dividers.size(); i<i_end; ++i) {
			denominator_max_factor = std::max(denominator_max_factor, factor_denominators[i]);
			denominator_max_index = i;
		}
		if ((1!=divider_max_factor)&&(1!=denominator_max_factor)) {
			factor_dividers.erase(factor_dividers.begin() + divider_max_index);
			factor_denominators.erase(factor_denominators.begin() + denominator_max_index);
			for (unsigned int n = divider_max_factor; n>denominator_max_factor;--n)
				MultiplyBy(n);
			for (unsigned int n = denominator_max_factor; n>divider_max_factor;--n)
				DivideBy(n);
			continue;
		}
		for (unsigned int n = divider_max_factor; n>1;--n)
			MultiplyBy(n);
		for (unsigned int n = denominator_max_factor; n>1;--n)
			DivideBy(n);
		return;
	}
}

/*void LargeFactorHelper::Factorize(void)
{
	for (unsigned int i = 0, i_end = dividers.size(); i!=i_end; ++i)
		for (unsigned int pr = 0, pr_end = prime_list.size(); pr!=pr_end; ++pr) {
			while (dividers[i]%prime_list[pr]==0) {
				++prime_dividers[pr];
				dividers[i]/=prime_list[pr];
			}
			if (1==dividers[i])
				break;
		}
	for (unsigned int i = 0, i_end = denominators.size(); i!=i_end; ++i)
		for (unsigned int pr = 0, pr_end = prime_list.size(); pr!=pr_end; ++pr) {
			while (denominators[i]%prime_list[pr]==0) {
				++prime_denominators[pr];
				denominators[i]/=prime_list[pr];
			}
			if (1==denominators[i])
				break;
		}
	denominators.erase(denominators.begin(),denominators.end());
	dividers.erase(dividers.begin(),dividers.end());
}*/

void LargeFactorHelper::Simplify(void)
{
	for (unsigned int pr = 0, pr_end = prime_list.size(); pr!=pr_end; ++pr) {
		if (prime_denominators[pr]>=prime_dividers[pr]){
			prime_denominators[pr] -= prime_dividers[pr];
			prime_dividers[pr] = 0;
		} else {
			prime_dividers[pr] -= prime_denominators[pr];
			prime_denominators[pr] = 0;
		}
	}
}

LargeFactorHelper::LargeFactorHelper(unsigned int N_max)
{
	n_max =	N_max > 2 ? N_max : 2;
	for (unsigned int i = 2; i<=n_max; ++i) {
		bool prime = true;
		for (unsigned int pr = 0; pr<prime_list.size(); ++pr) {
			if (0==i%prime_list[pr]) {
				prime = false;
				break;
			}
		}
		if (prime)
			prime_list.push_back(i);
	}
}

void LargeFactorHelper::MultiplyByFactorial(unsigned int n)
{
	if (n>n_max) {
		std::cout<<"LargeFactorHelper::MultiplyByFactorial() error: "<<n<<" is larger than n_max = "<<n_max<<std::endl;
		return;
	}
	factor_dividers.push_back(n);
}

void LargeFactorHelper::DivideByFactorial(unsigned int n)
{
	if (n>n_max) {
		std::cout<<"LargeFactorHelper::DivideByFactorial() error: "<<n<<" is larger than n_max = "<<n_max<<std::endl;
		return;
	}
	factor_denominators.push_back(n);
}

void LargeFactorHelper::MultiplyBy(unsigned int n, unsigned int power)
{
	if (power<1)
		return;
	if (n>n_max) {
		std::cout<<"LargeFactorHelper::MultiplyBy() error: "<<n<<" is larger than n_max = "<<n_max<<std::endl;
		return;
	}
	for (unsigned int pr = 0, pr_end = prime_list.size(); pr!=pr_end; ++pr) {
		while (n%prime_list[pr]==0) {
			prime_dividers[pr]+=power;
			n/=prime_list[pr];
		}
		if (1==n)
			break;
	}
}

void LargeFactorHelper::DivideBy(unsigned int n, unsigned int power)
{
	if (power<1)
		return;
	if (n>n_max) {
		std::cout<<"LargeFactorHelper::MultiplyBy() error: "<<n<<" is larger than n_max = "<<n_max<<std::endl;
		return;
	}
	for (unsigned int pr = 0, pr_end = prime_list.size(); pr!=pr_end; ++pr) {
		while (n%prime_list[pr]==0) {
			prime_denominators[pr]+=power;
			n/=prime_list[pr];
		}
		if (1==n)
			break;
	}
}

long double LargeFactorHelper::Output(void);

void LargeFactorHelper::Print(void);

void LargeFactorHelper::Clear(void);


long double argon_cross_elastic_diff (long double E, long double theta) {
	if (E <= 1) {

	}
	return 0;
}


