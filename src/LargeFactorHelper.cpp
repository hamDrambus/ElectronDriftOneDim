#include "LargeFactorHelper.h"

void LargeFactorHelper::SimplifyFactorial(void)
{
	while (true) {
		unsigned int divider_max_factor = 1;
		unsigned int divider_max_index = 0;
		unsigned int denominator_max_factor = 1;
		unsigned int denominator_max_index = 0;
		for (unsigned int i = 0, i_end = factor_dividers.size(); i<i_end; ++i) {
			if (divider_max_factor<factor_dividers[i]) {
				divider_max_index = i;
				divider_max_factor = factor_dividers[i];
			}
		}
		for (unsigned int i = 0, i_end = factor_denominators.size(); i<i_end; ++i) {
			if (denominator_max_factor<factor_denominators[i]) {
				denominator_max_index = i;
				denominator_max_factor = factor_denominators[i];
			}
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
		//no more pairs of factorials:
		for (unsigned int i = 0, i_end = factor_dividers.size(); i<i_end; ++i) {
			divider_max_factor = factor_dividers[i];
			for (unsigned int n = divider_max_factor; n>1;--n)
				MultiplyBy(n);
		}
		for (unsigned int i = 0, i_end = factor_denominators.size(); i<i_end; ++i) {
			denominator_max_factor = factor_denominators[i];
			for (unsigned int n = denominator_max_factor; n>1;--n)
				DivideBy(n);
		}
		return;
	}
}

void LargeFactorHelper::Simplify(void)
{
	for (unsigned int pr = 0, pr_end = prime_denominators.size(); pr!=pr_end; ++pr) {
		if (prime_denominators[pr]>=prime_dividers[pr]){
			prime_denominators[pr] -= prime_dividers[pr];
			prime_dividers[pr] = 0;
		} else {
			prime_dividers[pr] -= prime_denominators[pr];
			prime_denominators[pr] = 0;
		}
	}
}

unsigned long int LargeFactorHelper::global_n_max = 1;
std::vector<unsigned long int> LargeFactorHelper::g_prime_list;
LargeFactorHelper::LargeFactorHelper():n_max(1), sign_(plus)
{}

void LargeFactorHelper::ExtendPrimeList(unsigned long int N_max)
{
	if (N_max>global_n_max) {
		for (unsigned long int i = global_n_max+1; i<=N_max; ++i) {
			bool prime = true;
			for (unsigned int pr = 0; pr<g_prime_list.size(); ++pr) {
				if (0==i%g_prime_list[pr]) {
					prime = false;
					break;
				}
			}
			if (prime)
				g_prime_list.push_back(i);
		}
		global_n_max = N_max;
	}
	if (N_max > n_max) {
		unsigned long int size = 0, size_end = g_prime_list.size();
		for (size =0; size<size_end; ++size) {//
			if (N_max==g_prime_list[size]) {
				prime_denominators.resize(size+1,0);//only extends
				prime_dividers.resize(size+1,0);
				n_max = N_max;
				return;
			}
			if (N_max<g_prime_list[size]) {
				prime_denominators.resize(size,0);//only extends
				prime_dividers.resize(size,0);
				n_max = N_max;
				return;
			}
		}
		prime_denominators.resize(size_end,0);//only extends
		prime_dividers.resize(size_end,0);
		n_max = N_max;
	}
}

void LargeFactorHelper::MultiplyByFactorial(unsigned long int n)
{
	if (sign_==zero)
		return;
	ExtendPrimeList(n);
	factor_dividers.push_back(n);
}

void LargeFactorHelper::DivideByFactorial(unsigned long int n)
{
	if (sign_==zero)
		return;
	ExtendPrimeList(n);
	factor_denominators.push_back(n);
}

void LargeFactorHelper::MultiplyBy(unsigned long int n, unsigned int power)
{
	if ((power<1)||(sign_==zero))
		return;
	ExtendPrimeList(n);
	if (0==n) {
		MultiplyBy0();
		return;
	}
	for (unsigned int pr = 0, pr_end = prime_dividers.size(); pr!=pr_end; ++pr) {
		while (n%g_prime_list[pr]==0) {
			prime_dividers[pr]+=power;
			n/=g_prime_list[pr];
		}
		if (1==n)
			break;
	}
	if (1!=n)
		std::cout<<"LargeFactorHelper::MultiplyBy(): Error: number was not fully factorized!"<<std::endl;
}

void LargeFactorHelper::MultiplyBy0 (void)
{
	Clear();
	sign_ = zero;
}

void LargeFactorHelper::DivideBy(unsigned long int n, unsigned int power)
{
	if ((power<1)||(sign_==zero))
		return;
	ExtendPrimeList(n);
	if (0==n)
		return;
	for (unsigned int pr = 0, pr_end = prime_denominators.size(); pr!=pr_end; ++pr) {
		while (n%g_prime_list[pr]==0) {
			prime_denominators[pr]+=power;
			n/=g_prime_list[pr];
		}
		if (1==n)
			break;
	}
	if (1!=n)
		std::cout<<"LargeFactorHelper::DivideBy(): Error: number was not fully factorized!"<<std::endl;
}

long double LargeFactorHelper::Output(void)
{
	SimplifyFactorial();
	Simplify();
	long double out = 1;
	for (unsigned long int i=0, end_ = prime_dividers.size(); i!=end_; ++i) {
		for (unsigned int k = 1; k<=prime_dividers[i];++k) {
			out*=g_prime_list[i];
		}
		for (unsigned int k = 1; k<=prime_denominators[i];++k) {
			out/=g_prime_list[i];
		}
	}
	return out*sign_;
}

LargeFactorHelper& LargeFactorHelper::operator+= (const LargeFactorHelper& added)
{
	if (added.sign_==0)
		return *this;
	if (this->sign_==0) {
		(*this) = added;
		return *this;
	}
	SimplifyFactorial();
	Simplify();
	LargeFactorHelper copy = added;
	copy.SimplifyFactorial();
	copy.Simplify();
	unsigned long int left_size = this->prime_denominators.size();
	unsigned long int right_size = copy.prime_denominators.size();
	unsigned long int max_size = std::max(left_size, right_size);
	//3^2 5    3 5 11
	//----- + --------
	//2^2 7   2^2 13^2
	for (unsigned long int i=0; i<max_size; ++i) {
		if (i>=left_size) {
			this->DivideBy(g_prime_list[i], copy.prime_denominators[i]);
			this->MultiplyBy(g_prime_list[i], copy.prime_denominators[i]);
			continue;
		}
		if (i>=right_size) {
			copy.MultiplyBy(g_prime_list[i], this->prime_denominators[i]);
			continue;
		}
		if (copy.prime_denominators[i]>this->prime_denominators[i]) {
			unsigned int factor = copy.prime_denominators[i] - this->prime_denominators[i];
			this->DivideBy(g_prime_list[i], factor);
			this->MultiplyBy(g_prime_list[i], factor);
		}
		else
			copy.MultiplyBy(g_prime_list[i], this->prime_denominators[i] - copy.prime_denominators[i]);
	}
	//3^2 5 13^2     3 5 7 11
	//----------  + ----------
	//2^2 7 13^2    2^2 7 13^2
	left_size = this->prime_dividers.size();
	right_size = copy.prime_dividers.size();
	max_size = std::max(left_size, right_size);
	unsigned long int left_factor = 1, right_factor = 1;
	for (unsigned long int i=0; i<max_size; ++i) { //excluding common multiplier from two dividers
		if (i>=left_size) {
			for (unsigned int k=1; k<=copy.prime_dividers[i]; ++k)
				right_factor *=g_prime_list[i];
			continue;
		}
		if (i>=right_size) {
			for (unsigned int k=1; k<=this->prime_dividers[i]; ++k)
				left_factor *=g_prime_list[i];
			this->prime_dividers[i] = 0; //== common_factor[i]
			continue;
		}
		if (copy.prime_dividers[i]>this->prime_dividers[i]) {
			for (unsigned int k=1; k<=(copy.prime_dividers[i] - this->prime_dividers[i]); ++k)
				right_factor *=g_prime_list[i];
			//this->prime_dividers[i] == common_factor[i];
		} else {
			for (unsigned int k=1; k<=(this->prime_dividers[i] - copy.prime_dividers[i]); ++k)
				left_factor *=g_prime_list[i];
			this->prime_dividers[i] = copy.prime_dividers[i]; //== common_factor[i]
		}
	}
	//3^2 5 13^2    3 5 7 11 	  3 5 (3 13^2)   3 5 (7 11)      3 5          3 5 (7 11)
	//---------- + ---------- =   -----------  + ---------- --> ----------- + ----------
	//2^2 7 13^2   2^2 7 13^2      2^2 7 13^2    2^2 7 13^2     2^2 7 13^2    2^2 7 13^2
	//then left ratio (this->) is multiplied by (left_factor + right_factor) = (3 13^2 + 7 11)
	unsigned long int sum =0;
	if (this->sign_==copy.sign_) { //the following code caused by the usage of unsigned integers
		sum = left_factor+right_factor;
	} else {
		if (this->sign_==plus) {
			if (left_factor<right_factor) {
				sum = right_factor - left_factor;
				this->sign_ = minus;
			} else
				sum = left_factor - right_factor;
		} else {
			if (left_factor<right_factor) {
				sum = right_factor - left_factor;
				this->sign_ = plus;
			} else
				sum = left_factor - right_factor;
		}
	}
	this->MultiplyBy(sum); //excluded "left factor" from the divider in the previous loop (set it to "common factor")
	return *this;
}

LargeFactorHelper& LargeFactorHelper::operator-= (const LargeFactorHelper& added)
{
	if (added.sign_==0)
		return *this;
	if (this->sign_==0) {
		(*this) = added;
		if (this->sign_ == plus)
			this->sign_ = minus;
		else
			this->sign_ = plus;
		return *this;
	}
	// A - B == -(-A + B). (can't change sign of added since it is 'const', so sign of this is changed instead)
	if (this->sign_ == plus)
		this->sign_ = minus;
	else
		this->sign_ = plus;
	(*this)+=added;
	if (this->sign_!=zero) {
		if (this->sign_ == plus)
			this->sign_ = minus;
		else
			this->sign_ = plus;
	}
	return *this;
}

void LargeFactorHelper::Print(void)
{
	if (sign_==zero) {
		std::cout<<"0"<<std::endl;
		return;
	}
	if (sign_==plus)
		std::cout<<"  ";
	else
		std::cout<<"   ";
	int char_count_up = 0;
	int char_count_down = 0;
	for (int i = 0, end_ = factor_dividers.size(); i!=end_; ++i) {
		std::cout<<factor_dividers[i]<<"! ";
		char_count_up +=3 + (int) log10(factor_dividers[i]);
	}
	for (int i = 0, end_ = prime_dividers.size(); i!=end_; ++i) {
		if (0!=prime_dividers[i]) {
			std::cout<<g_prime_list[i]<<"^"<<prime_dividers[i]<<" ";
			char_count_up +=4 + (int) log10(prime_dividers[i]) + (int) log10(g_prime_list[i]);
		}
	}
	for (int i = 0, end_ = factor_denominators.size(); i!=end_; ++i)
		char_count_down +=3 + (int) log10(factor_denominators[i]);
	for (int i = 0, end_ = prime_denominators.size(); i!=end_; ++i)
		if (0!= prime_denominators[i])
			char_count_down +=4 + (int) log10(prime_denominators[i]) + (int) log10(g_prime_list[i]);
	std::cout<<std::endl;
	if (sign_==plus)
		std::cout<<"= ";
	else
		std::cout<<"= _";
	for (int i=0; i<std::max(char_count_up, char_count_down); ++i)
		std::cout<<"-";
	std::cout<<std::endl;
	if (sign_==plus)
		std::cout<<"  ";
	else
		std::cout<<"   ";
	for (int i = 0, end_ = factor_denominators.size(); i!=end_; ++i)
		std::cout<<factor_denominators[i]<<"! ";
	for (int i = 0, end_ = prime_denominators.size(); i!=end_; ++i)
		if (0!= prime_denominators[i])
			std::cout<<g_prime_list[i]<<"^"<<prime_denominators[i]<<" ";
	std::cout<<std::endl;
}

void LargeFactorHelper::Clear(void)
{
	for (unsigned int i=0, end_ = prime_denominators.size(); i!=end_; ++i) {
		prime_denominators[i]=0;
		prime_dividers[i]=0;
	}
	factor_denominators.erase(factor_denominators.begin(), factor_denominators.end());
	factor_dividers.erase(factor_dividers.begin(), factor_dividers.end());
	sign_=plus;
}
