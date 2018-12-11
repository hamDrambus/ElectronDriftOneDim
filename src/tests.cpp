#include "tests.h"

void test_polynomial_fit (void)
{
	std::string fname_raw = "tests/sinx_to_x_raw.txt";
	std::string fname_fit1 = "tests/sinx_to_x_fit1.txt";
	std::string fname_fit2 = "tests/sinx_to_x_fit2.txt";
	std::string fname_fit3 = "tests/sinx_to_x_fit3.txt";
	std::vector<double> xs, ys;
	xs.resize(50);
	ys.resize(50);
	std::ofstream str;
	str.open(fname_raw,std::ios_base::trunc);
	for (int i=0; i!=50;++i) { //not sorted
		double x = i*M_PI*3/100 + 0.02;
		xs[ i>=25 ? i-25 : i+25 ] = x;
		double y = sin(x)/x;
		ys[ i>=25 ? i-25 : i+25 ] = y;
		str<<x<<"\t"<<y<<std::endl;
	}
	str.close();
	std::vector<double> ys1 = ys;
	ys1.pop_back();
	//TESTING PARAMETERS:
	DataVector data1(3, 4);
	data1.initialize(xs, ys1, 3, 4);
	std::cout<<"data(0.5)(xs, ys1, 3, 4): "<<data1(0.5)<<std::endl;
	data1.initialize(xs, ys, 3, 3);
	std::cout<<"data(0.5)(xs, ys, 3, 3): "<<data1(0.5)<<std::endl;
	data1.setNused(4);
	std::cout<<"data(0.5)(xs, ys, 3, 4): "<<data1(0.5)<<std::endl;
	data1.initialize(xs, ys, 50, 51);
	std::cout<<"data(0.5)(xs, ys, 50, 51): "<<data1(0.5)<<std::endl;
	data1.initialize(xs, ys, 49, 50);
	std::cout<<"data(0.5)(xs, ys, 49, 50): "<<data1(0.5)<<std::endl;
	data1.setOrder(0);
	data1.setNused(1);
	std::cout<<"data(0.5)(xs, ys, 0, 1): "<<data1(0.5)<<std::endl;
	data1.setOrder(0);
	data1.setNused(50);
	std::cout<<"data(0.5)(xs, ys, 0, 50): "<<data1(0.5)<<std::endl;
	data1.setOrder(0);
	data1.setNused(51);
	std::cout<<"data(0.5)(xs, ys, 0, 51): "<<data1(0.5)<<std::endl;

	//TESTING QUALITY:
	data1.setOrder(4);
	data1.setNused(5);
	str.open(fname_fit1, std::ios_base::trunc);
	for (int i = 0; i<150; ++i) {
		double x = i*M_PI*3/250;
		double y = data1(x, x);
		str<<x<<"\t"<<y<<std::endl;
	}
	str.close();
	data1.setOrder(6);
	data1.setNused(15);
	data1.set_out_value(0.5);
	str.open(fname_fit2, std::ios_base::trunc);
	for (int i = 0; i<150; ++i) {
		double x = i*M_PI*3/250;
		double y = data1(x, x);
		str<<x<<"\t"<<y<<std::endl;
	}
	str.close();
	str.close();
	data1.setOrder(2);
	data1.setNused(10);
	data1.unset_out_value();
	str.open(fname_fit3, std::ios_base::trunc);
	for (int i = 0; i<150; ++i) {
		double x = i*M_PI*3/250;
		double y = data1(x, x);
		str<<x<<"\t"<<y<<std::endl;
	}
	str.close();
	std::string name = "tests/test_fit.sc";
	str.open(name, std::ios_base::trunc);
	str<<"plot \"tests/sinx_to_x_raw.txt\" u 1:2 title \"raw sin(x)/x\""<<std::endl;
	str<<"replot \"tests/sinx_to_x_fit1.txt\" u 1:2 title \"fit1\""<<std::endl;
	str<<"replot \"tests/sinx_to_x_fit2.txt\" u 1:2 title \"fit2\""<<std::endl;
	str<<"replot \"tests/sinx_to_x_fit3.txt\" u 1:2 title \"fit3\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_phase_shift_fit (void)
{
	std::string fname_McEachran = "tests/phase_shifts_McEachran_";
	std::string fname_MERT = "tests/phase_shifts_MERT.txt";
	std::string fname_phase_fit = "tests/phase_shifts_fit_exp.txt";

	std::ofstream str;
	for (unsigned int l=0; l<ArExper.phase_shifts_.size();++l) {
		std::string fname = fname_McEachran + std::to_string(l) + ".txt";
		str.open(fname, std::ios_base::trunc);
		str<<"E[eV]\tphase shift "<<l<<std::endl;
		for (std::size_t i = 0, i_end = ArExper.phase_shifts_[l].size(); i!=i_end; ++i)
			str<< pow(ArExper.phase_shifts_[l].getX(i)/a_h_bar_2e_m_e_SIconst, 2)<<"\t"<<ArExper.phase_shifts_[l].getY(i)<<std::endl;
		str.close();
	}

	int err;
	double E, k;
	{
		str.open(fname_phase_fit, std::ios_base::trunc);
		str<<"E[eV]\tphase shifts 0, 1, ... "<<std::endl;
		EnergyScanner EnRange(EnergyScanner::PlotDiffXS);
		while (true) {
			E = EnRange.Next(err);
			if (0!=err)
				break;
			k = sqrt(E)*a_h_bar_2e_m_e_SIconst;
			str<<E<<"\t";
			for (std::size_t l = 0, l_end = ArExper.phase_shifts_.size(); l!=l_end; ++l)
				str<<ArExper.phase_shifts_[l](k,k)<<"\t";
			str<<std::endl;
		}
		str.close();
	}

	{
		str.open(fname_MERT, std::ios_base::trunc);
		str<<"E[eV]\tphase shifts 0, 1, ... "<<std::endl;
		EnergyScanner EnRange(EnergyScanner::PlotDiffXS);
		while (true) {
			E = EnRange.Next(err);
			if (0!=err)
				break;
			if (E>THRESH_E_XS_)
				break;
			k = sqrt(E)*a_h_bar_2e_m_e_SIconst;
			str<<E<<"\t";
			for (std::size_t l = 0, l_end = ArExper.phase_shifts_.size(); l!=l_end; ++l) {
				long double tan, sin, cos;
				argon_phase_values_MERT5(k, l, tan, sin, cos);
				tan = std::atan(tan);
				str<<tan<<"\t";
			}
			str<<std::endl;
		}
		str.close();
	}

	for (std::size_t l = 0, l_end = ArExper.phase_shifts_.size(); l!=l_end; ++l) {
		std::string name = std::string("tests/test_phase_shift_") + std::to_string(l) + ".sc";
		str.open(name, std::ios_base::trunc);
		str<<"set logscale x"<<std::endl;
		str<<"plot '"<<fname_McEachran + std::to_string(l) + ".txt" <<"' u 1:2 title 'McEachran_"<<l<<"'"<<std::endl;
		str<<"replot '"<<fname_MERT <<"' u 1:"<<2+l<<" w line lc rgb \"#FF0000\" title 'MERT_"<<l<<"'"<<std::endl;
		str<<"replot '"<<fname_phase_fit <<"' u 1:"<<2+l<<" w line lc rgb \"#000000\" title 'fit_"<<l<<"'"<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}
}

void test_legendre_polynomial(void)
{
	LegendrePolynom Pl;
	std::cout<<"Pl(3,0.36) =\t "<<Pl(0.36, 3)<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.42336"<<std::endl;
	std::cout<<"Pl(4,0.36) =\t "<<Pl(0.36, 4)<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.0375168"<<std::endl;
	std::cout<<"Pl(20,0.36) =\t "<<Pl(0.36, 20)<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.0542800664"<<std::endl;
	std::cout<<"Pl(50,0.5) =\t "<<Pl(0.5, 50)<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.031059099239"<<std::endl;
	std::cout<<"Pl(49,0.5) =\t "<<Pl(0.5, 49)<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.086292778960940"<<std::endl;
	std::cout<<"Pl(48,0.5) =\t "<<Pl(0.5, 48)<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.118866275929531"<<std::endl;
	std::cout<<"Pl(47,0.5) =\t "<<Pl(0.5, 47)<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.032013921114504"<<std::endl;
	std::cout<<"Pl(6,-1) =\t "<<Pl(-1.0, 6)<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	std::cout<<"Pl(6,0) =\t "<<Pl(0, 6)<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.3125"<<std::endl;
	std::cout<<"Pl(6,1) =\t "<<Pl(1, 6)<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	std::cout<<"Pl(7,-1) =\t "<<Pl(-1.0, 7)<<std::endl;
	std::cout<<"Wolfram alpha:\t -1"<<std::endl;
	std::cout<<"Pl(7,0) =\t "<<Pl(0, 7)<<std::endl;
	std::cout<<"Wolfram alpha:\t 0"<<std::endl;
	std::cout<<"Pl(7,1) =\t "<<Pl(1, 7)<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
}

void test_legendre_intregral (void)
{
	int Ncalls = 10;
	auto start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_1_0 (20, 13);
	}*/
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 13, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 13, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_high_dx = end-start;
	/*start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 13, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(13,x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_1_0 (20, 13)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<Int_PlPl (20, 13, -1, 0, 1e-5)<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<Int_PlPl (20, 13, -1, 0, 1e-7)<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (20, 13, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.00217109"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_1_0 (20, 18);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 18, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 18, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 18, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(18,x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_1_0 (20, 18)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<Int_PlPl (20, 18, -1, 0, 1e-5)<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<Int_PlPl (20, 18, -1, 0, 1e-7)<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (20, 18, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.0"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_1_0 (18, 18);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (18, 18, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (18, 18, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (18, 18, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(18,x)Pl(18,x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_1_0 (18, 18)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<Int_PlPl (18, 18, -1, 0, 1e-5)<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<Int_PlPl (18, 18, -1, 0, 1e-7)<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (18, 18, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.027027"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_0_1 (20, 15);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 15, 0, 1, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 15, 0, 1, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 15, 0, 1, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(15,x)dx [0.0 ; 1.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_0_1 (20, 15)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<Int_PlPl (20, 15, 0, 1, 1e-5)<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<Int_PlPl (20, 15, 0, 1, 1e-7)<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (20, 15, 0, 1, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.00307571"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf_1_0 (20, 15);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(15,x)(1-x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_transf_1_0 (20, 15)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<Int_PlPl_transf (20, 15, -1, 0, 1e-5)<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<Int_PlPl_transf (20, 15, -1, 0, 1e-7)<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl_transf (20, 15, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: -0.00307571"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf_0_1 (20, 15);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, 0, 1, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, 0, 1, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, 0, 1, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(15,x)dx [0.0 ; 1.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_transf_0_1 (20, 15)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<Int_PlPl_transf (20, 15, 0, 1, 1e-5)<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<Int_PlPl_transf (20, 15, 0, 1, 1e-7)<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl_transf (20, 15, 0, 1, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.00307571"<<std::endl;
	std::cout<<"========================"<<std::endl;
}

void test_colored_interval (void)
{
	ColoredInterval ran1 (-1,0, 1);
	ColoredInterval ran2 (0,1, 1);
	ColoredInterval ran3 (1,4, 2);
	ColoredInterval ran4 (2,3, 0.5);
	ColoredInterval ran5 (2.8,3, 0.1);
	ColoredInterval ran6 (5,6, 0.1);
	ColoredInterval ran7 (-2,10, 5);
	ColoredRange range1, range2;
	std::cout<<"Sum of following intervals [from; color; to]:"<<std::endl;
	std::cout<<"[-1; 1; 0] + [0; 1; 1] + [2; 0.5; 3] + [-2; 5; 10]"<<std::endl;
	std::cout<<"Result: "<<std::endl;
	range1 = ran1 + ran2 + ran4 + ran7;
	range1.Print(std::cout);
	std::cout<<"========================"<<std::endl;
	std::cout<<"Sum of following intervals [from; color; to]:"<<std::endl;
	std::cout<<"[-1; 1; 0] + [0; 1; 1] + [2; 0.5; 3] + [-2; 5; 10] + "<<std::endl;
	std::cout<<"+ [1; 2; 4] + [2.8; 0.1; 3] + [5; 0.1; 6]"<<std::endl;
	std::cout<<"Result: "<<std::endl;
	range2 = range1;
	range2 = range2 + ran3 + ran5 + ran6;
	range2.Print(std::cout);
	std::cout<<"========================"<<std::endl;
	std::cout<<"Scan over following intervals [from; step; to]:"<<std::endl;
	std::cout<<"[0; 0.5; 1] + [1; 1; 3] + [6; 0.5; 7]"<<std::endl;
	ColoredInterval ran8 (0,1, 0.5);
	ColoredInterval ran9 (1,3, 1);
	ColoredInterval ran10 (6,7, 0.5);
	range1 = ran8 + ran9 +ran10;
	for (int i=0, end_=range1.NumOfIndices(); i!=end_; ++i) {
		std::cout<<range1.Value(i)<<"; ";
	}
	std::cout<<std::endl;
	std::cout<<"========================"<<std::endl;
}

void test_factor_helper (void)
{
	LargeFactorHelper helper;
	std::cout<<"396*500*22/2^5/11"<<std::endl;
	helper.MultiplyBy(396);
	helper.MultiplyBy(500);
	helper.MultiplyBy(22);
	helper.DivideBy(2, 5);
	helper.DivideBy(11);
	helper.Print();
	std::cout<<"Direct = ";
	long double val = 396.0*500*22/pow(2,5)/11;
	std::cout<<val<<std::endl;
	std::cout<<"Helper = "<< helper.Output()<<std::endl;
	std::cout<<"Wolfram = 12375"<<std::endl;
	helper.Clear();
	//=================================================
	std::cout<<"C(9,20)/80"<<std::endl;
	helper.MultiplyByFactorial(20);
	helper.DivideByFactorial(9);
	helper.DivideByFactorial(11);
	helper.DivideBy(80);
	helper.Print();
	std::cout<<"Direct = ";
	val = 1;
	for (int i=12;i<=20;++i)
		val*=i;
	for (int i=1;i<=9;++i)
		val/=i;
	val/=80;
	std::cout<<val<<std::endl;
	std::cout<<"Helper = "<< helper.Output()<<std::endl;
	std::cout<<"Wolfram = 2099.5"<<std::endl;
	helper.Clear();
	//=================================================
	std::cout<<"C(9,20)*C(20, 60)/80/3^6/10!"<<std::endl;
	helper.MultiplyByFactorial(20);
	helper.DivideByFactorial(9);
	helper.DivideByFactorial(11);
	helper.MultiplyByFactorial(60);
	helper.DivideByFactorial(20);
	helper.DivideByFactorial(40);
	helper.DivideBy(80);
	helper.DivideBy(3,6);
	helper.DivideByFactorial(10);
	helper.Print();
	std::cout<<"Direct = ";
	val = 1;
	for (int i=12;i<=20;++i)
		val*=i;
	for (int i=2;i<=9;++i)
		val/=i;
	for (int i=41;i<=60;++i)
		val*=i;
	for (int i=2;i<=20;++i)
		val/=i;
	val/=80;
	val/=pow(3,6);
	for (int i=2;i<=10;++i)
		val/=i;
	std::cout<<val<<std::endl;
	std::cout<<"Helper = "<< helper.Output()<<std::endl;
	std::cout<<"Wolfram = 3.32682902e+9"<<std::endl;
	helper.Clear();
	//=================================================
	LargeFactorHelper helper1;
	std::cout<<"30*27*17/125/121/8 + 3^4*13*5/17/11/1024:"<<std::endl;
	helper.MultiplyBy(30);
	helper.MultiplyBy(27);
	helper.MultiplyBy(17);
	helper.DivideBy(125);
	helper.DivideBy(121);
	helper.DivideBy(8);
	helper.Print();
	helper1.MultiplyBy(3,4);
	helper1.MultiplyBy(13);
	helper1.MultiplyBy(5);
	helper1.DivideBy(17);
	helper1.DivideBy(11);
	helper1.DivideBy(1024);
	helper1.Print();
	std::cout<<"Direct = ";
	val = 30.0*27*17/125/121/8 + 9.0*9*13*5/17/11/1024;
	std::cout<<val<<std::endl;
	helper+=helper1;
	std::cout<<"Helper = "<< helper.Output()<<std::endl;
	helper.Clear();
	helper1.Clear();
	std::cout<<"27/81 - 6/18:"<<std::endl;
	helper.MultiplyBy(27);
	helper.DivideBy(81);
	helper.Print();
	helper1.MultiplyBy(6);
	helper1.DivideBy(18);
	helper1.Print();
	std::cout<<"Direct = ";
	val = 27.0/81 - 6.0/18;
	std::cout<<val<<std::endl;
	helper-=helper1;
	std::cout<<"Helper = "<< helper.Output()<<std::endl;
	helper.Clear();
	helper1.Clear();
}

void test_diff_tot_cross (void)
{
	std::string fname_diff = "tests/diff_cross_elastic_10eV.txt";
	std::string fname_tot = "tests/total_elastic_from_diff.txt";
	std::ofstream str;
	str.open(fname_diff, std::ios_base::trunc);
	str<<"theta[deg]\tXS[1e-20m^2]"<<std::endl;
	for (int i=0; i<600; ++i) {
		double th = i*(M_PI)/599;
		str<<th*180/M_PI<<"\t"<<argon_cross_elastic_diff(10.0, th)<<std::endl;
	}
	str.close();

	str.open(fname_tot, std::ios_base::trunc);
	str<<"E[eV]\tXS from diff [1e-20m^2]\tXS tot [1e-20m^2]\tXS tot PS"<<std::endl;
	int err;
	EnergyScanner eScan(EnergyScanner::PlotDiffXS);
	while (true) {
		double E = eScan.Next(err);
		if (0!=err)
			break;
		long double integral = 0;
		for (int j=0;j<10001; ++j)
			integral+=(M_PI/10000.0)*argon_cross_elastic_diff(E, j*M_PI/10000.0)*sin(j*M_PI/10000.0);
		str<<E<<"\t"<<integral<<"\t"<<argon_cross_elastic(E)<<"\t"<<argon_cross_elastic_from_phases(E)<<std::endl;
	}
	str.close();

	std::string name = "tests/test_diff_XS_elastic.sc";
	str.open(name, std::ios_base::trunc);
	str<<"plot \""<<fname_diff<<"\" u 1:2 title \"Diff. XS\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
	name = "tests/test_diff_XS_elastic_by_total.sc";
	str.open(name, std::ios_base::trunc);
	str<<"set logscale x"<<std::endl;
	str<<"set logscale y"<<std::endl;
	str<<"plot \""<<fname_tot<<"\" u 1:2 title \"total from diff. XS\""<<std::endl;
	str<<"replot \"data/ArScatteringCross.dat\" u 1:2 lc rgb \"#000000\" title \"total XS experiment\""<<std::endl;
	str<<"replot \""<<fname_tot<<"\" u 1:4 title \"total XS from phase shifts\""<<std::endl;
	str<<"replot \""<<fname_tot<<"\" u 1:3 w lines title \"total XS\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_backward_scatter_prob (void)
{
	std::string fname = "tests/backward_scattering_prob.txt";
	std::ofstream str;
	int err;
	{
		str.open(fname, std::ios_base::trunc);
		str<<"E[eV]\tW_backward"<<std::endl;
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<argon_back_scatter_prob(E)<<std::endl;
		}
		str.close();
	}
	std::string name = "tests/test_backward_scatter.sc";
	str.open(name, std::ios_base::trunc);
	//str<<"set logscale x"<<std::endl;
	str<<"plot \""<<fname<<"\" u 1:2 title \"W_backward\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_TM_forward (void)
{
	std::string fname = "tests/TM_forward.txt";
	std::ofstream str;
	int err;
	{
		str.open(fname, std::ios_base::trunc);
		str<<"E[eV]\tTM_forward"<<std::endl;
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<argon_TM_forward(E)<<std::endl;
		}
		str.close();
	}
	std::string name = "tests/test_TM_forward.sc";
	str.open(name, std::ios_base::trunc);
	//str<<"set logscale x"<<std::endl;
	str<<"plot \""<<fname<<"\" u 1:2 title \"TM_forward\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_TM_backward (void)
{
	std::string fname = "tests/TM_backward.txt";
	std::ofstream str;
	int err;
	{
		str.open(fname, std::ios_base::trunc);
		str<<"E[eV]\tTM_backward"<<std::endl;
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<argon_TM_backward(E)<<std::endl;
		}
		str.close();
	}
	std::string name = "tests/test_TM_backward.sc";
	str.open(name, std::ios_base::trunc);
	//str<<"set logscale x"<<std::endl;
	str<<"plot \""<<fname<<"\" u 1:2 title \"TM_backward\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_total_cross_all (void)
{
	std::ofstream str, str1;

	std::string fname_XS = "tests/total_elastic_and_resonances_XS.txt";
	std::string fname_XS_ext = "tests/excitation_XS.txt";
	std::string fname_XS_ext2 = "tests/excitation_XS2.txt";
	int err;
	{
		EnergyScanner EnRange(EnergyScanner::PlotAllXS);
		str.open(fname_XS, std::ios_base::trunc);
		str<<std::scientific;
		str<<"E[eV]\tXS elastic+Feshbach resonances total [1e-20 m^2]"<<std::endl;
		str1.open(fname_XS_ext, std::ios_base::trunc);
		str1<<std::scientific;
		str1<<"E[eV]\tXS S ext.[1e-20 m^2]\tXS P ext.\tXS sum ext.\tXS ion."<<std::endl;
		while (true) {
			double E = EnRange.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<ArTables.XS_elastic(E)+ArTables.XS_resonance_3o2(E)+ArTables.XS_resonance_1o2(E)<<std::endl;
			double XS_S=0, XS_P=0, XS_EXT=0, XS_ION=0;
			for (int j =0, end_ = ArExper.ionizations.size(); j!=end_; ++j) {
				XS_ION+=ArExper.ionizations[j](E);
			}
			for (int j =0, end_ = ArExper.excitations.size(); j!=end_; ++j) {
				std::string name = ArExper.excitations[j].get_name();
				if (name.find("S")!=std::string::npos)
					XS_S+=ArExper.excitations[j](E);
				if (name.find("P")!=std::string::npos)
					XS_P+=ArExper.excitations[j](E);
				XS_EXT+=ArExper.excitations[j](E);
			}
			XS_S = std::max(1e-4, XS_S); //for logscale
			XS_P = std::max(1e-4, XS_P);
			XS_EXT = std::max(1e-4, XS_EXT);
			XS_ION = std::max(1e-4, XS_ION);
			str1<<E<<"\t"<<XS_S<<"\t"<<XS_P<<"\t"<<XS_EXT<<"\t"<<XS_ION<<std::endl;
		}
		str.close();
		str1.close();
	}

	{
		EnergyScanner EnRange(EnergyScanner::PlotInelasticXS);
		str.open(fname_XS_ext2, std::ios_base::trunc);
		str<<std::scientific;
		str<<"E[eV]\tXS S ext.[1e-20 m^2]\tXS P ext.\tXS sum ext.\tXS ion."<<std::endl;
		while (true) {
			double E = EnRange.Next(err);
			if (0!=err)
				break;
			double XS_S=0, XS_P=0, XS_EXT=0, XS_ION=0;
			for (int j =0, end_ = ArExper.ionizations.size(); j!=end_; ++j) {
				XS_ION+=ArExper.ionizations[j](E);
			}
			for (int j =0, end_ = ArExper.excitations.size(); j!=end_; ++j) {
				std::string name = ArExper.excitations[j].get_name();
				if (name.find("S")!=std::string::npos)
					XS_S+=ArExper.excitations[j](E);
				if (name.find("P")!=std::string::npos)
					XS_P+=ArExper.excitations[j](E);
				XS_EXT+=ArExper.excitations[j](E);
			}
			XS_S = std::max(1e-4, XS_S); //for logscale
			XS_P = std::max(1e-4, XS_P);
			XS_EXT = std::max(1e-4, XS_EXT);
			XS_ION = std::max(1e-4, XS_ION);
			str<<E<<"\t"<<XS_S<<"\t"<<XS_P<<"\t"<<XS_EXT<<"\t"<<XS_ION<<std::endl;
		}
		str.close();
	}

	std::string name = "tests/test_XS_all.sc";
	str.open(name, std::ios_base::trunc);
	str<<"set logscale x"<<std::endl;
	str<<"set logscale y"<<std::endl;
	str<<"plot \""<<fname_XS<<"\" u 1:2 w lines lc rgb \"#000000\" title \"XS elastic + Fishbach resonance\""<<std::endl;
	str<<"replot \""<<fname_XS_ext<<"\" u 1:2 w lines lc rgb \"#0000CC\" title \"XS S excitation\""<<std::endl;
	str<<"replot \""<<fname_XS_ext<<"\" u 1:3 w lines lc rgb \"#10CA73\" title \"XS P excitation\""<<std::endl;
	str<<"replot \""<<fname_XS_ext<<"\" u 1:4 w lines lc rgb \"#D90F2B\" title \"XS sum excitation\""<<std::endl;
	str<<"replot \""<<fname_XS_ext<<"\" u 1:5 w lines lc rgb \"#D90FD0\" title \"XS ionization\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);

	name = "tests/test_XS_ext.sc";
	str.open(name, std::ios_base::trunc);
	str<<"set logscale x"<<std::endl;
	str<<"set logscale y"<<std::endl;
	str<<"plot \""<<fname_XS_ext2<<"\" u 1:2 w lines lc rgb \"#0000CC\" title \"XS S excitation\""<<std::endl;
	str<<"replot \""<<fname_XS_ext2<<"\" u 1:3 w lines lc rgb \"#10CA73\" title \"XS P excitation\""<<std::endl;
	str<<"replot \""<<fname_XS_ext2<<"\" u 1:4 w lines lc rgb \"#D90F2B\" title \"XS sum excitation\""<<std::endl;
	str<<"replot \""<<fname_XS_ext2<<"\" u 1:5 w lines lc rgb \"#D90FD0\" title \"XS ionization\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);

}

//Dependent on previous tests.
//This function tests that interpolation/fit of ArDataTables works properly.
void test_data_table (void)
{
	int err;
	std::ofstream str;
	{
		EnergyScanner eScan(EnergyScanner::PlotElasticResXS);
		std::string fname_XS = "tests/table_total_elastic_XS.txt";
		std::string fname_XS1 = "tests/total_elastic_from_diff.txt";
		str.open(fname_XS, std::ios_base::trunc);
		str<<"E[eV]\tXS elastic and resonance total [1e-20 m^2]"<<std::endl;
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<ArTables.XS_elastic(E) + ArTables.XS_resonance_3o2(E)+ ArTables.XS_resonance_1o2(E)<<std::endl;
		}
		str.close();
		std::string name = "tests/test_table_XS.sc";
		str.open(name, std::ios_base::trunc);
		str<<"set logscale x"<<std::endl;
		str<<"plot \""<<fname_XS1<<"\" u 1:3 title \"XS from function\""<<std::endl;
		str<<"replot \""<<fname_XS<<"\" u 1:2 w lines title \"XS from table\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		std::string fname_back_prob = "tests/table_back_prob.txt";
		std::string fname_back_prob1 = "tests/backward_scattering_prob.txt";
		str.open(fname_back_prob, std::ios_base::trunc);
		str<<"E[eV]\tP elastic\tP resonance 3/2\tP resonance 1/2"<<std::endl;
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<ArTables.P_backward_elastic(E)<<"\t"<<ArTables.P_backward_resonance_3o2(E)<<"\t"<<ArTables.P_backward_resonance_1o2(E)<<std::endl;
		}
		str.close();
		std::string name = "tests/test_table_backward_prob.sc";
		str.open(name, std::ios_base::trunc);
		str<<"plot \""<<fname_back_prob1<<"\" u 1:2 title \"P from function\""<<std::endl;
		str<<"replot \""<<fname_back_prob<<"\" u 1:2 w lines title \"P from table elastic\""<<std::endl;
		str<<"replot \""<<fname_back_prob<<"\" u 1:3 title \"P from table resonance 3/2\""<<std::endl;
		str<<"replot \""<<fname_back_prob<<"\" u 1:4 title \"P from table resonance 1/2\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		std::string fname_TM_back = "tests/table_TM_backward.txt";
		std::string fname_TM_back1 = "tests/TM_backward.txt";
		str.open(fname_TM_back, std::ios_base::trunc);
		str<<"E[eV]\tTM backward elastic\tTM backward resonance 3/2\tTM backward resonance 1/2"<<std::endl;
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<ArTables.TM_backward_elastic(E)<<"\t"<<ArTables.TM_backward_resonance_3o2(E)<<"\t"<<ArTables.TM_backward_resonance_1o2(E)<<std::endl;
		}
		str.close();
		std::string name = "tests/test_table_TM_backward.sc";
		str.open(name, std::ios_base::trunc);
		str<<"plot \""<<fname_TM_back1<<"\" u 1:2 title \"TM backward from function\""<<std::endl;
		str<<"replot \""<<fname_TM_back<<"\" u 1:2 w lines title \"TM backward from table elastic\""<<std::endl;
		str<<"replot \""<<fname_TM_back<<"\" u 1:3 title \"TM backward from table resonance 3/2\""<<std::endl;
		str<<"replot \""<<fname_TM_back<<"\" u 1:4 title \"TM backward from table resonance 1/2\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		std::string fname_TM_for = "tests/table_TM_forward.txt";
		std::string fname_TM_for1 = "tests/TM_forward.txt";
		str.open(fname_TM_for, std::ios_base::trunc);
		str<<"E[eV]\tTM forward elastic\tTM forward resonance 3/2\tTM forward resonance 1/2"<<std::endl;
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<ArTables.TM_forward_elastic(E)<<"\t"<<ArTables.TM_forward_resonance_3o2(E)<<"\t"<<ArTables.TM_forward_resonance_1o2(E)<<std::endl;
		}
		str.close();
		std::string name = "tests/test_table_TM_forward.sc";
		str.open(name, std::ios_base::trunc);
		str<<"plot \""<<fname_TM_for1<<"\" u 1:2 title \"TM forward from function\""<<std::endl;
		str<<"replot \""<<fname_TM_for<<"\" u 1:2 w lines title \"TM forward from table elastic\""<<std::endl;
		str<<"replot \""<<fname_TM_for<<"\" u 1:3 title \"TM forward from table resonance 3/2\""<<std::endl;
		str<<"replot \""<<fname_TM_for<<"\" u 1:4 title \"TM forward from table resonance 1/2\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		EnergyScanner eScan(EnergyScanner::PlotResonances);
		std::string fname_XS = "tests/table_resonance_XS.txt";
		std::string fname_XS1 = "tests/resonance_XS.txt";
		str.open(fname_XS, std::ios_base::trunc);
		str<<"E[eV]\tXS resonances + elastic total [1e-20 m^2]"<<std::endl;
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<ArTables.XS_elastic(E)+ArTables.XS_resonance_3o2(E)+ArTables.XS_resonance_1o2(E)<<std::endl;
		}
		str.close();
		std::string name = "tests/test_table_resonance_XS.sc";
		str.open(name, std::ios_base::trunc);
		str<<"plot \""<<fname_XS1<<"\" u 1:2 title \"Resonance XS from function\""<<std::endl;
		str<<"replot \""<<fname_XS<<"\" u 1:2 title \"Resonance XS from table\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		std::string fname_XS = "data_derived/total_cross_section_integral.dat";
		std::string name = "tests/test_table_XS_integral.sc";
		str.open(name, std::ios_base::trunc);
		str<<"plot \""<<fname_XS<<"\" u 1:2 w lines title \"XS integral\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}
}

void test_resonance_cross (void)
{
	std::ofstream str;
	{
		EnergyScanner eScan(EnergyScanner::PlotResonances);
		std::string fname_XS = "tests/resonance_XS.txt";
		str.open(fname_XS, std::ios_base::trunc);
		str<<"E[eV]\tXS resonance total [1e-20 m^2]"<<std::endl;
		int err;
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str<<E<<"\t"<<argon_cross_elastic(E) + argon_cross_resonance_3o2(E)+argon_cross_resonance_1o2(E)<<std::endl;
		}
		str.close();
		std::string name = "tests/test_resonance_XS.sc";
		str.open(name, std::ios_base::trunc);
		str<<"plot \""<<fname_XS<<"\" u 1:2 title \"Resonance XS from function\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}
}

void test_all (void)
{
	/*
	std::cout<<"Testing polynimial fit:"<<std::endl;
	test_polynomial_fit ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing phase shifts fit:"<<std::endl;
	test_phase_shift_fit ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*/
	std::cout<<"Testing factor helping class:"<<std::endl;
	test_factor_helper ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing legendre polynomials:"<<std::endl;
	test_legendre_polynomial ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing colored intervals:"<<std::endl;
	test_colored_interval ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	/*
	std::cout<<"Testing integrals of legendre polynomials:"<<std::endl;
	test_legendre_intregral ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*//*
	std::cout<<"Testing differential cross section:"<<std::endl;
	test_diff_tot_cross ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing backward scatter probability:"<<std::endl;
	test_backward_scatter_prob ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing forward momentum transfer factor:"<<std::endl;
	test_TM_forward ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing backward momentum transfer factor:"<<std::endl;
	test_TM_backward ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*/
	std::cout<<"Testing resonance cross section:"<<std::endl;
	test_resonance_cross ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing total cross sections:"<<std::endl;
	test_total_cross_all ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing Ar data tables:"<<std::endl;
	test_data_table ();
	std::cout<<"==============================================="<<std::endl;

	std::cout<<"Testing finished."<<std::endl<<std::endl<<std::endl;
}
