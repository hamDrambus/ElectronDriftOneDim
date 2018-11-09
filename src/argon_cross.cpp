#include "argon_cross.h"

ArExperimental ArExper;
ArDataTables ArTables;


EnergyScanner::EnergyScanner(int type): i(0), type_(type)
{}
long double EnergyScanner::Next(int& err)
{
	if (0==type_) {
		if (i<700) { //from 5e-3 to 13.0 eV
			long double out = (i<425 ? (5e-3 + i*(0.5-5e-3)/425) : ( i<575 ? (0.5 + (i-425)*(3.0-0.5)/150) : (3.0 + (i-575)*(13.0-3.0)/124 ) ) );
			++i;
			err = 0;
			return out;
		}
		err = 1;
		Reset();
		return DBL_MAX;
	} else {
		if (i<200) { //from 10.85 to 11.7 eV
			long double out = 10.85 + i*(11.7 - 10.85)/199;
			++i;
			err = 0;
			return out;
		}
		err = 1;
		Reset();
		return DBL_MAX;
	}
}
void EnergyScanner::Reset(void)
{
	i = 0;
}

ArExperimental::ArExperimental(void): total_elastic_cross(3, 5) /*fit by 3rd order polynomial*/
{
	std::ifstream inp;
	inp.open("data/ArScatteringCross.dat");
	std::string line, word;
	while (!inp.eof()) {
		std::getline(inp, line);
		if (line.size()>=2)
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double k = a_h_bar_2e_m_e_SIconst*sqrt(std::stod(word));
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double XS = std::stod(word); //Cross section is stored in 1e-20 m^2 in both files and this program.
		total_elastic_cross.push(k, XS);
	}
	inp.close();
	for (int l=0;l<6;++l)
		if (l<4)
			phase_shifts_.push_back(DataVector(3,4)); /*interpolation by 3rd order polynomial because there's data for 11eV*/
		else
			phase_shifts_.push_back(DataVector(3,5)); /*fit by 3rd order polynomial TODO: tests, tests*/
	inp.open("data/McEachranArPhaseShifts.dat");
	while (!inp.eof()) {
		std::getline(inp, line);
		if (line.size()>=2) {
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
			if ((line[0]=='\t')&&(line[1]=='\t')) //negative phase shifts are ignored
				continue;
		}
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double k = std::stod(word);
		//now process string "double\tdouble\tdouble\t..."
		//\t\t means that the double is omitted
		std::vector<double> vals;
		std::vector<bool> is_d;//is double presented in the table
		vals.resize(9,0);
		is_d.resize(9, false);
		for (int l=0;l<9;++l) {
			if (line.empty())
				break;
			if (line[0]=='\t') { //value is omitted in the table
				line.erase(line.begin());
				continue;
			}
			//now line does not start with \t
			word = strtoken(line,"\t"); //removes one \t after value
			vals[l] = std::stod(word);
			is_d[l] = true;
		}
		for (int l=0;l<6;++l) {
			if (is_d[l])
				phase_shifts_[l].push(k, vals[l]);
		}
	}
	inp.close();
}

unsigned int ArExperimental::max_L (long double k)
{
	return phase_shifts_.size()-1;
}

long double ArExperimental::phase_shift (long double k, unsigned int l)
{
	if (l>max_L(k))
		return 0;
	return phase_shifts_[l](k, k);
}

void ArDataTables::read_data (std::ifstream &inp, DataVector &data, long double y_factor)
{
	std::string line, word;
	while (!inp.eof()&&inp.is_open()) {
		std::getline(inp, line);
		if (line.size()>=2)
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double E = std::stod(word);
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double XS = std::stod(word)*y_factor;
		data.push_back(E, XS);
	}
}

ArDataTables::ArDataTables():
	total_cross_elastic_fname("data_derived/total_cross_section_elastic.dat"),
	total_cross_resonance_fname("data_derived/total_cross_section_resonance.dat"),
	back_scatter_elastic_prob_fname("data_derived/backward_scatter_elastic_prob.dat"),
	back_scatter_resonance_prob_fname("data_derived/backward_scatter_resonance_prob.dat"),
	TM_backward_elastic_fname("data_derived/TM_backward_elastic.dat"),
	TM_backward_resonance_fname("data_derived/TM_backward_resonance.dat"),
	TM_forward_elastic_fname("data_derived/TM_forward_elastic.dat"),
	TM_forward_resonance_fname("data_derived/TM_forward_resonance.dat"),
	total_cross_integral_fname("data_derived/total_cross_section_integral.dat"),
	total_cross_elastic_(3,4), //interpolation with 3rd order polynomial
	total_cross_resonance_(3,4),
	back_scatter_elastic_prob_(3,4),
	back_scatter_resonance_prob_(3,4),
	TM_backward_elastic_(3,4),
	TM_backward_resonance_(3,4),
	TM_forward_elastic_(3,4),
	TM_forward_resonance_(3,4),
	total_cross_integral_(1,2)
{
	std::cout<<"Constructing Ar data tables"<<std::endl;
	ensure_file(total_cross_elastic_fname);
	ensure_file(total_cross_resonance_fname);
	ensure_file(back_scatter_elastic_prob_fname);
	ensure_file(back_scatter_resonance_prob_fname);
	ensure_file(TM_backward_elastic_fname);
	ensure_file(TM_backward_resonance_fname);
	ensure_file(TM_forward_elastic_fname);
	ensure_file(TM_forward_resonance_fname);
	ensure_file(total_cross_integral_fname);
	total_cross_resonance_.enable_out_value(0);
	back_scatter_resonance_prob_.enable_out_value(0);
	TM_backward_resonance_.enable_out_value(0);
	TM_forward_resonance_.enable_out_value(0);

	std::ifstream inp;
	std::ofstream str;
	EnergyScanner allRange(0), resRange(1);
	int err;
	inp.open(total_cross_elastic_fname);
	read_data(inp, total_cross_elastic_); //Cross section is stored in 1e-20 m^2 in both files and this program
	inp.close();
	if (total_cross_elastic_.size()<total_cross_elastic_.getNused()) {
		std::cout<<"Calculating total elastic cross section..."<<std::endl;
		while (0!=total_cross_elastic_.size())
			total_cross_elastic_.erase(0);
		str.open(total_cross_elastic_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tXS elastic [1e-20 m^2]"<<std::endl;
		double E, cross;
		while (true) {
			E = allRange.Next(err);
			if (0!=err)
				break;
			cross = argon_cross_elastic(E);
			str<<E<<"\t"<<cross<<std::endl;
			total_cross_elastic_.push_back(E, cross);
		}
		str.close();
	}

	inp.open(total_cross_resonance_fname);
	read_data(inp, total_cross_resonance_);
	inp.close();
	if (total_cross_resonance_.size()<total_cross_resonance_.getNused()) {
		std::cout<<"Calculating total resonance cross section..."<<std::endl;
		while (0!=total_cross_resonance_.size())
			total_cross_resonance_.erase(0);
		str.open(total_cross_resonance_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tXS resonance [1e-20 m^2]"<<std::endl;
		double E, cross;
		while (true) {
			E = resRange.Next(err);
			if (0!=err)
				break;
			cross = argon_cross_resonance(E);
			str<<E<<"\t"<<cross<<std::endl;
			total_cross_resonance_.push_back(E, cross);
		}
		str.close();
	}

	inp.open(back_scatter_elastic_prob_fname);
	read_data(inp, back_scatter_elastic_prob_);
	inp.close();
	if (back_scatter_elastic_prob_.size()<back_scatter_elastic_prob_.getNused()) {
		std::cout<<"Calculating backward elastic scatter probability..."<<std::endl;
		while (0!=back_scatter_elastic_prob_.size())
			back_scatter_elastic_prob_.erase(0);
		str.open(back_scatter_elastic_prob_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tBack scatter elastic probability"<<std::endl;
		double E, prob;
		while (true) {
			E = allRange.Next(err);
			if (0!=err)
				break;
			prob = argon_back_scatter_prob(E);
			str<<E<<"\t"<<prob<<std::endl;
			back_scatter_elastic_prob_.push_back(E, prob);

		}
		str.close();
	}

	inp.open(back_scatter_resonance_prob_fname);
	read_data(inp, back_scatter_resonance_prob_);
	inp.close();
	if (back_scatter_resonance_prob_.size()<back_scatter_resonance_prob_.getNused()) {
		std::cout<<"Calculating backward resonance scatter probability..."<<std::endl;
		while (0!=back_scatter_resonance_prob_.size())
			back_scatter_resonance_prob_.erase(0);
		str.open(back_scatter_resonance_prob_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tBack scatter resonance probability"<<std::endl;
		double E, prob;
		while (true) {
			E = resRange.Next(err);
			if (0!=err)
				break;
			prob = argon_back_resonance_prob(E);
			str<<E<<"\t"<<prob<<std::endl;
			back_scatter_resonance_prob_.push_back(E, prob);
		}
		str.close();
	}

	inp.open(TM_backward_elastic_fname);
	read_data(inp, TM_backward_elastic_);
	inp.close();
	if (TM_backward_elastic_.size()<TM_backward_elastic_.getNused()) {
		std::cout<<"Calculating backward elastic momentum transfer factor..."<<std::endl;
		while (0!=TM_backward_elastic_.size())
			TM_backward_elastic_.erase(0);
		str.open(TM_backward_elastic_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tTM backward elastic factor"<<std::endl;
		double E, prob;
		while (true) {
			E = allRange.Next(err);
			if (0!=err)
				break;
			prob = argon_TM_backward(E);
			str<<E<<"\t"<<prob<<std::endl;
			TM_backward_elastic_.push_back(E, prob);
		}
		str.close();
	}

	inp.open(TM_backward_resonance_fname);
	read_data(inp, TM_backward_resonance_);
	inp.close();
	if (TM_backward_resonance_.size()<TM_backward_resonance_.getNused()) {
		std::cout<<"Calculating backward resonance momentum transfer factor..."<<std::endl;
		while (0!=TM_backward_resonance_.size())
			TM_backward_resonance_.erase(0);
		str.open(TM_backward_resonance_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tTM backward resonance factor"<<std::endl;
		double E, prob;
		while (true) {
			E = resRange.Next(err);
			if (0!=err)
				break;
			prob = argon_TM_backward_resonance(E);
			str<<E<<"\t"<<prob<<std::endl;
			TM_backward_resonance_.push_back(E, prob);
		}
		str.close();
	}

	inp.open(TM_forward_elastic_fname);
	read_data(inp, TM_forward_elastic_);
	inp.close();
	if (TM_forward_elastic_.size()<TM_forward_elastic_.getNused()) {
		std::cout<<"Calculating forward elastic momentum transfer factor..."<<std::endl;
		while (0!=TM_forward_elastic_.size())
			TM_forward_elastic_.erase(0);
		str.open(TM_forward_elastic_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tTM forward elastic factor"<<std::endl;
		double E, prob;
		while (true) {
			E = allRange.Next(err);
			if (0!=err)
				break;
			prob = argon_TM_forward(E);
			str<<E<<"\t"<<prob<<std::endl;
			TM_forward_elastic_.push_back(E, prob);
		}
		str.close();
	}

	inp.open(TM_forward_resonance_fname);
	read_data(inp, TM_forward_resonance_);
	inp.close();
	if (TM_forward_resonance_.size()<TM_forward_resonance_.getNused()) {
		std::cout<<"Calculating forward resonance momentum transfer factor..."<<std::endl;
		while (0!=TM_forward_resonance_.size())
			TM_forward_resonance_.erase(0);
		str.open(TM_forward_resonance_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tTM forward resonance factor"<<std::endl;
		double E, prob;
		while (true) {
			E = resRange.Next(err);
			if (0!=err)
				break;
			prob = argon_TM_forward_resonance(E);
			str<<E<<"\t"<<prob<<std::endl;
			TM_forward_resonance_.push_back(E, prob);
		}
		str.close();
	}

	inp.open(total_cross_integral_fname);
	read_data(inp, total_cross_integral_);
	inp.close();
	if (total_cross_integral_.size()<total_cross_integral_.getNused()) {
		std::cout<<"Calculating integral of total cross section..."<<std::endl;
		while (0!=total_cross_integral_.size())
			total_cross_integral_.erase(0);
		str.open(total_cross_integral_fname, std::ios_base::trunc);
		//Fill
		str<<"//E[eV]\tInt{XS_elastic}[1e-20 m^2 eV]"<<std::endl;
		double Int = 0;
		double E = -EN_MAXIMUM_;
		double dE = 1e-5;
		while (E<=EN_MAXIMUM_) {
			str<<E<<"\t"<<Int<<std::endl;
			total_cross_integral_.push_back(E, Int);
			Int+=(XS_elastic(std::fabs(E))+XS_resonance(std::fabs(E)))*dE;
			E+=dE;
		}
		str.close();
	}
	std::cout<<"Finished constructing Ar data tables"<<std::endl;
}

double ArDataTables::XS_elastic(double E)
{
	if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value. TODO: maybe make linear increase of XS to 0
		E = EN_MINIMUM;
	return total_cross_elastic_(E, E);
}

double ArDataTables::XS_resonance(double E)
{
	if (E<EN_MINIMUM)
		E = EN_MINIMUM;
	return total_cross_resonance_(E, E);
}
double ArDataTables::P_backward_elastic(double E)
{
	if (E<EN_MINIMUM)
		E = EN_MINIMUM;
	//return 0.31;
	return back_scatter_elastic_prob_(E, E);
}
double ArDataTables::P_backward_resonance(double E)
{	
	if (E<EN_MINIMUM)
		E = EN_MINIMUM;
	//return 0.31;
	return back_scatter_resonance_prob_(E, E);
}
double ArDataTables::TM_backward_elastic(double E)
{
	if (E<EN_MINIMUM)
		E = EN_MINIMUM;
	return TM_backward_elastic_(E, E); 
}
double ArDataTables::TM_backward_resonance(double E)
{
	if (E<EN_MINIMUM)
		E = EN_MINIMUM;
	return TM_backward_resonance_(E, E);
}
double ArDataTables::TM_forward_elastic(double E)
{
	if (E<EN_MINIMUM)
		E = EN_MINIMUM;
	return TM_forward_elastic_(E, E); 
}
double ArDataTables::TM_forward_resonance(double E)
{
	if (E<EN_MINIMUM)
		E = EN_MINIMUM;
	return TM_forward_resonance_(E, E); 
}

//always from -EN_MAXIMUM_
double ArDataTables::XS_integral(double E)
{
	return total_cross_integral_(E, E);
}
//finds E value corresponding to Int value of integral. Int==XS_integral(returned value)
double ArDataTables::XS_integral_find(double Int, Event &event)
{
	event.deb_solver_y_left=0;
	event.deb_solver_y_right=0;
	if (Int<0) {
		event.deb_solver_E_left=-EN_MAXIMUM_;
		event.deb_solver_E_right=-EN_MAXIMUM_;
		event.process = Event::Overflow;
		return -EN_MAXIMUM_;
	}
	std::size_t sz = total_cross_integral_.size();
	for (std::size_t i=1; i<sz; ++i) {
		if (total_cross_integral_.getY(i)>Int) { //Int is between y[i-1] and y[i]
			//linear formula
			event.deb_solver_y_left=total_cross_integral_.getY(i-1);
			event.deb_solver_y_right=total_cross_integral_.getY(i);
			event.deb_solver_E_left=total_cross_integral_.getX(i-1);
			event.deb_solver_E_right=total_cross_integral_.getX(i);
			return event.deb_solver_E_left + (Int-event.deb_solver_y_left)*(event.deb_solver_E_right-event.deb_solver_E_left)/(event.deb_solver_y_right-event.deb_solver_y_left);
		}
	}
	event.deb_solver_E_left=EN_MAXIMUM_;
	event.deb_solver_E_right=EN_MAXIMUM_;
	event.process = Event::Overflow;
	return EN_MAXIMUM_;
}

void ArDataTables::setOrder(int order)
{
	total_cross_elastic_.setOrder(order);
	total_cross_resonance_.setOrder(order);
	back_scatter_elastic_prob_.setOrder(order);
	back_scatter_resonance_prob_.setOrder(order);
	TM_backward_elastic_.setOrder(order);
	TM_backward_resonance_.setOrder(order);
	TM_forward_elastic_.setOrder(order);
	TM_forward_resonance_.setOrder(order);
}

void ArDataTables::setNused(int N)
{
	total_cross_elastic_.setNused(N);
	total_cross_resonance_.setNused(N);
	back_scatter_elastic_prob_.setNused(N);
	back_scatter_resonance_prob_.setNused(N);
	TM_backward_elastic_.setNused(N);
	TM_backward_resonance_.setNused(N);
	TM_forward_elastic_.setNused(N);
	TM_forward_resonance_.setNused(N);
}

int ArDataTables::getOrder(void)
{
	return total_cross_elastic_.getOrder();
}

int ArDataTables::getNused(void)
{
	return total_cross_elastic_.getNused();
}

//k is in atomic units
void argon_phase_values_exp(long double k, unsigned int l, long double &tan, long double &sin, long double &cos)
{
	double angle = ArExper.phase_shift(k, l);
	tan = std::tan(angle);
	sin = std::sin(angle);
	cos = std::cos(angle);
}

//k is in atomic units
void argon_phase_values_MERT5(long double k, unsigned int l, long double &tan, long double &sin, long double &cos)
{
	//see Kurokawa Phys. Rev. A84 2011, MERT5+ fit http://dx.doi.org/10.1103/PhysRevA.84.062717
	double A = -1.365;
	double D = 80.5;
	double F = -153;
	double G = 31.0;
	double A1 = 8.8;
	double H = 29.7;
	double alpha_d = 11.08;
	double alpha_q = 0.0;
	if (0==l) {
		tan = -A*k*(1+4*alpha_d*k*k*log(k)/3)-M_PI*alpha_d*k*k/3 + D*pow(k,3) + F* pow(k,4);
		tan/=(1+G*pow(k,3));
	} else {
		unsigned int l2 = 2*l;
		long double al = M_PI/((l2+3)*(l2+1)*(l2-1));
		long double bl = M_PI*(15*pow(l2+1,4)-140*pow(l2+1,2)+128)/(pow((l2+3)*(l2+1)*(l2-1),3)*(l2+5)*(l2-3));
		long double cl = 3*al/((l2+5)*(l2-3));
		tan = al*alpha_d*k*k+(bl*pow(alpha_d,2)+cl*alpha_q)*pow(k,4);
		if (1==l)
			tan+=H*pow(k,5) - A1*pow(k,3);
	}
	sin = fabs(tan)*sqrt(1/(1+tan*tan));
	cos = (tan>0?1.0:-1.0)*sqrt(1/(1+tan*tan));
}

//E in eV
long double argon_cross_elastic_diff (long double E, long double theta) {
	if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value. TODO: maybe make linear increase of XS to 0
		E= EN_MINIMUM;
	//different formulas are used for E<1eV and E>1eV!
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double cross = 0;
	long double cos_th = cos(theta);
	unsigned int L_MAX = (E<THRESH_E_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			cross+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*P1(cos_th, l)*P2(cos_th, f);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			//P1 and P2 because each of them has separate cache.
		}
	}
	cross*=2*M_PI/pow(k,2);
	return cross*a_bohr_SIconst*a_bohr_SIconst; //const is multiplied bt 1e10
}

long double argon_cross_elastic (long double E)
{
	if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value. TODO: maybe make linear increase of XS to 0
		E = EN_MINIMUM;
	if (E <THRESH_E_XS_) {
		long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
		// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
		long double cross = 0;
		unsigned int L_MAX = L_MAX_;
		for (unsigned int l=0; l<=L_MAX; ++l) {
			long double sin_phase_l = 0;
			long double cos_phase_l = 1;
			long double tan_l = 0;
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
			cross+=(2*l+1)*sin_phase_l*sin_phase_l;
		}
		cross*=4*M_PI/pow(k,2);
		return cross*a_bohr_SIconst*a_bohr_SIconst;
	} else {
		return ArExper.total_elastic_cross(a_h_bar_2e_m_e_SIconst*sqrt(E), a_h_bar_2e_m_e_SIconst*sqrt(E));
	}
	return 0;
}

//for testing only
long double argon_cross_elastic_from_phases (long double E)
{
	if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value. TODO: maybe make linear increase of XS to 0
		E = EN_MINIMUM;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	long double cross = 0;
	unsigned int L_MAX = (E<THRESH_E_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		cross+=(2*l+1)*sin_phase_l*sin_phase_l;
	}
	cross*=4*M_PI/pow(k,2);
	return cross*a_bohr_SIconst*a_bohr_SIconst;
}

long double argon_back_scatter_prob (long double E)
{
	if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value. TODO: maybe make linear increase of XS to 0
		E = EN_MINIMUM;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			W+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl(l,f, -1, 0, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
		}
		cross+=(2*l+1)*sin_phase_l*sin_phase_l;
	}
	return W/(2*cross);
}
//energy loss and input are in eV
long double argon_TM_forward (long double E)
{
	if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value. TODO: maybe make linear increase of XS to 0
		E = EN_MINIMUM;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			W+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl_transf(l,f, 0, 1, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl(l,f, 0, 1, 1e-5);
		}
	}
	return W/cross;
}

long double argon_TM_backward (long double E)
{
	if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value. TODO: maybe make linear increase of XS to 0
		E = EN_MINIMUM;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			W+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl_transf(l,f, -1, 0, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl(l,f, -1, 0, 1e-5);
		}
	}
	return W/cross;
}

long double argon_cross_resonance_diff (long double E, long double theta) //TODO: maybe need to implement
{
	return 0;
}

long double argon_cross_resonance (long double E)
{
	if ((E<10.85)||(E>11.7))
		return 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E);
	long double sin_phase_1 = 0;
	long double cos_phase_1 = 1;
	long double tan_1 = 0;
	argon_phase_values_exp(k, 1, tan_1, sin_phase_1, cos_phase_1);
	long double cot = 2*(E - En_3o2_ )/Width_3o2_;
	//resonance phase changes from 0 to -Pi when energy changes from -inf to +inf, hence the choice of signs
	long double sin_3o2 = - sqrt(1.0/(1+pow( cot , 2)));
	long double cos_3o2 = - cot * sin_3o2;
	long double cross =8*M_PI*(pow(sin_phase_1*cos_3o2 + cos_phase_1*sin_3o2, 2) - sin_phase_1*sin_phase_1);
	cot = 2*(E - En_1o2_ )/Width_1o2_;
	long double sin_1o2 = - sqrt(1.0/(1+pow( cot , 2)));
	long double cos_1o2 = - cot * sin_1o2;
	cross +=4*M_PI*(pow(sin_phase_1*cos_1o2 + cos_phase_1*sin_1o2, 2) - sin_phase_1*sin_phase_1);
	cross /=pow (k, 2);
	return cross * a_bohr_SIconst * a_bohr_SIconst;
}

long double argon_back_resonance_prob (long double E)
{
	return argon_back_scatter_prob(E);
}
long double argon_TM_forward_resonance (long double E)
{
	return argon_TM_forward(E);
}
long double argon_TM_backward_resonance (long double E)
{
	return argon_TM_backward(E);;
}
