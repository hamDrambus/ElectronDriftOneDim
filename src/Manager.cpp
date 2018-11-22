#include "Manager.h"

Manager::Manager(UInt_t RandomSeed) : is_ready_(false), eField_(0), Concentration_(0), Coefficient_(0), skip_counter_(0)
{
	random_generator_ = new TRandom1(RandomSeed); //TRandom3 is a default. TRandom1 is slower but better
	InitTree();
}

void Manager::InitTree (void)
{
	sim_data_ = new TTree("ElectronHistory", "ElectronHistory");
	sim_data_->Branch("energy_initial", &event_.En_start);
	sim_data_->Branch("energy_coll", &event_.En_collision);
	sim_data_->Branch("energy_final", &event_.En_finish);
	sim_data_->Branch("energy_average", &event_.En_avr);

	sim_data_->Branch("time_initial", &event_.time_start);
	sim_data_->Branch("time_delta", &event_.delta_time);
	sim_data_->Branch("time_delta_full", &event_.delta_time_full);
	sim_data_->Branch("process_type", &event_.process);

	sim_data_->Branch("position_initial",&event_.pos_start);
	sim_data_->Branch("position_final",&event_.pos_finish);
	sim_data_->Branch("position_delta",&event_.delta_x);

	/*sim_data_->Branch("deb_log_rand",&event_.deb_log_rand);
	sim_data_->Branch("deb_solver_y_left",&event_.deb_solver_y_left);
	sim_data_->Branch("deb_solver_y_right",&event_.deb_solver_y_right);
	sim_data_->Branch("deb_solver_E_left",&event_.deb_solver_E_left);
	sim_data_->Branch("deb_solver_E_right",&event_.deb_solver_E_right);*/
}

void Manager::Clear(void)
{
	is_ready_ = false;
	sim_data_->Delete();
	skip_counter_ = 0;
	InitTree();
}

void Manager::SetParameters(double Concentr /*in SI*/, double E /*in SI*/)
{
	Concentration_ = std::fabs(Concentr);
	eField_ = std::fabs(E);
	is_ready_ = true;
	if ((0 == Concentration_) || (0 == eField_)) {
		std::cout << "Manager::SetParameters(): Can't have 0 values" << std::endl;
		is_ready_ = false;
	}
	Coefficient_ = 1e20/* *e_charge_SIconst*/ * eField_ / Concentration_; // me/mu is set to 1.0
	//1e20 is because XS is expressed in 1e-20 m^2, not in m^2
}

void Manager::SetParameters(double T /*in K*/, double Pressure /*in SI*/, double E /*in SI*/)
{
	if (0 == T) {
		std::cout << "Manager::SetParameters(): Can't have 0 temperature" << std::endl;
		is_ready_ = false;
		return;
	}
	SetParameters(Pressure / (T*boltzmann_SIconst), E);
}

long double Manager::XS_integral(long double from, long double to)
{
	double E = from, E_abs;
	long double dx;
	long double Int = 0;
	while (E < to) {
		E_abs = std::fabs(E);
		//uses adaptive dx from the shape of XS(E).
		if (E_abs < 1e-2) {
			dx = 1e-3;
			goto l_out;
		}
		if (E_abs < 0.1) {
			dx = 1e-2;
			goto l_out;
		}
		if (E_abs < 1) {
			dx = 2e-2;
			goto l_out;
		}
		if (E_abs<10.7) {
			dx = 0.02 + (E_abs-1)*(0.1-0.02)/(9.7); //linear from 0.02 to 0.1 from 1 to 10.7 eV
			goto l_out;
		}
		dx = 0.001; //above 10.7
		l_out:
		E += dx;
		if (E > to)
			dx = dx - E + to;
		Int += dx*(ArTables.XS_elastic(E_abs) + ArTables.XS_resonance(E_abs));
	}
	return Int;
}

void Manager::Initialize(Event &event)
{
	if (!is_ready_)
		return;
	event.CrossSections.resize(ArExper.max_process_ID + 2, 0); //2==Elastic + Resonance
	event.CrossSectionsSum.resize(ArExper.max_process_ID + 2, 0);
	//event.En_start = 1 + 6*random_generator_->Uniform();
	event.En_start = 12;
	event.En_collision = 0;
	event.En_finish = 0;
	event.pos_start = 0;
	event.pos_finish = 0;
	event.delta_x = 0;
	event.time_start = 0;
	//event.velocity_start = true;
	//event.velocity_collision = true;
	//event.velocity_finish = true;
	event.delta_time = 0;
	event.delta_time_full = 0;
	event.process = Event::None;

	event_ = event;
}

//Solves equation -ln (R) = (N*me/(e*E*mu)) Int XS_total(E)dE from E_initial to E_collision
//where E final is unknown. Negative energy corresponds to electron moving against z axis. In this case is still dE>0
//because electron is decelerated E' = E (negative) + dE.
//So XS domain is symmetrically extended to negative energies for convenience.
//Since E[SI] = E[eV]*e, actual equation is -E ln(R)/N ~= Int XS_total
void Manager::Solve (long double LnR, Event &event)
{
	event.deb_log_rand = LnR;
	event.En_collision = EN_MAXIMUM_;
	long double I_max = XS_integral(event.En_start, EN_MAXIMUM_);
	long double e_finish = event.En_collision;
	double convergence_criteria = 5e-4*(std::min(std::fabs(event.En_start), std::fabs((double)e_finish)) );
	long double prev_solution = event.En_start;
	long double left = event.En_start, right = event.En_collision;
	long double f_left = -LnR, f_right = I_max - LnR, f_new;
	if (I_max < LnR) {
		event.process = Event::Overflow;
	} else {
		while (convergence_criteria < std::fabs(e_finish - prev_solution)) {
			prev_solution = e_finish;
			e_finish = (left*f_right - right*f_left) / (f_right - f_left);
			f_new = XS_integral(event.En_start, e_finish) - LnR;
			if (f_new < 0) {
				left = e_finish;
				f_left = f_new;
			} else {
				right = e_finish;
				f_right = f_new;
			}
			convergence_criteria = 5e-4*(std::min(std::fabs(event.En_start), std::fabs((double)e_finish)) );//std::max(2e-6, std::fabs(5e-4*event.En_finish));
		}
		event.En_collision = e_finish;
	}
	event.deb_solver_y_left = f_left;
	event.deb_solver_y_right = f_right;
	event.deb_solver_E_left = left;
	event.deb_solver_E_right = right;
}

void Manager::NewSolve (long double LnR, Event &event)
{
	event.deb_log_rand = LnR;
	double Int = ArTables.XS_integral(std::fabs(event.En_start));
	if (event.En_start < 0) {
		if (Int>LnR) {
			event.En_collision = -ArTables.XS_integral_find(Int-LnR, event);
		} else {
			event.En_collision = ArTables.XS_integral_find(LnR-Int, event); //e passes through 0 energy.
		}
	} else
		event.En_collision = ArTables.XS_integral_find(Int+LnR, event);
}

void Manager::DoStepLength(Event &event)
{
	if (!is_ready_)
		return;
	//event.En_start *= event.velocity_start ? 1 : -1;

	long double L = - log(random_generator_->Uniform());
	L *= Coefficient_; //Calculated once for fixed parameters;
	//solving L = XS_integral(Ei, Ec) for Ec===E collision.
	NewSolve(L, event);

	//Energy is in eV
	long double vel_0 = event.En_start > 0 ? sqrt(2.0*e_charge_SIconst*event.En_start / e_mass_SIconst) :
		-sqrt(-2.0 * e_charge_SIconst*event.En_start / e_mass_SIconst);
	long double vel_1 = event.En_collision > 0 ? sqrt(2.0*e_charge_SIconst*event.En_collision / e_mass_SIconst) :
		-sqrt(-2.0 * e_charge_SIconst*event.En_collision / e_mass_SIconst);
	event.delta_time = (vel_1 - vel_0)*e_mass_SIconst / (e_charge_SIconst * eField_); //in s
	event.delta_x = (std::fabs(event.En_collision) - std::fabs(event.En_start)) / eField_;
	event.pos_finish = event.pos_start + event.delta_x;

	/*event.En_start *= event.velocity_start ? 1 : -1;
	if (event.En_collision < 0) {
		event.velocity_collision = false;
		event.En_collision *= -1;
	} else
		event.velocity_collision = true;
	event.velocity_finish = event.velocity_collision;*/
	event.En_avr = std::fabs(event.En_start) + eField_*vel_0*event.delta_time / 4.0 + e_charge_SIconst*std::pow(eField_*event.delta_time, 2) / (6 * e_mass_SIconst);
	
	//!!!ENERGY CUT!!! TODO:remove
	event.En_collision = (std::fabs(event.En_collision)<EN_CUT_) ? ((event.En_collision>0) ? EN_CUT_ : EN_CUT_): event.En_collision;
	//!!!ENERGY CUT!!! TODO:remove
	event_ = event;
}

void Manager::DoScattering(Event &event)
{
	if (!is_ready_)
		return;
	for (int i=0, end_ = event.CrossSections.size(); i!=end_; ++i) {
		event.CrossSections[i] = ArTables.CrossSection(std::fabs(event.En_collision), i);//process type == index in cross section array.
		event.CrossSectionsSum[i] = std::max(event.CrossSections[i], 0.0) + ((i==0) ? 0.0 : event.CrossSectionsSum[i-1]);
	}
	for (int i=0, end_ = event.CrossSections.size(); i!=end_; ++i)
		event.CrossSectionsSum[i] /= event.CrossSectionsSum[end_-1];

	bool is_overflow = (event.process==Event::Overflow);
	double R1 = random_generator_->Uniform();
	for (int i=0, end_ = event.CrossSections.size(); i!=end_; ++i)
		if (R1<event.CrossSectionsSum[i]) {
			event.process = i;
			break;
		}

	double BackScatterProb;
	switch (event.process) {
		case (Event::Resonance): {
			BackScatterProb = ArTables.P_backward_resonance(std::fabs(event.En_collision));
			break;
		}
		case (Event::Elastic): {
			BackScatterProb = ArTables.P_backward_elastic(std::fabs(event.En_collision));
			break;
		}
		default: {
			BackScatterProb = 0.5;
		}
	}
	bool is_backward = (random_generator_->Uniform() < BackScatterProb);
	double EnergyLoss = 0;
	long double gamma_f = e_mass_eVconst/Ar_mass_eVconst;
	switch (event.process) {
		case (Event::Resonance): {
			double TM_factor = is_backward ? ArTables.TM_backward_resonance(std::fabs(event.En_collision))
					: ArTables.TM_forward_resonance(std::fabs(event.En_collision));
			EnergyLoss = 2 * TM_factor*event.En_collision*gamma_f /pow(1 + gamma_f, 2);
			break;
		}
		case (Event::Elastic): {
			double TM_factor = is_backward ? ArTables.TM_backward_elastic(std::fabs(event.En_collision))
					: ArTables.TM_forward_elastic(std::fabs(event.En_collision));
			EnergyLoss = 2 * TM_factor*event.En_collision*gamma_f /pow(1 + gamma_f, 2);
			break;
		}
		case (Event::Ionization): {
			InelasticProcess *p = ArExper.FindInelastic(event.process-Event::Ionization);
			if (NULL!=p) {
				EnergyLoss = (event.En_collision > 0) ? p->get_En_thresh() : -p->get_En_thresh();
				EnergyLoss += (event.En_collision - EnergyLoss)/2.0; //Consider that residual energy is equally divided between 2 electrons
			}
			break;
		}
		default: {
			InelasticProcess *p = ArExper.FindInelastic(event.process-Event::Ionization);
			if (NULL!=p)
				EnergyLoss = (event.En_collision > 0) ? p->get_En_thresh() : -p->get_En_thresh();
		}
	}
	event.En_finish = event.En_collision - EnergyLoss;
	if (is_backward)
		event.En_finish = -event.En_finish;

	event.delta_time_full = event.delta_time;
	if (event.process == Event::Resonance) {
		event.delta_time_full += random_generator_->Exp(resonance_time_const);
	}
	if (is_overflow)
		event.process = Event::Overflow;
	event_ = event;
}

void Manager::PostStepAction(Event &event)
{
	if (!is_ready_)
		return;
	if ((0==skip_counter_)||(IsFinished(event))) {
		sim_data_->Fill();
	}
	++skip_counter_;
	if (skip_counter_>SKIP_HISTORY_)
		skip_counter_=0;
}

void Manager::DoGotoNext(Event &event)
{
	if (Event::None!=event.process) { //== for the very first event
		event.time_start += event.delta_time_full;
		event.En_start = event.En_finish;
		//event.velocity_start = event.velocity_finish;
		event.pos_start = event.pos_finish;
		event.process = Event::None;
		event.delta_time = 0;
		event.delta_x = 0;
		event.delta_time_full = 0;
	}

	event_ = event;
}

void Manager::DoStep(Event &event)
{
	if (!is_ready_)
		return;
	DoGotoNext(event);
	DoStepLength(event);
	DoScattering(event);
	PostStepAction(event);
}

bool Manager::IsFinished(Event &event)
{
	if (!is_ready_)
		return true;
	return !(event.pos_finish < DRIFT_DISTANCE_);
	//return !((event.pos_finish < DRIFT_DISTANCE_)&&(event.pos_finish> -10*DRIFT_DISTANCE_));
}

void Manager::LoopSimulation(void)
{
	if (!is_ready_)
		return;
	Initialize(event_);
	while (!IsFinished(event_)) {
		DoStep(event_);
	}
}

void Manager::WriteHistory(std::string root_fname)
{
	TFile *file = new TFile(root_fname.c_str(), "RECREATE");
	file->cd();
	std::cout<<"Event number: "<<sim_data_->GetEntries()<<std::endl;
	sim_data_->Write("", TObject::kOverwrite);
	file->Close();
	std::ofstream str;
	root_fname.pop_back();
	root_fname.pop_back();
	root_fname.pop_back();
	root_fname.pop_back();
	root_fname.pop_back();
	str.open(root_fname+"process_legend.txt", std::ios_base::trunc);
	str<<"//tree->process value | meaning"<<std::endl;
	str<<Event::Overflow<<"\tOverflow (energy exceeds "<<EN_MAXIMUM_<<" eV)"<<std::endl;
	str<<Event::None<<"\tNone"<<std::endl;
	str<<Event::Elastic<<"\tElastic scattering"<<std::endl;
	str<<Event::Resonance<<"\tFeshbach resonances"<<std::endl;
	for (int proc = Event::Ionization, end_ = ArExper.max_process_ID + Event::Ionization; proc!=end_; ++proc) {
		InelasticProcess *p = ArExper.FindInelastic(proc-Event::Ionization);
		if (NULL!=p)
			str<<proc<<"\t"<<p->get_name()<<std::endl;
	}
	str.close();
}

void Manager::Test(void)
{
	int Order = ArTables.getOrder();
	int Nused = ArTables.getNused();

	std::string fname_XS_3_4 = "tests/Int_XS_3_4.txt";
	std::string fname_XS_high_3_4 = "tests/Int_XS_high_3_4.txt";
	std::string fname_XS_1_2 = "tests/Int_XS_1_2.txt";
	std::string fname_XS_high_1_2 = "tests/Int_XS_high_1_2.txt";

	ArTables.setOrder(3);
	ArTables.setNused(4);
	std::ofstream str;
	str.open(fname_XS_3_4, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	double E;
	for (unsigned int l=0; l<50001;++l) {
		E = -EN_MAXIMUM_ + l*2*EN_MAXIMUM_/50000.0;
		str<< E<<"\t"<<XS_integral(-EN_MAXIMUM_, E)<<std::endl;
	}
	str.close();

	str.open(fname_XS_high_3_4, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	double Int = 0;
	E = -EN_MAXIMUM_;
	double dE = 1e-5;
	while (E<EN_MAXIMUM_) {
		str<< E<<"\t"<<Int<<std::endl;
		Int+=ArTables.XS_elastic(std::fabs(E))*dE;
		E+=dE;
	}
	str.close();

	ArTables.setOrder(1);
	ArTables.setNused(2);
	str.open(fname_XS_1_2, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	for (unsigned int l=0; l<50001;++l) {
		E = -EN_MAXIMUM_ + l*2*EN_MAXIMUM_/50000.0;
		str<< E<<"\t"<<XS_integral(-EN_MAXIMUM_, E)<<std::endl;
	}
	str.close();

	str.open(fname_XS_high_1_2, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	Int = 0;
	E = -EN_MAXIMUM_;
	dE = 1e-5;
	while (E<EN_MAXIMUM_) {
		str<< E<<"\t"<<Int<<std::endl;
		Int+=ArTables.XS_elastic(std::fabs(E))*dE;
		E+=dE;
	}
	str.close();

	std::string name = std::string("tests/test_XS_integral.sc");
	str.open(name, std::ios_base::trunc);
	str<<"plot '"<<fname_XS_3_4<<"' u 1:2 title 'I(XS) 3,4'"<<std::endl;
	str<<"replot '"<<fname_XS_high_3_4 <<"' u 1:2 w line lc rgb \"#000000\" title 'I(XS) high acc. 3,4'"<<std::endl;
	str<<"replot '"<<fname_XS_1_2<<"' u 1:2 title 'I(XS) 1,2'"<<std::endl;
	str<<"replot '"<<fname_XS_high_1_2 <<"' u 1:2 w line lc rgb \"#FF0000\" title 'I(XS) high acc. 1,2'"<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);

	ArTables.setOrder(Order);
	ArTables.setNused(Nused);
}
