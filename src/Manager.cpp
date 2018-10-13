#include "Manager.h"

Manager::Manager(UInt_t RandomSeed = 42) : is_ready_(false)
{
	random_generator_ = new TRandom3(RandomSeed);
	sim_data_ = new TTree("ElectronHistory", "ElectronHistory");
	sim_data_->Branch("energy_initial", &event_.En_start);
	sim_data_->Branch("energy_final", &event_.En_finish);
	sim_data_->Branch("energy_average", &event_.En_avr);

	sim_data_->Branch("velocity_initial", &event_.velocity_start);
	sim_data_->Branch("velocity_final", &event_.velocity_finish);
	sim_data_->Branch("time_initial", &event_.time_start);
	sim_data_->Branch("time_delta", &event_.delta_time);
	sim_data_->Branch("process_type", &event_.process);
}

void Manager::Clear(void)
{
	is_ready_ = false;
	sim_data_->Delete();
	sim_data_ = new TTree("ElectronHistory", "ElectronHistory");
	sim_data_->Branch("energy_initial", &event_.En_start);
	sim_data_->Branch("energy_final", &event_.En_finish);
	sim_data_->Branch("energy_average", &event_.En_avr);

	sim_data_->Branch("velocity_initial", &event_.velocity_start);
	sim_data_->Branch("velocity_final", &event_.velocity_finish);
	sim_data_->Branch("time_initial", &event_.time_start);
	sim_data_->Branch("time_delta", &event_.delta_time);
	sim_data_->Branch("process_type", &event_.process);
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
	Coefficient_ = e_charge_SIconst * eField_ / Concentration_; // me/mu is set to 1.0
}

void Manager::SetParameters(double T /*in K*/, double Pressure /*in SI*/, double E /*in SI*/)
{
	if (0 == T) {
		std::cout << "Manager::SetParameters(): Can't have 0 teperature" << std::endl;
		is_ready_ = false;
		return;
	}
	SetParameters(Pressure / (T*boltzmann_SIconst), E);
}

long double Manager::XS_integral(double from, double to)
{
	double E = from, E_abs;
	double dx;
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
		dx = 0.02 + (E_abs-1)*(0.1-0.02)/(12.0); //linear from 0.02 to 0.1 from 1 to 13 eV
	l_out:
		E += dx;
		if (E > to)
			dx = dx - E - to;
		Int += dx*(ArTables.XS_elastic(E_abs) + ArTables.XS_resonance(E_abs));
	}
	return Int;
}

//Solves equation -ln (R) = (N*me/(e*E*mu)) Int XS_total(E)dE from E_initial to E_final
//where E final is unknown. Negative energy corresponds to electron moving agains z axis. In this case is still dE>0
//because electron is decelerated E' = E (negative) + dE.
//So XS domain is symmetrically extented to negative energies for convinience.
void Manager::DoStepLength(Event &event)
{
	if (!is_ready_)
		return;
	event.En_start *= event.velocity_start ? 1 : -1;

	long double L = - log(random_generator_->Uniform());
	L *= Coefficient_; //Calculated once for fixed parameters;
	//solving L = XS_integral(Ei, Ef) for Ef.
	long double I_max = XS_integral(event.En_start, EN_MAXIMUM_);
	if (I_max < L) {
		event.process = Event::Overflow;
		event.En_finish = EN_MAXIMUM_;
	} else {
		event.En_finish = EN_MAXIMUM_;
		double convergence_criteria = std::max(2e-6, std::fabs(5e-4*event.En_finish));
		double prev_solution = event.En_start;
		double left = event.En_start, right = event.En_finish;
		double f_left = -L, f_right = I_max - L, f_new;
		while (convergence_criteria < std::fabs(event.En_finish - prev_solution)) {
			prev_solution = event.En_finish;
			event.En_finish = (left*f_right - right*f_left) / (f_right - f_left);
			f_new = XS_integral(event.En_start, event.En_finish) - L;
			if (f_new < 0) {
				left = event.En_finish;
				f_left = f_new;
			} else {
				right = event.En_finish;
				f_right = f_new;
			}
			convergence_criteria = std::max(2e-6, std::fabs(5e-4*event.En_finish));
		}
	}

	//Energy is in eV
	double vel_0 = event.En_start > 0 ? sqrt(2.0*e_charge_SIconst*event.En_start / e_mass_SIconst) :
		-sqrt(-2.0 * e_charge_SIconst*event.En_start / e_mass_SIconst);
	double vel_1 = event.En_finish > 0 ? sqrt(2.0*e_charge_SIconst*event.En_finish / e_mass_SIconst) :
		-sqrt(-2.0 * e_charge_SIconst*event.En_finish / e_mass_SIconst);
	event.delta_time = (vel_1 - vel_0)*e_mass_SIconst / (e_charge_SIconst * eField_); //in s
	event.pos_finish = event.pos_start + (event.En_finish - event.En_start) / eField_;

	event.En_start *= event.velocity_start ? 1 : -1;
	if (event.En_finish < 0) {
		event.velocity_finish = false;
		event.En_finish *= -1;
	}
	event.En_avr = event.En_start + eField_*vel_0*event.delta_time / 4.0 + e_charge_SIconst*std::pow(eField_*event.delta_time, 2) / (6 * e_mass_SIconst);
	
	event_ = event;
}

void Manager::DoScattering(Event &event)
{
	if (!is_ready_)
		return;
	long double XS_elastic = ArTables.XS_elastic(event.En_finish);
	long double XS_resonance = ArTables.XS_resonance(event.En_finish);
	if (XS_resonance > 0) {
		if (random_generator_->Uniform() < (XS_resonance / (XS_resonance + XS_elastic)))
			event.process = Event::Resonance;
		else
			event.process = (event.process == Event::Overflow ? event.process : Event::Elastic);
	} else {
		event.process = Event::Resonance;
	}
	double BackScatterProb = (event.process == Event::Resonance ? ArTables.P_backward_resonance(event.En_finish) :
		ArTables.P_backward_elastic(event.En_finish));
	long double gamma = (event.En_finish + e_mass_eVconst) / e_mass_eVconst;
	long double TM_factor;
	if (random_generator_->Uniform() < BackScatterProb) {
		TM_factor = (event.process == Event::Resonance ? ArTables.TM_backward_resonance(event.En_finish) :
			ArTables.TM_backward_elastic(event.En_finish));
		event.velocity_finish = !event.velocity_finish;
	} else {
		TM_factor = (event.process == Event::Resonance ? ArTables.TM_forward_resonance(event.En_finish) :
			ArTables.TM_forward_elastic(event.En_finish));
	}
	event.En_finish -= 2 * TM_factor*event.En_finish*gamma /pow(1 + gamma, 2);
	event_ = event;
}

void Manager::DoGotoNext(Event &event)
{
	double delta_time = 0;
	if (event.process == Event::Resonance) {
		delta_time = random_generator_->Exp(resonance_time_const);
	}
	event.time_start += event.time_start + delta_time + event.delta_time;
	event.En_start = event.En_finish;
	event.velocity_start = event.velocity_finish;
	event.pos_start = event.pos_finish;
	event.process = Event::None;

	event_ = event;
}

void Manager::Initialize(Event &event)
{
	if (!is_ready_)
		return;
	event.En_start = 0;
	event.pos_start = 0;
	event.time_start = 0;
	event.velocity_start = true;
	event.process = Event::None;
}

void Manager::DoStep(Event &event)
{
	if (is_ready_)
		return;
	DoStepLength(event);
	DoScattering(event);
	PostStepAction(event);
	DoGotoNext(event);
}

void Manager::PostStepAction(Event &event)
{
	if (!is_ready_)
		return;
	sim_data_->Fill();
}

bool Manager::IsFinished(Event &event)
{
	if (!is_ready_)
		return true;
	return !(event.pos_finish < 1e-7); //100 pm.
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
	sim_data_->Write("", TObject::kOverwrite);
	file->Close();
}
