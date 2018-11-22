#ifndef ARGON_CROSS_H
#define ARGON_CROSS_H

/*	TODO: add comprehensive explanation of used data and calculations of cross sections
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <tgmath.h>
#include "global_definitions.h"
#include "PolynomialFit.h"
#include "LegendrePolynomials.h"
#include "ColoredInterval.h"

class EnergyScanner
{
protected:
	unsigned int i;
	enum ScanType: short {ElasticXS = 0, ResonanceXS = 1, ResonanceDiffXS = 2,
		DiffXS = 3,	InelasticXS = 4, ElasticResXS = 5, XSIntegral = 6, PlotElastic = 7,
		PlotResonance = 8, PlotDiff = 9, PlotInelastic = 10, PlotAllXS = 11} type_;
	ColoredRange energy_range_;
public:
	EnergyScanner(ScanType type);
	long double Next(int& err);
	void Reset(void);
};

class InelasticProcess
{
protected:
	std::string name_;
	DataVector exp_XS_;
	unsigned int ID_;
	double En_threshold_;
	double Oscillator_strength_;
public:
	InelasticProcess(std::string name, unsigned int ID, double En, double F, std::vector<double> &Ens, std::vector<double> &XSs);
	double operator ()(double E); //returns cross section in 1e-16 cm^2
	double BB_XS(double E);
	double Exp_XS(double E);
	double get_En_thresh (void) const;
	std::string get_name (void) const;
	unsigned int get_ID (void) const;
};

class ArExperimental
{
protected:
	void read_inelastic(std::ifstream &inp, std::vector<InelasticProcess> &to);
public:
	std::vector<DataVector> phase_shifts_; //TODO: account for phase plus and phase minus
	DataVector total_elastic_cross;
	std::vector<InelasticProcess> excitations;
	std::vector<InelasticProcess> ionizations;
	short max_process_ID;
	double E_Ionization;
	ArExperimental(void);
	InelasticProcess * FindInelastic(short ID);
	unsigned int max_L (long double k);
	long double phase_shift (long double k, unsigned int l);
};

class ArDataTables //loads data from default files if presented. If not then values are calculated and files are created.
{
protected:
	std::string total_cross_elastic_fname;
	std::string total_cross_resonance_fname;
	std::string back_scatter_elastic_prob_fname;
	std::string back_scatter_resonance_prob_fname;
	std::string TM_backward_elastic_fname;
	std::string TM_backward_resonance_fname;
	std::string TM_forward_elastic_fname;
	std::string TM_forward_resonance_fname;
	std::string total_cross_integral_fname;
	//atm backward scatter probability and TM's for resonance are taken the same as elastic values
	DataVector total_cross_elastic_;
	DataVector total_cross_resonance_;
	DataVector back_scatter_elastic_prob_;
	DataVector back_scatter_resonance_prob_;
	DataVector TM_backward_elastic_;
	DataVector TM_backward_resonance_;
	DataVector TM_forward_elastic_;
	DataVector TM_forward_resonance_;
	DataVector total_cross_integral_;
	void read_data (std::ifstream &inp, DataVector &data, long double y_factor = 1);
public:
	ArDataTables();
	long double CrossSection (double E, short type);
	long double TotalCrossSection (double E);
	long double XS_elastic(double E);
	long double XS_resonance(double E);
	double P_backward_elastic(double E);
	double P_backward_resonance(double E);
	double TM_backward_elastic(double E);
	double TM_backward_resonance(double E);
	double TM_forward_elastic(double E);
	double TM_forward_resonance(double E);
	long double XS_integral(double E);//always from -EN_MAXIMUM_
	long double XS_integral_find(long double Int, Event &event);
	//^finds E value corresponding to Int value of integral. Int==XS_integral(returned value)
	void setOrder(int order);
	void setNused(int N);
	int getOrder(void);
	int getNused(void);
};

extern ArExperimental ArExper;
extern ArDataTables ArTables;

void argon_phase_values_exp(long double k, unsigned int l, long double &tan, long double &sin, long double &cos);
void argon_phase_values_MERT5(long double k, unsigned int l, long double &tan, long double &sin, long double &cos);
//E in eV, theta in radians, output is in m
long double argon_cross_elastic_diff (long double E, long double theta);
long double argon_cross_elastic (long double E);
long double argon_cross_elastic_from_phases (long double E);
long double argon_back_scatter_prob (long double E);
long double argon_TM_forward (long double E);
long double argon_TM_backward (long double E);

long double argon_cross_resonance_diff (long double E, long double theta);
long double argon_cross_resonance (long double E);
long double argon_back_resonance_prob (long double E);
long double argon_TM_forward_resonance (long double E);
long double argon_TM_backward_resonance (long double E);

#endif

