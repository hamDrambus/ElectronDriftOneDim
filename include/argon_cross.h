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
#include "global_definitions.h"
#include "PolynomialFit.h"
#include "LegendrePolynomials.h"

class EnergyScanner
{
protected:
	unsigned int i;
	int type_; //0 - from 5e-3 to 13 eV. other - 10.85 - 11.7 eV - near resonance
public:
	EnergyScanner(int type = 0);
	long double Next(int& err);
	void Reset(void);
};

class ArExperimental
{
public:
	std::vector<DataVector> phase_shifts_; //TODO: account for phase plus and phase minus
	DataVector total_elastic_cross;
	ArExperimental(void);
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
	double XS_elastic(double E);
	double XS_resonance(double E);
	double P_backward_elastic(double E);
	double P_backward_resonance(double E);
	double TM_backward_elastic(double E);
	double TM_backward_resonance(double E);
	double TM_forward_elastic(double E);
	double TM_forward_resonance(double E);
	double XS_integral(double E);//always from -EN_MAXIMUM_
	double XS_integral_find(double Int, Event &event);
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

