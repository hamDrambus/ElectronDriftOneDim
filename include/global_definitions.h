#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include <string>
#include <iostream>
#if defined(__WIN32__)
#define NOMINMAX
#include <windows.h>
#include <direct.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif

#define L_MAX_ 10
#define EN_MINIMUM 8e-3
#define THRESH_E_ 0.8
//^MERT 5 is applied only below this energy
#define EN_MAXIMUM_ 18
//^when elastic XS is still significantly larger than inelastic
#define En_3o2_ 11.103
#define En_1o2_ 11.270
#define Width_3o2_ 2.3e-3
#define Width_1o2_ 2.2e-3
//^in eV

const long double e_charge_SIconst = 1.60217662e-19; //in coulombs (SI)
const long double e_mass_SIconst = 9.10938356e-31; //in kg (SI)
const long double e_mass_eVconst = 5.109989461e5; //in eV
const long double h_bar_SIconst = 1.054571800e-34; //in SI
const long double a_bohr_SIconst = 5.2917721092e-11; //in meters (SI)
const long double a_h_bar_2e_m_e_SIconst = 2.711063352e-1; //in SI a_bohr*sqrt(2*Me*e)/h_bar
const long double boltzmann_SIconst = 1.38064852e-23; //SI
const double resonance_time_const = 2.9918725e-13;//in s. = h_bar /Width

#if defined(__WIN32__)
#define INVOKE_GNUPLOT(a) system(("start \"\" \"%GNUPLOT%\\gnuplot.exe\" --persist \"" + a + "\"").c_str())
#else
#define INVOKE_GNUPLOT(a) system(("konsole -e gnuplot \"" + a +"\"").c_str());
#endif //__WIN32__

std::string strtoken(std::string &in, std::string break_symbs);
void ensure_folder(std::string folder);

#endif
