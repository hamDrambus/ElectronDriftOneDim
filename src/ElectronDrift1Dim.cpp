//============================================================================
// Name        : ElectronDrift1Dim.cpp
// Author      : Frolov Egor geffdroid@gmail.com
// Version     :
// Copyright   : No copyrights
// Description :
//============================================================================

#include <iostream>
#include "argon_cross.h"
#include "Manager.h"
#include "tests.h"

int main(int argn, char * argv[]) {
	ensure_folder("Output");
	std::string a;
	test_all();
	Manager manman;
	double Td = 3; //=E/N in 1e-21 in Si
	double pressure = 1.015e5;
	double temperature = 87;
	double concentration = pressure / (temperature*boltzmann_SIconst);
	double field = Td*1e-21 * concentration;
	manman.SetParameters(concentration, field);
	manman.LoopSimulation();
	manman.WriteHistory("Output/eData_3Td.root");
	std::cout<<"Enter something: ";
	std::cin>>a;
	return 0;
}
