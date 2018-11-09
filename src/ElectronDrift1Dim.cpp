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
	std::string a;
	std::string root_fname = "Output/eData_3Td_.root";
	double Td = 3; //=E/N in 1e-21 in Si
	unsigned int seed = 42;
	unsigned int num_of_electrons = 1000;
	double pressure = 1.015e5;
	double temperature = 87;
	if (argn!=1) {
		if (argn<5) {
			std::cout<<"Input error, the following parameters must be specified:"<<std::endl;
			std::cout<<"(int)<seed>\t(double)<Td>\t(int)<number of electrons>\t(string)<output file name>"<<std::endl;
			return 0;
		}
		if (argn>5)
			std::cout<<"Warning! Only first 4 parameters are used."<<std::endl;
		seed = std::stoi(argv[1]);
		Td = std::stod(argv[2]);
		num_of_electrons = std::abs(std::stoi(argv[3]));
		root_fname = argv[4];
	}
	ensure_file(root_fname);
	ensure_file("tests/t");
	//test_all();
	//std::cout<<"Enter something: ";
	//std::cin>>a;
	//return 0;
	Manager manman(seed);
	//manman.Test();
	double concentration = pressure / (temperature*boltzmann_SIconst);
	double field = Td*1e-21 * concentration;
	auto start = std::chrono::system_clock::now();
	manman.SetParameters(concentration, field);
	unsigned int incr = num_of_electrons>1000 ? 1 + num_of_electrons/100: 1;
	for (unsigned int i=0;i<num_of_electrons;++i) {
		manman.LoopSimulation();
		if (0==(i+1)%incr)
			std::cout<<i+1<<"/"<<num_of_electrons<<std::endl;
	}
	if (0!=num_of_electrons%incr) {
		std::cout<<num_of_electrons<<"/"<<num_of_electrons<<std::endl;
	}
	manman.WriteHistory(root_fname);
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end-start;
	std::cout<<"Elapsed time:\t"<<diff.count()<<" s"<<std::endl;
	//std::cout<<"Enter something: ";
	//std::cin>>a;
	return 0;
}
