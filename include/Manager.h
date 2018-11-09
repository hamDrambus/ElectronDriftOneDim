#ifndef MANAGER_H_
#define MANAGER_H_

#include <fstream>

#include <TRandom3.h>
#include <TRandom1.h>
#include <TFile.h>
#include <TTree.h>

#include "global_definitions.h"
#include "argon_cross.h"

class Manager
{
protected:
	TRandom * random_generator_;
	TTree * sim_data_;
	Event event_;
	std::vector<double> last_pos;
	//Event current_event;
	void DoStepLength (Event &event);
	void Solve (long double LnR, Event &event);
	void NewSolve (long double LnR, Event &event);
	void DoScattering(Event &event);
	void PostStepAction(Event &event);
	void DoGotoNext(Event &event);
	void InitTree (void);
	bool is_ready_;
	double Concentration_;
	double eField_;
	double Coefficient_;
	long double XS_integral(long double from, long double to);
public:
	void Initialize(Event &event);
	void Clear(void);
	void DoStep(Event &event);
	bool IsFinished(Event &event);
	void SetParameters(double Concetr /*in SI*/, double E /*in SI*/);
	void SetParameters(double T /*in K*/, double Pressure /*in SI*/, double E /*in SI*/);
	Manager(UInt_t RandomSeed = 42);
	void LoopSimulation(void);
	void WriteHistory(std::string root_fname);
	void Test(void);
};

#endif 
