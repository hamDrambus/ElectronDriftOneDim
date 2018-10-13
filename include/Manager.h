#ifndef MANAGER_H_
#define MANAGER_H_

#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>

#include "global_definitions.h"
#include "argon_cross.h"

struct Event
{
	double En_start;
	bool velocity_start; //1 - along z. 0 - against
	double En_finish;
	bool velocity_finish; //1 - along z. 0 - against
	double En_avr;
	long double pos_start;
	long double pos_finish;
	long double time_start;
	long double delta_time;
	enum ProcessType {None = 0, Elastic = 1, Resonance = 2, Overflow = 3}
	process;
	/*enum StateType { Start = 0, Collision = 1, Resonance = 2, Finished = 3}
	state;*/
};

class Manager
{
protected:
	TRandom * random_generator_;
	TTree * sim_data_;
	Event event_;
	//Event current_event;
	void DoStepLength (Event &event);
	void DoScattering(Event &event);
	void PostStepAction(Event &event);
	void DoGotoNext(Event &event);
	bool is_ready_;
	double Concentration_;
	double eField_;
	double Coefficient_;
	long double XS_integral(double from, double to);
public:
	void Initialize(Event &event);
	void Clear(void);
	void DoStep(Event &event);
	bool IsFinished(Event &event);
	void SetParameters(double Concetr /*in SI*/, double E /*in SI*/);
	void SetParameters(double T /*in K*/, double Pressure /*in SI*/, double E /*in SI*/);
	Manager(UInt_t RandomSeed = 42); //programm instance is required for concurrent launch of this code (without multithreading)
	void LoopSimulation(void);
	void WriteHistory(std::string root_fname);
};

#endif 
