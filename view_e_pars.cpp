 {
	TH1D* histE = new TH1D ("EnergyC [eV]","EnergyC [eV]",500,0, 15);
	TH1D* hist_dE = new TH1D ("Energy increase in resonance [eV]","Energy increase in resonance [eV]",300,-0.3, 0.3);
	TH1D* hist_dE_abs = new TH1D ("Energy abs increase in resonance [eV]","Energy abs increase in resonance [eV]",300,0.0, 0.3);
	TH1D* hist_Eloss = new TH1D ("Energy loss in resonance [eV]","Energy loss in resonance [eV]",300,0, 0.001);
	TH1D* hist_dT = new TH1D ("dTime [s]","dTime [s]",300,4e-16, 1e-12);
	TH1D* hist_dx = new TH1D ("dX [m]","dX [m]",300,-1e-6, 4e-6);
	TH1D* hist_V_drift = new TH1D ("Drift velocity [m/s]","Drift velocity [m/s]",300,0, 1e4);
	TH1D* hist_T_drift = new TH1D ("Drift time [s]","Drift time [s]",300,0, 2e-6);
	TH1D* histEfinalAvr = new TH1D ("Energy final avr","Energy final avr",300,0, 20);
	TH1D* histEfinal = new TH1D ("Energy final collision","Energy final collision",300,0, 15);
	
	int DEF_W = 900, DEF_H = 700;
	double E_at_time = 1e-11;
	double dt = 6e-12;
	double DRIFT_DISTANCE = 3e-3;
	std::string fname1("Output/v13.0/eData_0.3Td.root");
	std::string fname2("Output/v/eData_3Td_1.root");
	
	double En_start;
	double En_collision;
	double En_finish;
	//bool velocity_start; //1 - along z. 0 - against
	//bool velocity_finish; //1 - along z. 0 - against
	double En_avr;
	double pos_start;
	double pos_finish;
	double delta_x;
	double time_start;
	double delta_time;
	double delta_time_full;
	short process; //enum ProcessType : short {None = 0, Elastic = 1, Resonance = 2, Overflow = 3};
	
	long int num_of_elastic = 0, num_of_resonance=0, num_of_overflow = 0, num_of_electrons =0;
	long int num_of_ions = 0, num_of_excitations=0;
	for (int a = 0; a<2;++a) {
	    TFile * file = 0;
	    switch (a) 
	    {
		case 0: {
		    file = new TFile (fname1.c_str());
		    break;
		}
		case 1: {
		    file = new TFile (fname2.c_str());
		    break;
		}
	    }
	    if (0==file){
		std::cout<<"File not opened"<<std::endl;
		break;
	    }
	    if (!file->IsOpen()){
		std::cout<<"File not opened"<<std::endl;
		break;
	    }
	    TTree * tree = (TTree*) file->Get("ElectronHistory");
	    
	    tree->SetBranchAddress("energy_initial", &En_start);
	    tree->SetBranchAddress("energy_final", &En_finish);
	    tree->SetBranchAddress("energy_coll", &En_collision);
	    tree->SetBranchAddress("energy_average", &En_avr);
	    
	    //tree->SetBranchAddress("velocity_initial", &velocity_start);
	    //tree->SetBranchAddress("velocity_final", &velocity_finish);
	    tree->SetBranchAddress("time_initial", &time_start);
	    tree->SetBranchAddress("time_delta", &delta_time);
	    tree->SetBranchAddress("time_delta_full", &delta_time_full);
	    tree->SetBranchAddress("process_type", &process);
	    
	    tree->SetBranchAddress("position_initial",&pos_start);
	    tree->SetBranchAddress("position_final",&pos_finish);
	    tree->SetBranchAddress("position_delta",&delta_x);
	    
	    unsigned long int _end_ = tree->GetEntries();
	    for (unsigned long int i=0;i!=_end_;++i){
		    tree->GetEntry(i);
		    histE->Fill(std::fabs(En_collision));
		    hist_dT->Fill(delta_time_full);
		    hist_dx->Fill(delta_x);
		    
		    //if ((std::fabs(En_collision)<11.5)&&((std::fabs(En_collision)>10.8))) {
		    //if ((1==process)||(2==process)) {
			{
			hist_dE->Fill(std::fabs(En_collision) - std::fabs(En_start));
			hist_dE_abs->Fill(std::fabs(En_collision - En_start));
			hist_Eloss->Fill(std::fabs(En_collision)- std::fabs(En_finish));
		    }

		    switch (process) {
		      case 0:{
			  ++num_of_elastic;
			  break;
		      }
		      case 1:
		      case 2:{
			  ++num_of_resonance;
			  break;
		      }
		      case 3: {
			  ++num_of_ions;
			  break;
		      }
		      case -1:{
			  break;
		      }
		      case -2:{
			  ++num_of_overflow;
			  break;
		      }
		      default: {
			  ++num_of_excitations;
		      }
		    }
		    //if (!((pos_finish<DRIFT_DISTANCE)&&(pos_finish>-10*DRIFT_DISTANCE))) {
		    if (!(pos_finish<DRIFT_DISTANCE)) {
			++num_of_electrons;
			hist_V_drift->Fill(pos_finish/(time_start+delta_time_full));
			hist_T_drift->Fill(time_start+delta_time_full);
			histEfinalAvr->Fill(En_avr);
			histEfinal->Fill(std::fabs(En_collision));
		    }
	    }
	}
	TCanvas *c_ = new TCanvas ("e energy before collision", "e energy before collision", DEF_W, DEF_H);
	histE->Draw();
	TCanvas *c_2 = new TCanvas ("Delta time", "Delta time", DEF_W, DEF_H);
	hist_dT->Draw();
	TCanvas *c_3 = new TCanvas ("Drift velocity", "Drift velocity", DEF_W, DEF_H);
	hist_V_drift->Draw();
	TCanvas *c_3_5 = new TCanvas ("Drift time", "Drift time", DEF_W, DEF_H);
	hist_T_drift->Draw();
	TCanvas *c_4 = new TCanvas ("Delta x", "Delta x", DEF_W, DEF_H);
	hist_dx->Draw();
	TCanvas *c_5 = new TCanvas ("e energy after drift average", "e energy after drift average", DEF_W, DEF_H);
	histEfinalAvr->Draw();
	TCanvas *c_6 = new TCanvas ("e energy after drift before collision", "e energy after drift before collision", DEF_W, DEF_H);
	histEfinal->Draw();
	TCanvas *c_7 = new TCanvas ("e delta E in resonance [eV]", "e delta E in resonance [eV]", DEF_W, DEF_H);
	hist_dE->Draw();
	TCanvas *c_8 = new TCanvas ("e energy loss in resonance", "e energy loss in resonance", DEF_W, DEF_H);
	hist_Eloss->Draw();
	TCanvas *c_9 = new TCanvas ("e delta E abs in resonance [eV]", "e delta E abs in resonance [eV]", DEF_W, DEF_H);
	hist_dE_abs->Draw();
	
	TCanvas *c_12 = new TCanvas ("c Ec(t)", "c Ec(t)", DEF_W, DEF_H);
	TGraph *gr_Ec_t = new TGraph ();
	TCanvas *c_13 = new TCanvas ("c xf(t)", "c xf(t)", DEF_W, DEF_H);
	TGraph *gr_x_t = new TGraph ();
	
	std::cout<<"Drift distance: "<<DRIFT_DISTANCE<<std::endl;
	std::cout<<"Num of elastic: "<<num_of_elastic<<std::endl;
	std::cout<<"Num of resonance: "<<num_of_resonance<<std::endl;
	std::cout<<"Num of overflow: "<<num_of_overflow<<std::endl;
	std::cout<<"Num of ionizations: "<<num_of_ions<<std::endl;
	std::cout<<"Num of excitations: "<<num_of_excitations<<std::endl;
	std::cout<<"Num of electrons: "<<num_of_electrons<<std::endl;
	
}

//uses canvases and graphs defined above
void Plot_e (int n_of_e) {
    //std::cout<<"Plot_e("<<n_of_e<<")"<<std::endl;
	
    double En_start;
    double En_collision;
    double En_finish;
    //bool velocity_start; //1 - along z. 0 - against
    //bool velocity_finish; //1 - along z. 0 - against
    double En_avr;
    double pos_start;
    double pos_finish;
    double delta_x;
    double time_start;
    double delta_time;
    double delta_time_full;
    short process; //enum ProcessType : short {None = 0, Elastic = 1, Resonance = 2, Overflow = 3};
	
    long int num_of_electrons =0; 
    std::vector<double> Ec_, t_, x_; //for plotting Ec(t) 
    for (int a = 0; (a<2)&&(num_of_electrons<=n_of_e); ++a) {
	TFile * file = 0;
	switch (a) {
	case 0: {
	    file = new TFile (fname1.c_str());
	    break;
	}
	case 1: {
	    file = new TFile (fname2.c_str());
	    break;
	}
	}
	if (0==file){
	    std::cout<<"File not opened"<<std::endl;
	    break;
	}
	if (!file->IsOpen()){
	    std::cout<<"File not opened"<<std::endl;
	    break;
	}
	TTree * tree = (TTree*) file->Get("ElectronHistory");
	    
	tree->SetBranchAddress("energy_initial", &En_start);
	tree->SetBranchAddress("energy_final", &En_finish);
	tree->SetBranchAddress("energy_coll", &En_collision);
	tree->SetBranchAddress("energy_average", &En_avr);
	    
	//tree->SetBranchAddress("velocity_initial", &velocity_start);
	//tree->SetBranchAddress("velocity_final", &velocity_finish);
	tree->SetBranchAddress("time_initial", &time_start);
	tree->SetBranchAddress("time_delta", &delta_time);
	tree->SetBranchAddress("time_delta_full", &delta_time_full);
	tree->SetBranchAddress("process_type", &process);
	    
	tree->SetBranchAddress("position_initial",&pos_start);
	tree->SetBranchAddress("position_final",&pos_finish);
	tree->SetBranchAddress("position_delta",&delta_x);
	    
	unsigned long int _end_ = tree->GetEntries();
	for (unsigned long int i=0;(i!=_end_)&&(num_of_electrons<=n_of_e);++i){
	    tree->GetEntry(i);
	    if (n_of_e==num_of_electrons) {
		t_.push_back(time_start);
		Ec_.push_back(std::fabs(En_collision));
		x_.push_back(pos_finish);
	    }
	    //if (!((pos_finish<DRIFT_DISTANCE)&&(pos_finish>-10*DRIFT_DISTANCE))) {
	    if (!(pos_finish<DRIFT_DISTANCE)) {
		++num_of_electrons;
	    } 
	}
    }
    std::cout<<"Size: "<<t_.size()<<std::endl;
    gr_Ec_t->Set(t_.size());
    gr_x_t->Set(t_.size());
    for (int i=0, end_ = t_.size(); i!=end_; ++i) {
	gr_Ec_t->SetPoint(i, t_[i], Ec_[i]);
	gr_x_t->SetPoint(i, t_[i], x_[i]);
    }
    c_12->cd();
    gr_Ec_t->Draw("ALP");
    //c_12->Update();

    c_13->cd();
    gr_x_t->Draw("ALP");
    //c_13->Update();
}
 


