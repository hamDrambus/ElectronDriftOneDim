//REQUIRES .L Tdrift_fitfunc.cpp
 {
	double TIME = 2e-8;
	int NBIN = 300;
	TH1D* histT_3_0 = new TH1D ("Drigt time 3.0 Td [s]","Drigt time 3.0 Td [s]", NBIN, 0, TIME);
	histT_3_0->SetStats(false);
	
	double DRIFT_DISTANCE = 1e-5;
	std::string fname_3_0("Output/v1/eData_1Td");
	
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
	
	double max_val = 0;
	for (int nhist = 0; nhist<1;++nhist) {
	    TH1D* histT = NULL;
	    std::string fname;
	    switch (nhist) 
	    {
		case 0: {
		    histT = histT_3_0;
		    fname = fname_3_0;
		    break;
		}
	    }
	    if (0==histT){
		std::cout<<"No histogram"<<std::endl;
		continue;
	    }
	    for (int ver=0; ver<2; ++ver) {
		TFile * file = 0;
		std::string filename = fname;
		switch (ver) 
		{
		    case 0: {
			break;
		    }
		    case 1: {
			filename = filename + "_1";
			break;
		    }
		}
		filename= filename+".root";
		file = new TFile (filename.c_str());
		if (0==file ){
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
		    if (!(pos_finish<DRIFT_DISTANCE)) {
			histT->Fill(time_start+delta_time_full);
		    }
		}
		for (int bin = 1, bin_end = histT->GetNbinsX(); bin!=bin_end; ++bin) {
		    max_val = std::max(max_val, (double) histT->GetBinContent(bin));
		}
	    }
	}
	max_val*=1.1;
	TCanvas *c_ = new TCanvas ("Drift time spectra_", "Drift time spectra_");
	gStyle->SetOptStat("");
	TLegend *legend = new TLegend( 0.3, 0.5, 0.9, 0.9);
	//legend->SetHeader("");
	legend->SetMargin(0.4);
	TH2F* frame = new TH2F( "frame", "Drift time spectrun 1.0 Td", 500, 0, TIME, 500, 0, max_val);
	frame->GetXaxis()->SetTitle("t [s]");
	frame->GetYaxis()->SetTitle("");
	frame->Draw();
	
	histT_3_0->SetLineWidth(2);
	histT_3_0->SetLineColor(2);
	histT_3_0->Draw("csame");
	
	TF1 *func1 = new TF1("fit1",fit_f_1,0,TIME, 3);
	func1->SetParLimits(0, 0, 1e-10); //signal start time
	func1->SetParLimits(1, 1, 1000); //Amplitude at signal start time
	func1->SetParLimits(2, 2e-10, 8e-9); //time const
	TF1 *func2 = new TF1("fit2",fit_f_2,0,TIME, 5);
	func2->SetParLimits(0, 0, 1e-10); //signal start time
	func2->SetParLimits(1, 1, 1000); //Amplitude at signal start time
	func2->SetParLimits(2, 2e-10, 6e-9); //time const
	func2->SetParLimits(3, 1, 1000); //Amplitude at signal start time
	func2->SetParLimits(4, 1e-9, 5e-8); //time const
	
	histT_3_0->Fit(func1, "NV");
	histT_3_0->Fit(func2, "NV");
	
	func1->SetLineColor(3);
	func2->SetLineColor(4);
	func1->Draw("same");
	func2->Draw("same");
	
	legend->AddEntry(histT_3_0, (std::string("1.0 Td, <t>[ns]=")+std::to_string(1e9*histT_3_0->GetMean())).c_str(), "l");
	legend->AddEntry(func1, (std::string("fit, tau[ns]=")+std::to_string(1e9*func1->GetParameter(2))).c_str(), "l");
	legend->AddEntry(func2, (std::string("fit, tau1[ns]=")+std::to_string(1e9*func2->GetParameter(2)) + " tau2[ns]="
	  + std::to_string(1e9*func2->GetParameter(4))).c_str(), "l");
	
	frame->Draw("sameaxis");
	legend->Draw("same");
	c_->Update();
}



