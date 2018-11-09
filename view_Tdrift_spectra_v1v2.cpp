 {
	double TIME = 1.8e-8;
	int NBIN = 300;
	TH1D* histT_3_0_1 = new TH1D ("Drigt time 3.0 Td [s]","Drigt time 3.0 Td [s]", NBIN, 0, TIME);
	TH1D* histT_3_0_2 = new TH1D ("Drigt time 3.0 Td [s] Pb=0.31","Drigt time 3.0 Td [s] Pb=0.31", NBIN, 0, TIME);
	histT_3_0_1->SetStats(false);
	histT_3_0_2->SetStats(false);
	
	double DRIFT_DISTANCE = 1e-5;
	std::string fname_3_0_1("Output/v1/eData_3Td");
	std::string fname_3_0_2("Output/v2/eData_3Td");
	
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
	for (int nhist = 0; nhist<2;++nhist) {
	    TH1D* histT = NULL;
	    std::string fname;
	    switch (nhist) 
	    {
		case 0: {
		    histT = histT_3_0_1;
		    fname = fname_3_0_1;
		    break;
		}
		case 1: {
		    histT = histT_3_0_2;
		    fname = fname_3_0_2;
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
	TLegend *legend = new TLegend( 0.5, 0.7, 0.9, 0.9);
	//legend->SetHeader("");
	legend->SetMargin(0.4);
	TH2F* frame = new TH2F( "frame", "Drift time spectra", 500, 0, TIME, 500, 0, max_val);
	frame->GetXaxis()->SetTitle("t [s]");
	frame->GetYaxis()->SetTitle("");
	frame->Draw();
	
	histT_3_0_1->SetLineWidth(2);
	histT_3_0_1->SetLineColor(2);
	histT_3_0_1->Draw("csame");
	histT_3_0_2->SetLineWidth(2);
	histT_3_0_2->SetLineColor(12);
	histT_3_0_2->Draw("csame");
	
	legend->AddEntry(histT_3_0_1, (std::string("3.0 Td, <t>[ns]=")+std::to_string(1e9*histT_3_0_1->GetMean())).c_str(), "l");
	legend->AddEntry(histT_3_0_2, (std::string("3.0 Td, Pb=0.31, <t>[ns]=")+std::to_string(1e9*histT_3_0_2->GetMean())).c_str(), "l");
	
	frame->Draw("sameaxis");
	legend->Draw("same");
	c_->Update();
}



