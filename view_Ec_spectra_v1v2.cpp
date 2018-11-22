 {
	
	TH1D* histE_3_0_1 = new TH1D ("EnergyC 3.0 Td [eV] Ld=1e-6","EnergyC 3.0 Td [eV] Ld=1e-6",300,0, 14);
	TH1D* histE_3_0_2 = new TH1D ("EnergyC 3.0 Td [eV] Ld=1e-5","EnergyC 3.0 Td [eV] Ld=1e-5",300,0, 14);
	TH1D* histE_3_0_3 = new TH1D ("EnergyC 3.0 Td [eV] Ld=1e-4","EnergyC 3.0 Td [eV] Ld=1e-4",300,0, 14);
	TH1D* histE_3_0_4 = new TH1D ("EnergyC 3.0 Td [eV] Ld=1e-3","EnergyC 3.0 Td [eV] Ld=1e-3",300,0, 14);
	histE_3_0_1->SetStats(false);
	histE_3_0_2->SetStats(false);
	histE_3_0_3->SetStats(false);
	histE_3_0_4->SetStats(false);
	
	int DEF_W = 900, DEF_H = 700;
	std::vector<double> DRIFT_DISTANCE = {1e-6, 1e-5, 1e-4, 1e-3};
	std::string fname_3_0_1("Output/v8/eData_3Td");
	std::string fname_3_0_2("Output/v7/eData_3Td");
	std::string fname_3_0_3("Output/v9/eData_3Td");
	std::string fname_3_0_4("Output/v12/eData_3Td");
	
	double En_start;
	double En_collision;
	double En_finish;
	double En_avr;
	double pos_start;
	double pos_finish;
	double delta_x;
	double time_start;
	double delta_time;
	double delta_time_full;
	short process; //enum ProcessType : short {None = 0, Elastic = 1, Resonance = 2, Overflow = 3};

	
	double max_val = 0;
	for (int nhist = 0; nhist<4;++nhist) {
	    TH1D* histE = NULL;
	    std::string fname;
	    switch (nhist) 
	    {
		case 0: {
		    histE = histE_3_0_1;
		    fname = fname_3_0_1;
		    break;
		}
		case 1: {
		    histE = histE_3_0_2;
		    fname = fname_3_0_2;
		    break;
		}
		case 2: {
		    histE = histE_3_0_3;
		    fname = fname_3_0_3;
		    break;
		}
		case 3: {
		    histE = histE_3_0_4;
		    fname = fname_3_0_4;
		    break;
		}
	    }
	    if (0==histE){
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
		}
	    }
	    //normalization of spectrum; like in Buzulutzkov paper
	    double Norm =0;
	    for (int bin = 1, bin_end = histE->GetNbinsX(); bin!=bin_end; ++bin) {
		Norm+=histE->GetBinContent(bin)/**pow(histE->GetBinCenter(bin),0.5)*/*histE->GetBinWidth(bin);
	    }
	    for (int bin = 1, bin_end = histE->GetNbinsX(); bin!=bin_end; ++bin) {
		histE->SetBinContent(bin, histE->GetBinContent(bin)/(Norm));
		max_val = std::max(max_val, (double) histE->GetBinContent(bin));
	    }
	}
	max_val*=1.1;
	TCanvas *c_ = new TCanvas ("Collision E spectra_", "Collision E spectra_", DEF_W, DEF_H);
	gStyle->SetOptStat("");
	TLegend *legend = new TLegend( 0.55, 0.65, 0.9, 0.9);
	//legend->SetHeader("");
	legend->SetMargin(0.25);
	TH2F* frame = new TH2F( "frame", "Collision E spectra", 500, 0, 14, 500, 0, max_val);
	frame->GetXaxis()->SetTitle("Ee [eV]");
	frame->GetYaxis()->SetTitle("");
	frame->Draw();
	
	histE_3_0_1->SetLineWidth(3);
	histE_3_0_1->SetLineColor(2);
	histE_3_0_1->Draw("csame");
	histE_3_0_2->SetLineWidth(3);
	histE_3_0_2->SetLineColor(12);
	histE_3_0_2->Draw("csame");
	histE_3_0_3->SetLineWidth(3);
	histE_3_0_3->SetLineColor(4);
	histE_3_0_3->Draw("csame");
	histE_3_0_4->SetLineWidth(3);
	histE_3_0_4->SetLineColor(3);
	histE_3_0_4->Draw("csame");
	
	legend->AddEntry(histE_3_0_1, (std::string("3.0 Td, Ld=1e-6, <Ec>=")+std::to_string(histE_3_0_1->GetMean())).c_str(), "l");
	legend->AddEntry(histE_3_0_2, (std::string("3.0 Td, Ld=1e-5, <Ec>=")+std::to_string(histE_3_0_2->GetMean())).c_str(), "l");
	legend->AddEntry(histE_3_0_3, (std::string("3.0 Td, Ld=1e-4, <Ec>=")+std::to_string(histE_3_0_3->GetMean())).c_str(), "l");
	legend->AddEntry(histE_3_0_4, (std::string("3.0 Td, Ld=1e-3, <Ec>=")+std::to_string(histE_3_0_4->GetMean())).c_str(), "l");
	
	frame->Draw("sameaxis");
	legend->Draw("same");
	c_->Update();
}



