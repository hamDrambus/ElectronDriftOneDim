 {
	double EN_MIN_=0;
	double EN_MAX_=14;
	int NN = 300;
	TH1D* histE_0 = new TH1D ("EnergyC 7.0 Td [eV] normal res.","EnergyC 7.0 Td [eV] normal res.",NN,EN_MIN_, EN_MAX_);
	TH1D* histE_1 = new TH1D ("EnergyC 7.0 Td [eV] W=3e-2 eV","EnergyC 7.0 Td [eV] W=3e-2 eV",NN,EN_MIN_, EN_MAX_);
	TH1D* histE_2 = new TH1D ("EnergyC 7.0 Td [eV] W=3e-2 eV, Loss=50x","EnergyC 7.0 Td [eV] W=3e-2 eV, Loss=50x",NN,EN_MIN_, EN_MAX_);
	TH1D* histE_3 = new TH1D ("EnergyC 7.0 Td [eV] W=3e-2 eV, Loss=100x","EnergyC 7.0 Td [eV] W=3e-2 eV, Loss=100x",NN,EN_MIN_, EN_MAX_);
	TH1D* histE_4 = new TH1D ("EnergyC 7.0 Td [eV] W=3e-2 eV, Loss=150x","EnergyC 7.0 Td [eV] W=3e-2 eV, Loss=150x",NN,EN_MIN_, EN_MAX_);
	histE_0->SetStats(false);
	histE_1->SetStats(false);
	histE_2->SetStats(false);
	histE_3->SetStats(false);
	histE_4->SetStats(false);
	
	int DEF_W = 900, DEF_H = 700;
	std::vector<double> DRIFT_DISTANCE = {3e-3, 3e-3, 3e-3, 3e-3};
	std::string fname_0("Output/v13.1/eData_0.3Td");
	std::string fname_1("Output/v13.1/eData_1.0Td");
	std::string fname_2("Output/v13.1/eData_3.0Td");
	std::string fname_3("Output/v13.1/eData_7.0Td");
	std::string fname_4("Output/v15.5/eData_7.0Td");
	
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
		    histE = histE_0;
		    fname = fname_0;
		    break;
		}
		case 1: {
		    histE = histE_1;
		    fname = fname_1;
		    break;
		}
		case 2: {
		    histE = histE_2;
		    fname = fname_2;
		    break;
		}
		case 3: {
		    histE = histE_3;
		    fname = fname_3;
		    break;
		}
		case 4: {
		    histE = histE_4;
		    fname = fname_4;
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
	    for (int bin = 1, bin_end = histE->GetNbinsX()+1; bin!=bin_end; ++bin) {
		Norm+=histE->GetBinContent(bin)/**pow(histE->GetBinCenter(bin),0.5)*/*histE->GetBinWidth(bin);
	    }
	    for (int bin = 1, bin_end = histE->GetNbinsX()+1; bin!=bin_end; ++bin) {
		histE->SetBinContent(bin, histE->GetBinContent(bin)/(Norm));
		max_val = std::max(max_val, (double) histE->GetBinContent(bin));
	    }
	}
	max_val*=1.1;
	gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("");
	TCanvas *c_ = new TCanvas ("Collision E spectra_", "Collision E spectra_", DEF_W, DEF_H);
	c_->SetGrid();
	c_->SetTicks();
	TLegend *legend = new TLegend( 0.55, 0.65, 0.9, 0.9);
	//legend->SetHeader("");
	legend->SetMargin(0.25);
	TH2F* frame = new TH2F( "frame", "Collision E spectra 1D", 500, EN_MIN_, EN_MAX_, 500, 0, max_val);
	frame->GetXaxis()->SetTitle("Ee [eV]");
	frame->GetYaxis()->SetTitle("");
	frame->Draw();
	
	histE_0->SetLineWidth(2);
	histE_0->SetLineColor(2);
	histE_0->Draw("csame");
	histE_1->SetLineWidth(2);
	histE_1->SetLineColor(12);
	histE_1->Draw("csame");
	histE_2->SetLineWidth(2);
	histE_2->SetLineColor(4);
	histE_2->Draw("csame");
	histE_3->SetLineWidth(2);
	histE_3->SetLineColor(3);
	histE_3->Draw("csame");
	histE_4->SetLineWidth(2);
	//histE_4->SetLineColor(6);
	//histE_4->Draw("csame");
	
	legend->AddEntry(histE_0, (std::string("0.3 Td")).c_str(), "l");
	legend->AddEntry(histE_1, (std::string("1.0 Td")).c_str(), "l");
	legend->AddEntry(histE_2, (std::string("3.0 Td")).c_str(), "l");
	legend->AddEntry(histE_3, (std::string("7.0 Td")).c_str(), "l");
	//legend->AddEntry(histE_4, (std::string("7.0 Td, W=3e-2 eV, Loss=150x")).c_str(), "l");
	
	frame->Draw("sameaxis");
	legend->Draw("same");
	c_->Update();
}



