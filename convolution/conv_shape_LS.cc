// LS data convolution
vector<long double> computeConvolution(const vector<long double>& a, const vector<long double>& b) {
    vector<long double> result(a.size() + b.size() - 1, 0);

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            result[i + j] += a[i] * b[j];
        }
    }

    return result;
}

void conv_shape_LS(){

	TH1::AddDirectory(kFALSE);

	const int hist_point = 4092;
	int fadc0[hist_point];
	int cnt = 0 ; 

	TChain * chain = new TChain("sadc");
	chain->Add("/home/s2022220206/Final_exam/file/FADC_data_1.root");
	chain->SetBranchAddress("fadc0",&fadc0);

	int total=chain->GetEntries();

	TH1D * his = new TH1D("his","his",hist_point,0,hist_point);

	double ped, Qtotal, Qtail;
	double ped_count;
	double peak_bin;
	long double value ;

	vector<vector<long double>> vecOfVecs;	

	//event loop  vector of vector 
	for(int i = 0; i < total ; i++){
		vector<long double> tempList;

		chain->GetEntry(i);
		ped = 0 ; Qtotal = 0 ;
		ped_count = 0; peak_bin = 0;
		cnt = 0 ;

		// pedestal loop
		for(int jj = 0 ; jj < hist_point; jj++){
			his->SetBinContent(jj+1,fadc0[jj]);
			if (jj >=0 && jj < 1100){
			ped += fadc0[jj];
			ped_count++;
			}
		}
		ped /= (double)ped_count;
		peak_bin = his->GetMaximumBin();

		//Qtotal loop
		for(int jj = peak_bin - 10 ; jj < peak_bin + 30 ; jj++){
			Qtotal += fadc0[jj] - ped;
		}

		//cout << "pedestal : " << ped << ", Q total : " << Qtotal << endl;
		// calculate the bin content using pedestal and Qtotal
		for(int jj = peak_bin - 10 ; jj < peak_bin + 30 ; jj++){
			value = (long double)((fadc0[jj] - ped)/Qtotal ) ;
			tempList.push_back(value); // listoflist input number
		
		}
		
		vecOfVecs.push_back(tempList);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
	TFile compare("/home/s2022220206/Final_exam/LS/LS_accu.root");
	TH1D *pdf = (TH1D*)compare.Get("hist");

	vector<long double> pdf_content(40) ;

	for (int i = 0 ; i < 40 ; ++i) pdf_content[i] = pdf->GetBinContent(i);

///////////////////////////////////////////////////////////////////////////////////////////////////

	TCanvas *c = new TCanvas("canvas","Convolution graph",800,600);
	TGraph * grp = new TGraph();
//	TH1D *conv_hist = new TH1D("conv_hist","conv_hist",100,0,1000000);
	
	long double max ;
	Int_t binmax ; 
	
	for(int event = 0 ; event < total ; event++){	
			
		vector<long double>result = computeConvolution(pdf_content , vecOfVecs[event]);
//		max = *max_element(result.begin(),result.end());
//		conv_hist->Fill(max);
		for(int i = 0 ; i < 50 ; i++){
//			if (result[i] == 0) cout << event << endl;
			grp->SetPoint(i,i,result[i]);
		}

		grp->GetYaxis()->SetRangeUser(-0.001,0.15);
		grp->SetMarkerStyle(7);
		grp->GetXaxis()->SetTitle("bin numb");
		grp->GetYaxis()->SetTitle("Convolution");

		c->cd();
		grp->Draw("AP");
		c->Update();
		c->Modified();
		gSystem->ProcessEvents();
	}

//	conv_hist->Draw("hist");
		
}

