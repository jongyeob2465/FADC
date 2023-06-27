//this file is for LS accumulation 

vector<double> computeConvolution(const vector<double>& a, const vector<double>& b){
	vector<double> result(a.size() + b.size() - 1, 0);

	for (size_t i = 0; i < a.size(); i++){
		for (size_t j = 0; j < b.size(); j++)
			result[i + j] += a[i] * b[j];
	}

	return result;
}

double getMaxConvolution(const vector<double>& a, const vector<double>& b){
	vector<double> result = computeConvolution(a,b);
	double max = *max_element(result.begin(),result.end());

	return max ;
}

void template_LS(){
	TH1::AddDirectory(kFALSE);

	const int hist_point = 4092;
	int fadc0[hist_point] , cnt ;
	double ped, Qtotal, Qtail;
	double ped_count, peak_bin;
	double value ;

	TChain * chain = new TChain("sadc");
	chain->Add("/home/s2022220206/Final_exam/file/FADC_data_1.root");
	chain->SetBranchAddress("fadc0",&fadc0);

	int total=chain->GetEntries();
	vector<double> Qtotal_list(total,0) ;

	TH1D * his = new TH1D("his","his",hist_point,0,hist_point); // pulse shape  

	vector<vector<double>> vecOfVecs(total,vector<double>(40));
	
	// event loop
	for(int i = 0; i < total ; i++){

		vector<double> tempList(40);
		chain->GetEntry(i);
		ped = 0 ; Qtotal = 0 ;
		ped_count = 0; peak_bin = 0;
		cnt = 0 ;		
		value = 0 ;

		// pedestal loop
		for(int jj = 0 ; jj < hist_point; jj++){
			his->SetBinContent(jj+1,fadc0[jj]);

			if (jj >=0 && jj < 1100){
			ped += fadc0[jj];
			ped_count++;
			}
		}
		ped /= (double)ped_count; // normalize pedestal
		peak_bin = his->GetMaximumBin();

		//Qtotal loop
		for(int jj = peak_bin - 10 ; jj < peak_bin + 30 ; jj++)	Qtotal += fadc0[jj] - ped;		
		Qtotal_list[i] = Qtotal ;

		// calculate the bin content using pedestal
		for(int jj = peak_bin - 10 ; jj < peak_bin + 30 ; jj++){
			value = (double)(fadc0[jj] - ped) ;
			tempList[cnt] = value ; 
			cnt ++ ;
		}	

		vecOfVecs[i] = tempList;
	}

	TH1D * total_charge = new TH1D("total_charge","total_charge",100,0,20000); // total charge distribution
	TH1D * max_conv = new TH1D("max_conv","max_conv",150,0.,0.15); // total charge distribution
	TH1D * LS_accu = new TH1D("LS_accu","LS_accu",40,0,40); // pdf
	vector<double> temp(40);
	double sum ;

	//make template
	for (int j = 0 ; j < 40 ; j ++){
		sum = 0 ;		
		for (int i = 0 ; i < total ; i++) sum += vecOfVecs[i][j];
		LS_accu->SetBinContent(j+1,sum);
	}
	LS_accu->Scale(1.0/LS_accu->Integral());

	for (int i = 0 ; i < 40 ; i++) temp[i] = LS_accu->GetBinContent(i+1);

	vector<double> norm(40,0);
	vector<double> conv_list(total,0) ;

	// calculate convolution
	for (int i = 0 ; i < total ; i++){
		for (int j = 0 ; j < 40 ; j++) norm[j] = vecOfVecs[i][j]/Qtotal_list[i] ;
		conv_list[i] = getMaxConvolution(temp,norm);
//		cout << conv_list[i] << endl;
	}
		
	// fill histogram
	for (int i = 0 ; i < total ; i++){
		total_charge->Fill(Qtotal_list[i]);
		max_conv->Fill(conv_list[i]);
	}

	// skimming
	for (int i = 0 ; i < total ; i++){
		if ((conv_list[i] < 0.05) || (Qtotal_list[i] < 1400)){
//			cout << "max convolution value : " <<	conv_list[i] << " , Qtotal : " << Qtotal_list[i] << ", event number : " <<  i << endl;
			for (int j = 0 ; j < 40 ; j++) vecOfVecs[i][j] = 0 ;
		}
	}

	//accumulate
	for (int j = 0 ; j < 40 ; j ++){
		sum = 0 ;		
		for (int i = 0 ; i < total ; i++) sum += vecOfVecs[i][j];
		LS_accu->SetBinContent(j+1,sum);
	}
	LS_accu->Scale(1.0/LS_accu->Integral());

	total_charge->GetXaxis()->SetTitle("totalQ");
	total_charge->GetYaxis()->SetTitle("number of events");

	max_conv->GetXaxis()->SetTitle("max convolution value");
	max_conv->GetYaxis()->SetTitle("number of events");


	TFile o("/home/s2022220206/Final_exam/LS/LS_accu.root","RECREATE");
	LS_accu->Write("hist");
	total_charge->Write("totalQ");
	max_conv->Write("max_conv");
	o.Close();
}

