// LS data likelihood calculation
void cal_LS(){

	TH1::AddDirectory(kFALSE);

	const int hist_point = 4092;
	int fadc0[hist_point];
	int cnt = 0 ;
	int PS_like = 0 , LS_like = 0 ;
	int charge_cut = 1500 ; 
	double content[40];
	double ped, Qtotal, Qtail;
	double ped_count;
	double peak_bin;
	double total_LS , total_PS ;

	TFile LS("LS/LS_accu.root") ; // LS pdf
	TH1D * LS_accu  = (TH1D*)LS.Get("hist") ;

	TFile PS("PS/PS_accu.root") ; // PS pdf
	TH1D * PS_accu = (TH1D*)PS.Get("hist") ;

	TChain * chain = new TChain("sadc");
	chain->Add("/home/s2022220206/Final_exam/file/FADC_data_1.root");
	chain->SetBranchAddress("fadc0",&fadc0);

	int total=chain->GetEntries();

	TH1D * his = new TH1D("his","his",hist_point,0,hist_point);
	TH1D * likelihood = new TH1D("likelihood","likelihood",200,-350,50);
	TH1D * likelihood_LS = new TH1D("likelihood_LS","likelihood_LS",1000,-50000,0);
	TH1D * likelihood_PS = new TH1D("likelihood_PS","likelihood_PS",1000,-50000,0);

	ofstream file1("LS_like.txt");
	ofstream file2("PS_like.txt");

	cout << "total " << total << endl;
	cout << "Qtotal cut : " << charge_cut << endl;

	//event loop
	for(int i = 0; i < total ; i++){

		chain->GetEntry(i);

		ped = 0 ; Qtotal = 0 ; Qtail = 0;
		ped_count = 0; peak_bin = 0;
		cnt = 0 ; total_LS = 0 ; total_PS = 0 ;
		
		//pedestal loop
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
			content[cnt] = fadc0[jj] - ped ;
			cnt++ ;	
		}
	
		if (Qtotal < charge_cut ) continue ;
		// likelihood calculation
		for (int jj = 0 ; jj < 40 ; jj++){
			 total_LS += content[jj] * log10( LS_accu->GetBinContent(jj+1));
			}

		for (int jj = 0 ; jj < 40 ; jj++){
			 total_PS += content[jj] * log10(PS_accu->GetBinContent(jj+1));
			}
		
		likelihood_LS->Fill(total_LS);
		likelihood_PS->Fill(total_PS);
		likelihood->Fill(total_PS - total_LS);
		if ((total_PS - total_LS) > 0) {
			file1 << i << endl;
			PS_like += 1 ;
		}
		else{
			file2 << i << endl;
			LS_like += 1 ;
		}
	}
	//event loop end
	cout << "Qtotal cut passed events : " << PS_like + LS_like << endl;
	cout << "PS like event number : " << PS_like << endl;
	cout << "LS like event number : " << LS_like << endl;
	
	file1.close();
	file2.close();

	TCanvas * c1 = new TCanvas("c1","c1");
	likelihood->GetXaxis()->SetTitle("log(PS) - log(LS)");
	likelihood->GetYaxis()->SetTitle("number of event");
	likelihood->Draw("hist");

	TFile file("result_LS.root","RECREATE");
	likelihood->Write();
//	likelihood_LS->Write();
//	likelihood_PS->Write();
	file.Close();

}
