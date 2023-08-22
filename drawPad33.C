void drawPad33(){
	
	gROOT->LoadMacro("sPhenixStyle.C");
	SetsPhenixStyle();

	TFile *f11011 = new TFile("inputFiles/Run11011_ebdc07_clusterizer.root","READ");
	TFile *f11028 = new TFile("inputFiles/Run11028_ebdc07_clusterizer.root","READ");

	TTree *t11 = (TTree*)f11011->Get("hitTree");
	TTree *t28 = (TTree*)f11028->Get("hitTree");

	double X11 = 0.0;
	double Y11 = 0.0;
	double ADC11 = 0.0;
	int layer11 = 0.0;

	double X28 = 0.0;
	double Y28 = 0.0;
	double ADC28 = 0.0;
	int layer28 = 0.0;

	TH2D *h11 = new TH2D("h11","Run11011 sector19 south hits, 1k events;#phi;Pad Row;ADC",148,-2.92069942013426,-2.31528833584873,18,21.5,40.5);
	TH2D *h28 = new TH2D("h28","Run11028 sector19 south hits, 1k events;#phi;Pad Row;ADC",148,-2.92069942013426,-2.31528833584873,18,21.5,40.5);

	t11->SetBranchAddress("X",&X11);
	t11->SetBranchAddress("Y",&Y11);
	t11->SetBranchAddress("ADC",&ADC11);
	t11->SetBranchAddress("layer",&layer11);

	t28->SetBranchAddress("X",&X28);
	t28->SetBranchAddress("Y",&Y28);
	t28->SetBranchAddress("ADC",&ADC28);
	t28->SetBranchAddress("layer",&layer28);

	for(int i=0; i<t11->GetEntries(); i++){
		t11->GetEntry(i);
		if(layer11 < 23 || layer11 > 38) continue;
		h11->Fill(atan2(Y11,X11),1.0*layer11,ADC11);
	}

	for(int i=0; i<t28->GetEntries(); i++){
		t28->GetEntry(i);
		if(layer28 < 23 || layer28 > 38) continue;
		h28->Fill(atan2(Y28,X28),1.0*layer28,ADC28);
	}

	TFile *f = new TFile("inputFiles/truthSectors.root","READ");

	TTree *tree = (TTree*)f->Get("truthTree");

	unsigned int id = 0;
	TVector3 *vec = nullptr;

	tree->SetBranchAddress("truthID",&id);
	tree->SetBranchAddress("truthPos",&vec);

	vector<double> row, phi;

	for(int i=0; i<tree->GetEntries(); i++){
		tree->GetEntry(i);

		int pet = id/10000;

		if(pet != 33 && pet != 34) continue;

		int rowV = (id - (pet*10000))/100;

		double rowPhi = atan2(vec->Y(),vec->X());
		rowPhi -= TMath::Pi()/18;

		if(rowV >= 16 && rowV <= 23 && rowPhi < -2.35 && rowPhi > -2.9){
			row.push_back(23.5 + 2*(rowV-16));
			phi.push_back(rowPhi);
		}

	}

	gStyle->SetOptStat(0);

	TGraph *gr = new TGraph(phi.size(),&phi[0],&row[0]);
	gr->SetMarkerStyle(53);
	gr->SetMarkerSize(1.25);
	gr->GetXaxis()->SetTitle("#phi");
	//gr->GetXaxis()->SetTitleOffset(1.2f);	
	gr->GetYaxis()->SetTitle("Stripe Row");
	//gr->GetYaxis()->SetTitleOffset(1.2f);
	gr->SetTitle("");

	TLatex sPHENIX;
	sPHENIX.SetTextFont(42);
	sPHENIX.SetTextAlign(12);
	sPHENIX.SetTextSize(0.035);
	
	TCanvas *c1 = new TCanvas("c1","",1000,1000);
	c1->SetTopMargin(0.1);
	c1->SetRightMargin(0.2f);
	c1->SetLeftMargin(0.15);
	h11->GetZaxis()->SetMaxDigits(1);
	h11->GetZaxis()->SetTitleOffset(1.4);
	h11->Draw("COLZ");
	gr->Draw("PSAME");
	sPHENIX.DrawLatexNDC(0.175,0.85,"#bf{#it{sPHENIX}} Preliminary");

	TPaveText *h11Title = new TPaveText(0.25,0.92,0.75,0.98,"NDC");
	h11Title->SetFillStyle(0);
	h11Title->SetBorderSize(0);
	h11Title->SetTextFont(42);
	h11Title->SetTextSize(0.05);
	h11Title->AddText("Run11011 sector19 south hits, 1k events");
	h11Title->Draw("same");
	
	c1->SaveAs("Run11011_hits_truth.svg");
	c1->SaveAs("Run11011_hits_truth.png");

	c1->Clear();
	h28->GetZaxis()->SetMaxDigits(2);
	h28->GetZaxis()->SetTitleOffset(1.4);
	h28->Draw("COLZ");
	gr->Draw("PSAME");

	TPaveText *h28Title = new TPaveText(0.25,0.92,0.75,0.98,"NDC");
	h28Title->SetFillStyle(0);
	h28Title->SetBorderSize(0);
	h28Title->SetTextFont(42);
	h28Title->SetTextSize(0.05);
	h28Title->AddText("Run11028 sector19 south hits, 1k events");
        h28Title->Draw("same");
	sPHENIX.DrawLatexNDC(0.175,0.85,"#bf{#it{sPHENIX}} Preliminary");
	
	c1->SaveAs("Run11028_hits_truth.svg");
	c1->SaveAs("Run11028_hits_truth.png");

}
