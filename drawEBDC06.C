void drawEBDC06(){
	
	gROOT->LoadMacro("sPhenixStyle.C");
	SetsPhenixStyle();

	TFile *f = new TFile("inputFiles/ebdc06_0000_clusterizer_new.root","READ");

	TTree *t = (TTree*)f->Get("hitTree");

	double X = 0.0;
	double Y = 0.0;
	double ADC = 0.0;
	int layer = 0.0;


	double binw = 0.0054541539;

	//TH2D *h = new TH2D("h","Run11011 sector18 south hits, 1k events;#phi;Pad Row;ADC",37,2.8775927,3.0255927,18,21.5,39.5);
	TH2D *h = new TH2D("h",";#phi;Pad Row;ADC",40,-0.42579939,-0.26579939,19,21.5,40.5);

	cout << "center " << h->GetXaxis()->GetBinCenter(21) << endl;

	t->SetBranchAddress("X",&X);
	t->SetBranchAddress("Y",&Y);
	t->SetBranchAddress("ADC",&ADC);
	t->SetBranchAddress("layer",&layer);

	for(int i=0; i<t->GetEntries(); i++){
		t->GetEntry(i);
		if(layer < 23 || layer > 38) continue;
		double phiV = atan2(Y,X);
		//if(phiV < 0) phiV += 2*TMath::Pi();
		if(phiV > -0.42179939) h->Fill(phiV,1.0*layer,ADC);
	}

	for(int i=1; i<=40; i++){
		for(int j=1; j<=19; j++){
			//if(h->GetBinContent(i,j) < (0.2e6)) h->SetBinContent(i,j,0.0);
		}
	}


	TFile *ft = new TFile("inputFiles/truthSectors.root","READ");

	TTree *tree = (TTree*)ft->Get("truthTree");

	unsigned int id = 0;
	TVector3 *vec = nullptr;

	tree->SetBranchAddress("truthID",&id);
	tree->SetBranchAddress("truthPos",&vec);

	vector<double> row, phi;

	for(int i=0; i<tree->GetEntries(); i++){
		tree->GetEntry(i);

		int pet = id/10000;

		if(pet < 18) continue;

		int rowV = (id - (pet*10000))/100;

		//cout << "rowV: " << rowV << endl;

		double rowPhi = atan2(vec->Y(),vec->X());
		//if(rowPhi < 0) rowPhi += 2*TMath::Pi();
		rowPhi += TMath::Pi()/18 - 0.005;

		//std::cout << "rowPhi: " << rowPhi << endl;

		if(rowV >= 16 && rowV <= 22){
			//cout << "good row" << endl;
			//if(rowPhi > 2.885 && rowPhi < 3.02){
			if(rowPhi > -0.42 && rowPhi < -0.27){
				//cout << "filling vectors" << endl;
				row.push_back(23.5 + 2*(rowV-16));
				phi.push_back(rowPhi);
			}
		}

	}

	gStyle->SetOptStat(0);

	TGraph *gr = new TGraph(phi.size(),&phi[0],&row[0]);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.5);
	gr->SetMarkerColor(kRed);
	gr->GetXaxis()->SetTitle("#phi");
	//gr->GetXaxis()->SetTitleOffset(1.2f);	
	gr->GetYaxis()->SetTitle("Stripe Row");
	//gr->GetYaxis()->SetTitleOffset(1.2f);
	gr->SetTitle("");

	TLatex sPHENIX;
	sPHENIX.SetTextFont(42);
	sPHENIX.SetTextAlign(12);
	sPHENIX.SetTextSize(0.035);

	TLatex TPC;
	TPC.SetTextSize(42);
	TPC.SetTextAlign(12);
	TPC.SetTextSize(0.035);

	TLatex laser;
	laser.SetTextSize(42);
	laser.SetTextAlign(12);
	laser.SetTextSize(0.035);

	TLatex date;
	date.SetTextSize(35);
	date.SetTextAlign(12);
	date.SetTextSize(0.035);

	TCanvas *c1 = new TCanvas("c1","",1000,1000);
	c1->SetTopMargin(0.055);
	c1->SetRightMargin(0.2f);
	c1->SetLeftMargin(0.15);
	//cout << "label size: " << h->GetXaxis()->GetLabelSize()  << endl;
	h->GetXaxis()->SetLabelSize(0.04);
	h->GetZaxis()->SetMaxDigits(2);
	h->GetZaxis()->SetTitleOffset(1.4);
	h->Draw("COLZ");
	//gr->Draw("PSAME");
	sPHENIX.DrawLatexNDC(0.175,0.915,"#bf{#it{sPHENIX}} Preliminary");
	laser.DrawLatexNDC(0.175,0.885,"Diffuse Laser Test");
	TPC.DrawLatexNDC(0.175,0.85,"TPC hits sector 23 R2 south");
	date.DrawLatexNDC(0.67,0.97,"#it{08/30/2023}");
	//Run11011 sector18 south hits, 1k events



	//TPaveText *hTitle = new TPaveText(0.25,0.92,0.75,0.98,"NDC");
	//hTitle->SetFillStyle(0);
	//hTitle->SetBorderSize(0);
	//hTitle->SetTextFont(42);
	//hTitle->SetTextSize(0.05);
	//hTitle->AddText("Run11011 sector18 R2 south hits, 1k events");
	//hTitle->AddText("Run11011 sector23 R2 south hits, 1k events");
	//hTitle->Draw("same");
	
	c1->SaveAs("Run11011_ebdc06_hits_truth_new.svg");
	c1->SaveAs("Run11011_ebdc06_hits_truth_new.png");

}
