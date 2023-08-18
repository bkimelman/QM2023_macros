void drawTruthPatternXY(bool pos=true){
	
	gROOT->LoadMacro("sPhenixStyle.C");
	SetsPhenixStyle();

	TFile *f = new TFile("inputFiles/truthSectors.root","READ");

	TTree *tree = (TTree*)f->Get("truthTree");

	unsigned int id = 0;
	TVector3 *vec = nullptr;

	tree->SetBranchAddress("truthID",&id);
	tree->SetBranchAddress("truthPos",&vec);

	vector<double> X1, X2, X3, Y1, Y2, Y3;

	for(int i=0; i<tree->GetEntries(); i++){
		tree->GetEntry(i);

		if( (pos && id >= 180000) || (!pos && id < 180000 ) ) continue;

		int sec = id/10000;
		int row = (id - (sec*10000))/100;

		if(row<16){
			X1.push_back(vec->X());
			Y1.push_back(vec->Y());
		}else if(row>23){
			X3.push_back(vec->X());
			Y3.push_back(vec->Y());
		}else{
			X2.push_back(vec->X());
			Y2.push_back(vec->Y());
		}
	}

	TGraph *g1 = new TGraph(X1.size(),&X1[0],&Y1[0]);
	g1->SetMarkerStyle(20);
	g1->SetMarkerSize(0.5);
	g1->SetMarkerColor(kBlue);

	TGraph *g2 = new TGraph(X2.size(),&X2[0],&Y2[0]);
	g2->SetMarkerStyle(20);
	g2->SetMarkerSize(0.5);
	g2->SetMarkerColor(kRed);

	TGraph *g3 = new TGraph(X3.size(),&X3[0],&Y3[0]);
	g3->SetMarkerStyle(20);
	g3->SetMarkerSize(0.5);
	g3->SetMarkerColor(kGreen+2);
	g3->GetXaxis()->SetTitle("X [cm]");
	//g3->GetXaxis()->SetTitleOffset(1.2f);	
	g3->GetYaxis()->SetTitle("Y [cm]");
	//g3->GetYaxis()->SetTitleOffset(1.2f);
	g3->SetTitle("");

	TCanvas *c1 = new TCanvas("c1","",1000,1000);
	g3->Draw("AP");
	g2->Draw("PSAME");
	g1->Draw("PSAME");

	c1->SaveAs("TruthPositions.svg");

}
