void drawInputDist(){
	
	gROOT->LoadMacro("sPhenixStyle.C");
	SetsPhenixStyle();

	//TFile *sDist = new TFile("/sphenix/user/rcorliss/distortion_maps/2023.02/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX99528306_5.distortion_map.hist.root","READ");
	TFile *sDist = new TFile("inputFiles/staticInput.root","READ");

	TH3F *sRad3D = (TH3F*)sDist->Get("hIntDistortionR_negz");
	TH3F *sPhi3D = (TH3F*)sDist->Get("hIntDistortionP_negz");

	//TFile *fDist = new TFile("/sphenix/user/rcorliss/distortion_maps/2023.02/TimeOrderedDistortions.root","READ");
	TFile *fDist = new TFile("inputFiles/fluctInput.root","READ");

	TTree *timeTree = (TTree*)fDist->Get("TimeDists");

	TH3F *fRad3D = nullptr;
	TH3F *fPhi3D = nullptr;

	timeTree->SetBranchAddress("hIntDistortionR_negz",&fRad3D);
	timeTree->SetBranchAddress("hIntDistortionPhi_negz",&fPhi3D);

	timeTree->GetEntry(0);

	TH2F *fluctR = new TH2F("fluctR",";#phi;R [cm];Radial Distortion [cm]",fRad3D->GetNbinsX(),fRad3D->GetXaxis()->GetBinLowEdge(1),fRad3D->GetXaxis()->GetBinLowEdge(fRad3D->GetNbinsX()+1),fRad3D->GetNbinsY(),fRad3D->GetYaxis()->GetBinLowEdge(1),fRad3D->GetYaxis()->GetBinLowEdge(fRad3D->GetNbinsY()+1));
	TH2F *fluctPhi = new TH2F("fluctPhi",";#phi;R [cm];#phi Distortion [rad]",fPhi3D->GetNbinsX(),fPhi3D->GetXaxis()->GetBinLowEdge(1),fPhi3D->GetXaxis()->GetBinLowEdge(fPhi3D->GetNbinsX()+1),fPhi3D->GetNbinsY(),fPhi3D->GetYaxis()->GetBinLowEdge(1),fPhi3D->GetYaxis()->GetBinLowEdge(fPhi3D->GetNbinsY()+1));

	int zBin = fRad3D->GetZaxis()->FindBin(-1.0);

	for(int i=1; i<=fRad3D->GetNbinsX(); i++){
		for(int j=1; j<=fRad3D->GetNbinsY(); j++){
			fluctR->SetBinContent(i,j, fRad3D->GetBinContent(i,j,zBin) - sRad3D->GetBinContent(i,j,zBin) );
			fluctPhi->SetBinContent(i,j, fPhi3D->GetBinContent(i,j,zBin) - sPhi3D->GetBinContent(i,j,zBin) );
		}
	}


	TCanvas *c1 = new TCanvas("c1","",1000,1000);
	c1->SetRightMargin(0.175);
	fluctR->GetZaxis()->SetMaxDigits(2);
	fluctR->Draw("COLZ");
	c1->SaveAs("radialDistortion_input.svg");

	c1->Clear();
	c1->SetRightMargin(0.185);
	fluctPhi->GetZaxis()->SetTitleOffset(1.4);
	fluctPhi->GetZaxis()->SetMaxDigits(1);
	fluctPhi->Draw("COLZ");
	c1->SaveAs("phiDistortion_input.svg");




}
