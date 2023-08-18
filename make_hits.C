#include "TH2_Interpolate.C"
#include "TH2_Interpolate_hist.C"

const int nRadii = 8;
const int nStripes_R1 = 6;
const int nStripes_R2 = 8;
const int nStripes_R3 = 12;

int nPads_R1 = 6 * 16;
int nPads_R2 = 8 * 16;
int nPads_R3 = 12 * 16;

double pi = TMath::Pi();

double mm = 1.0;
double cm = 10.0;

std::array<double, nRadii> R1_e = {{227.0902789 * mm, 238.4100043 * mm, 249.7297296 * mm, 261.049455 * mm, 272.3691804 * mm, 283.6889058 * mm, 295.0086312 * mm, 306.3283566 * mm}};
std::array<double, nRadii> R1 = {{317.648082 * mm, 328.9678074 * mm, 340.2875328 * mm, 351.6072582 * mm, 362.9269836 * mm, 374.246709 * mm, 385.5664344 * mm, 396.8861597 * mm}};
std::array<double, nRadii> R2 = {{421.705532 * mm, 442.119258 * mm, 462.532984 * mm, 482.9467608 * mm, 503.36069 * mm, 523.774416 * mm, 544.188015 * mm, 564.601868 * mm}};
std::array<double, nRadii> R3 = {{594.6048725 * mm, 616.545823 * mm, 638.4867738 * mm, 660.4277246 * mm, 682.3686754 * mm, 704.3096262 * mm, 726.250577 * mm, 748.1915277 * mm}};

std::array<int, nRadii> keepThisAndAfter = {{1, 0, 1, 0, 1, 0, 1, 0}};
std::array<int, nRadii> keepUntil_R1_e = {{4, 4, 5, 4, 5, 5, 5, 5}};
std::array<int, nRadii> keepUntil_R1 = {{5, 5, 6, 5, 6, 5, 6, 5}};
std::array<int, nRadii> keepUntil_R2 = {{7, 7, 8, 7, 8, 8, 8, 8}};
std::array<int, nRadii> keepUntil_R3 = {{11, 10, 11, 11, 11, 11, 12, 11}};

std::array<int, nRadii> nGoodStripes_R1_e = {};
std::array<int, nRadii> nGoodStripes_R1 = {};
std::array<int, nRadii> nGoodStripes_R2 = {};
std::array<int, nRadii> nGoodStripes_R3 = {};

double x_R1e[nStripes_R1][nRadii];
double y_R1e[nStripes_R1][nRadii];

double x_R1[nStripes_R1][nRadii];
double y_R1[nStripes_R1][nRadii];

double x_R2[nStripes_R2][nRadii];
double y_R2[nStripes_R2][nRadii];

double x_R3[nStripes_R3][nRadii];
double y_R3[nStripes_R3][nRadii];

vector<pair<TVector3,int>> hits2Index;
vector<pair<TVector3,unsigned int>> hits2;
vector<TVector3> hits;
vector<int> hitIndex;
vector<double> phi_gaps;
vector<double> r_gaps;

vector<double> phi_boundaries;
vector<double> r_boundaries;


double avgAngle = 0.0;

TH2D *r_phi_hist;
TH1D *phi_R1;
TH1D *phi_R2;
TH1D *phi_R3;


vector<double> hist_phi_gaps_R1;
vector<double> hist_phi_gaps_R2;
vector<double> hist_phi_gaps_R3;

double avgAngle_R1 = 0.0;
double avgAngle_R2 = 0.0;
double avgAngle_R3 = 0.0;


double padfrac_R1 = 0.5 * 5.59106385 * mm;
double padfrac_R2 = 0.5 * 10.13836283 * mm;
double padfrac_R3 = 0.5 * 10.90189537 * mm;

double str_width_R1_e[nStripes_R1][nRadii] ={};
double str_width_R1[nStripes_R1][nRadii] ={};
double str_width_R2[nStripes_R2][nRadii] ={};
double str_width_R3[nStripes_R3][nRadii] ={};

std::array<double, nRadii> widthmod_R1_e = {{1.493, 1.398, 1.334, 1.284, 1.243, 1.208, 1.178, 1.152}};
std::array<double, nRadii> widthmod_R1 = {{1.129, 1.109, 1.091, 1.076, 1.062, 1.050, 1.040, 1.030}};
std::array<double, nRadii> widthmod_R2 = {{1.015, 1.007, 1.002, 1.000, 1.001, 1.006, 1.013, 1.023}};
std::array<double, nRadii> widthmod_R3 = {{1.044, 1.064, 1.087, 1.115, 1.147, 1.186, 1.232, 1.288}};

bool compare_for_sort(const TVector3& v1, const TVector3& v2){
  if( v1.Perp() == v2.Perp() || fabs(v2.Perp() - v1.Perp()) < 0.00001){
    return v1.Phi() < v2.Phi();
  }else return v1.Perp() < v2.Perp();
}

bool compare_for_sort2(const pair<TVector3,unsigned int>& v1, const pair<TVector3,unsigned int>& v2){
  if( v1.first.Perp() == v2.first.Perp() || fabs(v2.first.Perp() - v1.first.Perp()) < 0.00001){
    return v1.first.Phi() < v2.first.Phi();
  }else return v1.first.Perp() < v2.first.Perp();
}

void calcVertices(int nStripes, int nPads, std::array<double, nRadii>& R, std::array<int, nRadii>& keepUntil, double x[][nRadii], double y[][nRadii], double padfrac, double str_width[][nRadii], const std::array<double, nRadii>& widthmod, std::array<int, nRadii>& nGoodStripes){

  double phi_mod = pi/6.0;
  int pr_mult = 3;
  int dw_mult = 8;
  double diffwidth = 0.6;
  double adjust = 0.015;

  double theta = 0.0;
  
  //double x[nStripes][nRadii];
  //double y[nStripes][nRadii];

  std::array<double, nRadii> spacing;

  for(int i=0; i<nRadii; i++){
    spacing[i] = 2.0 * ((dw_mult * diffwidth / R[i]) + (pr_mult * phi_mod /nPads) );
  }//end loop over radii

  for(int j=0; j<nRadii; j++){
    int i_out = 0;
    for(int i= keepThisAndAfter[j]; i<keepUntil[j]; i++){
      if( j%2 == 0){
       theta = i * spacing[j] + (spacing[j]/2) - adjust;
       x[i_out][j] = R[j] * cos(theta);
       y[i_out][j] = R[j] * sin(theta);
     }else{
       theta = (i+1) * spacing[j] - adjust;
       x[i_out][j] = R[j] * cos(theta);
       y[i_out][j] = R[j] * sin(theta);
     }


     cout << "widthmod[" << j << "]: " << widthmod[j] << endl;
     cout << "str_width[" << i << "][" << j << "]: " << str_width[i][j] << endl;

     str_width[i][j] = 1.0 * mm;

     TVector3 corner[4];
     corner[0].SetXYZ(-padfrac + 0.5*mm,-(widthmod[j] * str_width[i][j]) / 2, 0);
     corner[1].SetXYZ(padfrac - 0.5*mm,-(widthmod[j] * str_width[i][j]) / 2, 0);
     corner[2].SetXYZ(-padfrac + 0.5*mm,(widthmod[j] * str_width[i][j]) / 2, 0);
     corner[3].SetXYZ(padfrac - 0.5*mm,(widthmod[j] * str_width[i][j]) / 2, 0);


     for(int k=0; k<4; k++){
       cout << "corner[" << k << "] X=" << corner[k].X() << "  Y=" << corner[k].Y() << endl;
     }

     cout << "width=" << corner[3].X() - corner[0].X() << "   height=" << corner[3].Y() - corner[0].Y() << endl;


      //cout << "i: " << i << "   i_out: " << i_out << endl;

     i_out++;


    }//end loop over phi
    nGoodStripes[j] = i_out;
  }//end loop over radii

  /*
  for(int i=0; i<i_out; i++){
    vector<double> tmpX;
    vector<double> tmpY;
    for(int j=0; j<nRadii; j++){
      tmpX.push_back(x[i][j]);
      tmpY.push_back(y[i][j]);
    }
    cx.push_back(tmpX);
    cy.push_back(tmpY);
  }
  */
}


void getAllVert(){

  calcVertices(nStripes_R1,nPads_R1,R1_e,keepUntil_R1_e,x_R1e,y_R1e,padfrac_R1,str_width_R1_e,widthmod_R1_e,nGoodStripes_R1_e);
  for(int i=0; i<nRadii; i++) cout << "nGoodStripes R1e " << i << " R=" << R1_e[i] << ": " << nGoodStripes_R1_e[i] << endl;
    calcVertices(nStripes_R1,nPads_R1,R1,keepUntil_R1,x_R1,y_R1,padfrac_R1,str_width_R1,widthmod_R1,nGoodStripes_R1);
  for(int i=0; i<nRadii; i++) cout << "nGoodStripes R1 " << i << " R=" << R1[i] << ": " << nGoodStripes_R1[i] << endl;
    calcVertices(nStripes_R2,nPads_R2,R2,keepUntil_R2,x_R2,y_R2,padfrac_R2,str_width_R2,widthmod_R2,nGoodStripes_R2);
  for(int i=0; i<nRadii; i++) cout << "nGoodStripes R2 " << i << " R=" << R2[i] << ": " << nGoodStripes_R2[i] << endl;
    calcVertices(nStripes_R3,nPads_R3,R3,keepUntil_R3,x_R3,y_R3,padfrac_R3,str_width_R3,widthmod_R3,nGoodStripes_R3);
  for(int i=0; i<nRadii; i++) cout << "nGoodStripes R3 " << i << " R=" << R3[i] << ": " << nGoodStripes_R3[i] << endl;



}

TVector3 getHitPos(int petal, int module, int radius, int stripe){
  TVector3 pos;
  if(module == 0) pos.SetXYZ(x_R1e[stripe][radius] / cm,y_R1e[stripe][radius] / cm,0);
  else if(module == 1) pos.SetXYZ(x_R1[stripe][radius] / cm,y_R1[stripe][radius] / cm,0);
  else if(module == 2) pos.SetXYZ(x_R2[stripe][radius] / cm,y_R2[stripe][radius] / cm,0);
  else if(module == 3) pos.SetXYZ(x_R3[stripe][radius] / cm,y_R3[stripe][radius] / cm,0);
  
  if(petal > 0) pos.RotateZ(petal * (pi/9.0));
  
  return pos;

}

void getAllHits(){

  getAllVert();

  TFile *treeFile = new TFile("inputFiles/truthSectors.root","RECREATE");
  TTree *truthTree = new TTree("truthTree","truthTree");
  unsigned int truthID;
  TVector3 truthPos;
  int sector;
  
  truthTree->Branch("truthID",&truthID);
  truthTree->Branch("truthPos",&truthPos);
  truthTree->Branch("sector",&sector);

  fstream out("CM_HitTruth.txt",ios::out);
  
  for(int i=0; i<18; i++){
    for(int j=0; j<8; j++){
      for(int k=0; k<nGoodStripes_R1_e[j]; k++){
       hits.push_back( getHitPos(i,0,j,k));
       out << i*10000 + j*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << 1.0 << endl;
       out << (i+18)*10000 + j*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << -1.0 << endl;
       hitIndex.push_back(i*10000 + j*100 + k);
     }

     for(int k=0; k<nGoodStripes_R1[j]; k++){
       hits.push_back( getHitPos(i,1,j,k));
       out << i*10000 + (j+8)*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << 1.0 << endl;
       out << (i+18)*10000 + (j+8)*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << -1.0 << endl;
       hitIndex.push_back(i*10000 + (j+8)*100 + k);
     }

     for(int k=0; k<nGoodStripes_R2[j]; k++){
       hits.push_back( getHitPos(i,2,j,k));
       out << i*10000 + (j+16)*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << 1.0 << endl;
       out << (i+18)*10000 + (j+16)*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << -1.0 << endl;
       hitIndex.push_back(i*10000 + (j+16)*100 + k);
     }

     for(int k=0; k<nGoodStripes_R3[j]; k++){
       hits.push_back( getHitPos(i,3,j,k));
       out << i*10000 + (j+24)*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << 1.0 << endl;
       out << (i+18)*10000 + (j+24)*100 + k << "\t" << hits[hits.size()-1].X() << "\t" << hits[hits.size()-1].Y() << "\t" << -1.0 << endl;
       hitIndex.push_back(i*10000 + (j+24)*100 + k);
     }

   }
 }

 for(int i=0; i<hits.size(); i++){
  pair<TVector3,unsigned int> tmp_pair;
  tmp_pair.first = hits[i];
  tmp_pair.second = i;
  hits2.push_back(tmp_pair);

  pair<TVector3,int> tmp_index;
  tmp_index.first = hits[i];
  tmp_index.first.SetZ(1.0);
  tmp_index.second = hitIndex[i];
  hits2Index.push_back(tmp_index);

  truthID = tmp_index.second;
  truthPos = tmp_index.first;
  double phi = tmp_index.first.Phi();
  if(phi < 0.0) phi += 2.0*TMath::Pi();

  if(phi < 1.0*TMath::Pi()/12 || phi >= 23*TMath::Pi()/12){
    sector = 0;
    cout << "sector low: 23*pi/12   sector high: pi/12" << endl;
  }else{
    for(int p=0; p<12; p++){
     if(phi >= ((2*p)-1)*TMath::Pi()/12 && phi < ((2*p)+1)*TMath::Pi()/12){

      cout << "sector low: " << ((2*p)-1) << "*pi/12   sector high: " << ((2*p)+1) << "*pi/12" << endl;
      sector = p;
      continue;
    }
  }
}



cout << "phi=" << phi << "   sector=" << sector << endl;
truthTree->Fill();

tmp_index.first.SetZ(-1.0);
tmp_index.second = hitIndex[i] + (18*10000);
hits2Index.push_back(tmp_index);

truthID = tmp_index.second;
truthPos = tmp_index.first;
sector += 12;
truthTree->Fill();

}

truthTree->Write();
treeFile->Close();

}

void getPhiGaps(vector<TVector3>& reco_pos){
  double prev_phi = -TMath::Pi();
  for(int i=0; i<(int)reco_pos.size(); i++){
    if(	reco_pos[i].Phi() < prev_phi ) prev_phi	-= 2.0*TMath::Pi();

    if(	(reco_pos[i].Phi() - prev_phi) > 0.075) phi_gaps.push_back(0.5*(reco_pos[i].Phi() + prev_phi));
    prev_phi = reco_pos[i].Phi();

  }
  //  return phi_gaps;

}

void getRGaps(vector<TVector3>& reco_pos){

  double prev_r = 20;
  for(int i=0; i<(int)reco_pos.size(); i++){
    if( (reco_pos[i].Perp() < 50 && (reco_pos[i].Perp() - prev_r) > 2.2) || (reco_pos[i].Perp() > 50 &&	(reco_pos[i].Perp() - prev_r) >	2.5)) r_gaps.push_back(	0.5 * (reco_pos[i].Perp() + prev_r));
    prev_r = reco_pos[i].Perp();
  }
  //  return r_gaps;
}

void getGaps(){
  getAllHits();

  //  std::sort(hits.begin(),hits.end(),compare_for_sort);
  std::sort(hits2.begin(),hits2.end(),compare_for_sort2);

  for(int i=0; i<hits2.size(); i++){
    hits[i] = hits2[i].first;
  }
  
  getPhiGaps(hits);
  getRGaps(hits);
  
  double phi_tmp[18] = {0.0};
  int nPhi[18] = {0};
  double r_tmp[3] = {0.0};
  int nR[3] = {0};

  for(int i=0; i<phi_gaps.size(); i++){
    for(int p=0; p<18; p++){
      if(nPhi[p] == 0 || fabs( (phi_tmp[p]/nPhi[p]) - phi_gaps[i] ) < 0.2 ){
       phi_tmp[p] += phi_gaps[i];
       nPhi[p]++;
       break;
     }
   }
 }

 for(int i=0; i<r_gaps.size(); i++){
  for(int p=0; p<3; p++){
    if(nR[p] == 0 || fabs( (r_tmp[p]/nR[p]) - r_gaps[i] ) < 10 ){
     r_tmp[p] += r_gaps[i];
     nR[p]++;
     break;
   }
 }
}

for(int i=0; i<18; i++){
  phi_boundaries.push_back(phi_tmp[i] / nPhi[i]);
}

for(int i=0; i<3; i++){
  r_boundaries.push_back(r_tmp[i] / nR[i]);
}

}

void getAvgAngle(){

  double angles[18];
  double diff[18];
  double sum = 0.0;
  
  for(int i=0; i<18; i++){
    angles[i] = -pi + (i*pi/9.0);
    diff[i] = phi_boundaries[i] - angles[i];
    sum += diff[i];
  }

  avgAngle = sum/18.0;
  
}

void makeHitHist(){

  r_phi_hist = new TH2D("r_phi_hist","",360,-TMath::Pi(),TMath::Pi(),500,0.0,100.0);

  for(int i=0; i<hits.size(); i++){
    r_phi_hist->Fill(hits[i].Phi(), hits[i].Perp());
  }

  phi_R1 = r_phi_hist->ProjectionX("",151,206);
  phi_R2 = r_phi_hist->ProjectionX("",206,290);
  phi_R3 = r_phi_hist->ProjectionX("",290,500);

  return 0;
}

void getHistAngles(){

  double angles[18];
  double diff[18];
  double sum = 0.0;
  
  for(int i=0; i<18; i++){
    angles[i] = -pi + (i*pi/9.0);
  }
  
  for(int i=2; i<=phi_R1->GetNbinsX()+1; i++){
    if(phi_R1->GetBinContent(i) > 0 && phi_R1->GetBinContent(i-1) == 0.0){
      if(hist_phi_gaps_R1.size() == 0) hist_phi_gaps_R1.push_back(phi_R1->GetBinCenter(i));
      else if(phi_R1->GetBinCenter(i) - hist_phi_gaps_R1[hist_phi_gaps_R1.size()-1] > (TMath::Pi()/36.)) hist_phi_gaps_R1.push_back(phi_R1->GetBinCenter(i));
    }
  }


  for(int i=2; i<=phi_R2->GetNbinsX()+1; i++){
    if(phi_R2->GetBinContent(i) > 0 && phi_R2->GetBinContent(i-1) == 0.0){
      if(hist_phi_gaps_R2.size() == 0) hist_phi_gaps_R2.push_back(phi_R2->GetBinCenter(i));
      else if(phi_R2->GetBinCenter(i) - hist_phi_gaps_R2[hist_phi_gaps_R2.size()-1] > (TMath::Pi()/36.)) hist_phi_gaps_R2.push_back(phi_R2->GetBinCenter(i));
    }
  }
  
  
  for(int i=2; i<=phi_R3->GetNbinsX()+1; i++){
    if(phi_R3->GetBinContent(i) > 0 && phi_R3->GetBinContent(i-1) == 0.0){
      if(hist_phi_gaps_R3.size() == 0) hist_phi_gaps_R3.push_back(phi_R3->GetBinCenter(i));
      else if(phi_R3->GetBinCenter(i) - hist_phi_gaps_R3[hist_phi_gaps_R3.size()-1] > (TMath::Pi()/36.)) hist_phi_gaps_R3.push_back(phi_R3->GetBinCenter(i));
    }
  }

  for(int i=0; i<18; i++){
    diff[i] = hist_phi_gaps_R1[i] - angles[i];
    sum += diff[i];
  }
  avgAngle_R1 = sum/18.;

  sum = 0.0;
  for(int i=0; i<18; i++){
    diff[i] = hist_phi_gaps_R2[i] - angles[i];
    sum += diff[i];
  }
  avgAngle_R2 = sum/18.;

  sum = 0.0;
  for(int i=0; i<18; i++){
    diff[i] = hist_phi_gaps_R3[i] - angles[i];
    sum += diff[i];
  }
  avgAngle_R3 = sum/18.;

  
  return 0;
}

vector<double> getPhiAngles_smooth( TH2F *r_phi ){

  vector<double> avgRot;
  
  TH1D *proj[3];

  proj[0] = r_phi->ProjectionX("proj_R1",151,206);
  proj[1] = r_phi->ProjectionX("proj_R2",206,290);
  proj[2] = r_phi->ProjectionX("proj_R3",290,499);

  for(int R=0; R<3; R++){

    vector<double> phiVec;
    
    for(int i=1; i<=proj[R]->GetNbinsX(); i++){
      if(proj[R]->GetBinContent(i) == 0) continue;
      int n = 0;
      while(n < proj[R]->GetBinContent(i)){
       phiVec.push_back(proj[R]->GetBinCenter(i));
       n++;
     }
   }

   TKDE *kernel = new TKDE(proj[R]->GetEntries(), &phiVec[0], -TMath::Pi(), TMath::Pi(), "", 0.1);

   TH1 *smoothedHist = kernel->GetFunction()->GetHistogram();

   TF1 *gausFit = new TF1("gausFit","[2]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[3]-[0])/[1],2)) + [4]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[5]-[0])/[1],2)) + [6]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[7]-[0])/[1],2)) + [8]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[9]-[0])/[1],2)) + [10]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[11]-[0])/[1],2)) + [12]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[13]-[0])/[1],2)) + [14]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[15]-[0])/[1],2)) + [16]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[17]-[0])/[1],2)) + [18]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[19]-[0])/[1],2)) + [20]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[21]-[0])/[1],2)) + [22]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[23]-[0])/[1],2)) + [24]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[25]-[0])/[1],2)) + [26]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[27]-[0])/[1],2)) + [28]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[29]-[0])/[1],2)) + [30]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[31]-[0])/[1],2)) + [32]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[33]-[0])/[1],2)) + [34]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[35]-[0])/[1],2)) + [36]*[1]*sqrt(2.0*TMath::Pi())*exp(-0.5*pow((x-[37]-[0])/[1],2))", -TMath::Pi(), TMath::Pi());

   gausFit->SetNpx(1000);

   for(int i=0; i<18; i++){
    double phiValue = -TMath::Pi() + ((2*i + 1)*TMath::Pi()/18);
    gausFit->SetParameter(2*i+2,smoothedHist->GetBinContent(smoothedHist->FindBin(phiValue)));
    gausFit->FixParameter(2*i+3, phiValue);
  }
  gausFit->SetParameter(0,0.0);
  gausFit->SetParameter(1,TMath::Pi()/18);

  smoothedHist->Fit("gausFit");

  avgRot.push_back(gausFit->GetParameter(0));

  }//end loop over modules


  return avgRot;
  
}

vector<double> getRGapsHist( TH2F *r_phi ){

  TH1D *proj = r_phi->ProjectionY("R_proj",1,360);
  
  std::vector<double> r_rows;
  
  for(int i=2; i<proj->GetNbinsX(); i++){
    if(proj->GetBinContent(i) > 0.15*proj->GetMaximum() && proj->GetBinContent(i) >= proj->GetBinContent(i-1) && proj->GetBinContent(i) >= proj->GetBinContent(i+1)) r_rows.push_back(proj->GetBinCenter(i));
  }
  
  for(int i=0; i<(int)r_rows.size()-1; i++){
    if(r_rows[i+1]-r_rows[i] > 0.75) continue;
    
    if(proj->GetBinContent(proj->FindBin(r_rows[i])) > proj->GetBinContent(proj->FindBin(r_rows[i+1]))) r_rows.erase(std::next(r_rows.begin(), i+1));
    else r_rows.erase(std::next(r_rows.begin(), i));
    
    i--;
  }
  
  return r_rows;
  
}

vector<double> getRGapsHist_gaps(vector<double> rad){

  vector<double> gaps;
  
  for(int i=0; i<rad.size()-1; i++){
    gaps.push_back(rad[i+1] - rad[i]);
  }

  return gaps;
  
}


vector<int> getRMatching(vector<double> hits, vector<double> clusts, vector<double> clustGaps){

  vector<int> RMatches;

  int R23Gap = -1;
  
  for(int c=0; c<clustGaps.size(); c++){
    if(clustGaps[c] > 2.5){
      R23Gap = c;
      break;
    }
  }

  for(int i=0; i<clusts.size(); i++){
    RMatches.push_back(i + 23 - R23Gap);
  }

  return RMatches;
  
}

int getClusterRMatch(vector<int> rMatches, vector<double> clustRs, double cr){

  int closest_clustR = -1;
  
  for(int i=0; i<rMatches.size(); i++){

    double lowGap = 0.0;
    double highGap = 0.0;

    if(rMatches[i] <= 14){
      lowGap = 0.565985;
      highGap = 0.565985;
    }else if(rMatches[i] == 15){
      lowGap = 0.565985;
      highGap = 1.2409686;
    }else if(rMatches[i] == 16){
      lowGap = 1.2409686;
      highGap = 1.020695;
    }else if(rMatches[i] >= 17 && rMatches[i] <= 22){
      lowGap = 1.020695;
      highGap = 1.020695;
    }else if(rMatches[i] == 23){
      lowGap = 1.020695;
      highGap = 1.5001502;
    }else if(rMatches[i] == 24){
      lowGap = 1.5001502;
      highGap = 1.09705;
    }else if(rMatches[i] >= 25){
      lowGap = 1.09705;
      highGap = 1.09705;
    }

    if( cr > (clustRs[i] - lowGap) && cr <= (clustRs[i] + highGap) ){
      closest_clustR = rMatches[i];
      break;
    }

  }

  return closest_clustR;
  
}

TH2F* doHistSubMap(TH2F *reco, TH3F *in){

  //  TH2F *diff = (TH2F*)reco->Clone();
  //  TH2F *diff = new TH2F("diff",";#phi;r (cm)",38,0.0,2.0*TMath::Pi(),18,20.0,78.0);
  TH2F *diff = new TH2F("diff",";#phi;r (cm)",800,-100,100,800,-80,80);

  //int zBin = -1;

  //  if(in->GetZaxis()->FindBin(1.0) == 0) zBin = in->GetZaxis()->FindBin(-1.0);
  //  else zBin = in->GetZaxis()->FindBin(1.0);

  double z = 1.0;
  if(in->GetZaxis()->FindBin(z) == 0) z = -1.0;


  //  TH2F *in2 = new TH2F("in2","",in->GetNbinsX(),in->GetXaxis()->GetBinLowEdge(1),in->GetXaxis()->GetBinLowEdge(in->GetNbinsX()),in->GetNbinsY(),in->GetYaxis()->GetBinLowEdge(1),in->GetYaxis()->GetBinLowEdge(in->GetNbinsY()));

  
  for(int i=0; i<hits.size(); i++){

    double phi = hits[i].Phi();

    if(phi < 0) phi += 2.0*TMath::Pi();

    double recoInt = TH2_Interpolate(reco,phi,hits[i].Perp());
    //double recoInt = TH2_Interpolate(in2,phi,hits[i].Perp());
    
    //    double recoInt = reco->Interpolate(phi,hits[i].Perp());
    double inInt = in->Interpolate(phi,hits[i].Perp(),z);

    cout << "Map phi: " << phi << "   R: " << hits[i].Perp() << "   reco interpolation: " << recoInt << "   in interpolation: " << inInt << endl;

    int binX = diff->GetXaxis()->FindBin(hits[i].X());
    int binY = diff->GetYaxis()->FindBin(hits[i].Y());

    if(recoInt == 0.0) continue;

  //    diff->SetBinContent(binX,binY, reco->Interpolate(hits[i].Phi(),hits[i].Perp()) - in->Interpolate(hits[i].Phi(),hits[i].Perp(),z));
    diff->SetBinContent(binX, binY, diff->GetBinContent(binX,binY) + recoInt - inInt);
  }

  return diff;
  
}

TH2F* doHistSub(TH2F *reco, TH3F *in){

  TH2F *diff = new TH2F("diff",";#phi;r (cm)",200,0.0,2.0*TMath::Pi(),58,20.0,78.0);
  
  int zBin = -1;

  if(in->GetZaxis()->FindBin(1.0) == 0) zBin = in->GetZaxis()->FindBin(-1.0);
  else zBin = in->GetZaxis()->FindBin(1.0);

  double z = 1.0;
  if(in->GetZaxis()->FindBin(z) == 0) z = -1.0;

  
  for(int i=1; i<=diff->GetNbinsX(); i++){
    for(int j=1; j<=diff->GetNbinsY(); j++){

      //      double recoInt = reco->Interpolate(diff->GetXaxis()->GetBinCenter(i),diff->GetYaxis()->GetBinCenter(j));
      double recoInt = TH2_Interpolate(reco,diff->GetXaxis()->GetBinCenter(i),diff->GetYaxis()->GetBinCenter(j));
      //double recoInt = reco->GetBinContent(i,j,zBin);
      double inInt = in->Interpolate(diff->GetXaxis()->GetBinCenter(i),diff->GetYaxis()->GetBinCenter(j),z);

      cout << "Phi: " << diff->GetXaxis()->GetBinCenter(i) << "   R: " << diff->GetYaxis()->GetBinCenter(j) << "   reco interpolation: " << recoInt << "   in interpolation: " << inInt << endl;

      if(recoInt == 0.0) continue;
      
      diff->SetBinContent(i,j, diff->GetBinContent(i,j) + recoInt - inInt);
    }
  }
  
  return diff;
  
}

TH2F* doHistRatMap(TH2F *reco, TH3F *in){

  TH2F *rat = new TH2F("rat",";#phi;r (cm)",800,-100,100,800,-80,80);

  double z = 1.0;
  if(in->GetZaxis()->FindBin(z) == 0) z = -1.0;


  for(int i=0; i<hits.size(); i++){

    double phi = hits[i].Phi();

    if(phi < 0) phi += 2.0*TMath::Pi();
    
    //    double recoInt = reco->Interpolate(phi,hits[i].Perp());
    double recoInt = TH2_Interpolate(reco,phi,hits[i].Perp());
    double inInt = in->Interpolate(phi,hits[i].Perp(),z);

    cout << "Map phi: " << phi << "   R: " << hits[i].Perp() << "   reco interpolation: " << recoInt << "   in interpolation: " << inInt << endl;

    int binX = rat->GetXaxis()->FindBin(hits[i].X());
    int binY = rat->GetYaxis()->FindBin(hits[i].Y());

    if(rat->GetBinContent(binX,binY) == 0.0) rat->SetBinContent(binX,binY, recoInt/inInt);
    else rat->SetBinContent(binX, binY, 0.5*(rat->GetBinContent(binX,binY) + (recoInt / inInt)) );
  }
  
  return rat;
  
}

TH2F* doHistRat(TH2F *reco, TH3F *in){

  TH2F *rat = new TH2F("rat",";#phi;r (cm)",200,0.0,2.0*TMath::Pi(),58,20.0,78.0);
  
  int zBin = -1;

  if(in->GetZaxis()->FindBin(1.0) == 0) zBin = in->GetZaxis()->FindBin(-1.0);
  else zBin = in->GetZaxis()->FindBin(1.0);

  double z = 1.0;
  if(in->GetZaxis()->FindBin(z) == 0) z = -1.0;

  
  for(int i=1; i<=rat->GetNbinsX(); i++){
    for(int j=1; j<=rat->GetNbinsY(); j++){

      //      double recoInt = reco->Interpolate(rat->GetXaxis()->GetBinCenter(i),rat->GetYaxis()->GetBinCenter(j));
      double recoInt = TH2_Interpolate(reco,rat->GetXaxis()->GetBinCenter(i),rat->GetYaxis()->GetBinCenter(j));
      double inInt = in->Interpolate(rat->GetXaxis()->GetBinCenter(i),rat->GetYaxis()->GetBinCenter(j),z);

      cout << "Phi: " << rat->GetXaxis()->GetBinCenter(i) << "   R: " << rat->GetYaxis()->GetBinCenter(j) << "   reco interpolation: " << recoInt << "   in interpolation: " << inInt << endl;
      
      rat->SetBinContent(i,j, recoInt/inInt);
    }
  }
  
  return rat;
  
}


void makeDistortionDiffMaps(TFile *inFile, TFile *recoFile){

  getAllHits();
  
  TH3F *inR_pos = (TH3F*)inFile->Get("hIntDistortionR_posz");
  TH3F *inP_pos = (TH3F*)inFile->Get("hIntDistortionP_posz");
  TH3F *inZ_pos = (TH3F*)inFile->Get("hIntDistortionZ_posz");
  TH3F *inR_neg = (TH3F*)inFile->Get("hIntDistortionR_negz");
  TH3F *inP_neg = (TH3F*)inFile->Get("hIntDistortionP_negz");
  TH3F *inZ_neg = (TH3F*)inFile->Get("hIntDistortionZ_negz");


  /*
  TH2F *recoR_pos = (TH2F*)recoFile->Get("hIntDistortionR_posz");
  TH2F *recoP_pos = (TH2F*)recoFile->Get("hIntDistortionP_posz");
  TH2F *recoZ_pos = (TH2F*)recoFile->Get("hIntDistortionZ_posz");
  TH2F *recoR_neg = (TH2F*)recoFile->Get("hIntDistortionR_negz");
  TH2F *recoP_neg = (TH2F*)recoFile->Get("hIntDistortionP_negz");
  TH2F *recoZ_neg = (TH2F*)recoFile->Get("hIntDistortionZ_negz");
  */

  
  TH2F *recoR_pos = TH2_Interpolate_hist((TH2F*)recoFile->Get("hIntDistortionR_posz"));
  TH2F *recoP_pos = TH2_Interpolate_hist((TH2F*)recoFile->Get("hIntDistortionP_posz"));
  TH2F *recoZ_pos = TH2_Interpolate_hist((TH2F*)recoFile->Get("hIntDistortionZ_posz"));
  TH2F *recoR_neg = TH2_Interpolate_hist((TH2F*)recoFile->Get("hIntDistortionR_negz"));
  TH2F *recoP_neg = TH2_Interpolate_hist((TH2F*)recoFile->Get("hIntDistortionP_negz"));
  TH2F *recoZ_neg = TH2_Interpolate_hist((TH2F*)recoFile->Get("hIntDistortionZ_negz"));
  
  TFile *diffFile = new TFile("distortionDiff.root","RECREATE");
  

  TH2F *diffMapR_pos = doHistSubMap(recoR_pos,inR_pos);
  diffMapR_pos->SetName("diffMapR_pos");
  diffMapR_pos->SetTitle("Reco-Input Distortion R z>0");

  TH2F *diffMapP_pos = doHistSubMap(recoP_pos,inP_pos);
  diffMapP_pos->SetName("diffMapP_pos");
  diffMapP_pos->SetTitle("Reco-Input Distortion P z>0");

  TH2F *diffMapZ_pos = doHistSubMap(recoZ_pos,inZ_pos);
  diffMapZ_pos->SetName("diffMapZ_pos");
  diffMapZ_pos->SetTitle("Reco-Input Distortion Z z>0");

  TH2F *diffMapR_neg = doHistSubMap(recoR_neg,inR_neg);
  diffMapR_neg->SetName("diffMapR_neg");
  diffMapR_neg->SetTitle("Reco-Input Distortion R z<0");

  TH2F *diffMapP_neg = doHistSubMap(recoP_neg,inP_neg);
  diffMapP_neg->SetName("diffMapP_neg");
  diffMapP_neg->SetTitle("Reco-Input Distortion P z<0");

  TH2F *diffMapZ_neg = doHistSubMap(recoZ_neg,inZ_neg);
  diffMapZ_neg->SetName("diffMapZ_neg");
  diffMapZ_neg->SetTitle("Reco-Input Distortion Z z<0");


  TH2F *diffR_pos = doHistSub(recoR_pos,inR_pos);
  diffR_pos->SetName("diffR_pos");
  diffR_pos->SetTitle("Reco-Input Distortion R z>0");

  TH2F *diffP_pos = doHistSub(recoP_pos,inP_pos);
  diffP_pos->SetName("diffP_pos");
  diffP_pos->SetTitle("Reco-Input Distortion P z>0");

  TH2F *diffZ_pos = doHistSub(recoZ_pos,inZ_pos);
  diffZ_pos->SetName("diffZ_pos");
  diffZ_pos->SetTitle("Reco-Input Distortion Z z>0");

  TH2F *diffR_neg = doHistSub(recoR_neg,inR_neg);
  diffR_neg->SetName("diffR_neg");
  diffR_neg->SetTitle("Reco-Input Distortion R z<0");

  TH2F *diffP_neg = doHistSub(recoP_neg,inP_neg);
  diffP_neg->SetName("diffP_neg");
  diffP_neg->SetTitle("Reco-Input Distortion P z<0");

  TH2F *diffZ_neg = doHistSub(recoZ_neg,inZ_neg);
  diffZ_neg->SetName("diffZ_neg");
  diffZ_neg->SetTitle("Reco-Input Distortion Z z<0");








  TH2F *ratMapR_pos = doHistRatMap(recoR_pos,inR_pos);
  ratMapR_pos->SetName("ratMapR_pos");
  ratMapR_pos->SetTitle("Reco/Input Distortion R z>0");

  TH2F *ratMapP_pos = doHistRatMap(recoP_pos,inP_pos);
  ratMapP_pos->SetName("ratMapP_pos");
  ratMapP_pos->SetTitle("Reco/Input Distortion P z>0");

  TH2F *ratMapZ_pos = doHistRatMap(recoZ_pos,inZ_pos);
  ratMapZ_pos->SetName("ratMapZ_pos");
  ratMapZ_pos->SetTitle("Reco/Input Distortion Z z>0");

  TH2F *ratMapR_neg = doHistRatMap(recoR_neg,inR_neg);
  ratMapR_neg->SetName("ratMapR_neg");
  ratMapR_neg->SetTitle("Reco/Input Distortion R z<0");

  TH2F *ratMapP_neg = doHistRatMap(recoP_neg,inP_neg);
  ratMapP_neg->SetName("ratMapP_neg");
  ratMapP_neg->SetTitle("Reco/Input Distortion P z<0");

  TH2F *ratMapZ_neg = doHistRatMap(recoZ_neg,inZ_neg);
  ratMapZ_neg->SetName("ratMapZ_neg");
  ratMapZ_neg->SetTitle("Reco/Input Distortion Z z<0");


  TH2F *ratR_pos = doHistRat(recoR_pos,inR_pos);
  ratR_pos->SetName("ratR_pos");
  ratR_pos->SetTitle("Reco/Input Distortion R z>0");

  TH2F *ratP_pos = doHistRat(recoP_pos,inP_pos);
  ratP_pos->SetName("ratP_pos");
  ratP_pos->SetTitle("Reco/Input Distortion P z>0");

  TH2F *ratZ_pos = doHistRat(recoZ_pos,inZ_pos);
  ratZ_pos->SetName("ratZ_pos");
  ratZ_pos->SetTitle("Reco/Input Distortion Z z>0");

  TH2F *ratR_neg = doHistRat(recoR_neg,inR_neg);
  ratR_neg->SetName("ratR_neg");
  ratR_neg->SetTitle("Reco/Input Distortion R z<0");

  TH2F *ratP_neg = doHistRat(recoP_neg,inP_neg);
  ratP_neg->SetName("ratP_neg");
  ratP_neg->SetTitle("Reco/Input Distortion P z<0");

  TH2F *ratZ_neg = doHistRat(recoZ_neg,inZ_neg);
  ratZ_neg->SetName("ratZ_neg");
  ratZ_neg->SetTitle("Reco/Input Distortion Z z<0");

  
  
  diffFile->cd();
  
  diffMapR_pos->Write();
  diffMapP_pos->Write();
  diffMapZ_pos->Write();

  diffMapR_neg->Write();
  diffMapP_neg->Write();
  diffMapZ_neg->Write();


  diffR_pos->Write();
  diffP_pos->Write();
  diffZ_pos->Write();

  diffR_neg->Write();
  diffP_neg->Write();
  diffZ_neg->Write();




  ratMapR_pos->Write();
  ratMapP_pos->Write();
  ratMapZ_pos->Write();

  ratMapR_neg->Write();
  ratMapP_neg->Write();
  ratMapZ_neg->Write();


  ratR_pos->Write();
  ratP_pos->Write();
  ratZ_pos->Write();

  ratR_neg->Write();
  ratP_neg->Write();
  ratZ_neg->Write();

  
  diffFile->Close();
  
}

void getHitsSelected(vector<double>& x, vector<double>& y, bool pos=true, int petal=-1, int row=-1, int stripe=-1){
  x.clear();
  y.clear();
  for(int i=0; i<hits2Index.size(); i++){
    int index = hits2Index[i].second;
    int index_petal = index/10000;
    if((pos && index_petal >= 18) || (!pos && index_petal < 18)) continue;
    if(petal != -1 && index_petal != petal) continue;
    int index_row = (index - (index_petal*10000))/100;
    if( row != -1 && index_row != row) continue;
    int index_stripe = index - (index_petal*10000) - (index_row*100);
    if(stripe != -1 && index_stripe != stripe) continue;

    x.push_back(hits2Index[i].first.X());
    y.push_back(hits2Index[i].first.Y());
  }
  return;
}

TH2D *h_empty;
TCanvas *canv = new TCanvas();
TGraph *gr;
void drawHitsSelected(bool pos=true, int petal=-1, int row=-1, int stripe=-1){
  vector<double> x;
  vector<double> y;

  if(!h_empty) h_empty = new TH2D("h_empty","",1,-90,90,1,-90,90);
  gStyle->SetOptStat(0);
  canv->Clear();
  h_empty->Draw();

  getHitsSelected(x,y,pos,petal,row,stripe);

  gr = new TGraph(x.size(),&x[0],&y[0]);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.25);
  gr->Draw("PSAME");

  return;
  
}
