#include "genrp_tree.C"
#include <TROOT.h>
#include <TRandom3.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TDirectory.h> 
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TChain.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

using namespace std;

// options
const Bool_t   ApplyElec  = true;
const Bool_t   ApplyElas  = true;
const Bool_t   PlotKine   = true;
const Bool_t   PlotBBCal  = true;
const Bool_t   PlotHCal   = true;
const Bool_t   PlotPol    = true;

//-----------------------------------------------------------------------------------------------------------------------------
// sclose and zclose calculation:
//-----------------------------------------------------------------------------------------------------------------------------

void calc_sclose_zclose( TVector3 Track1_Coord, TVector3 Track2_Coord, TVector3 Track1_Slope, TVector3 Track2_Slope, double &sclose, double &zclose ){

  double x1 = Track1_Coord.X() - Track1_Coord.Z()*Track1_Slope.X();
  double y1 = Track1_Coord.Y() - Track1_Coord.Z()*Track1_Slope.Y();
  double x2 = Track2_Coord.X() - Track2_Coord.Z()*Track2_Slope.X();
  double y2 = Track2_Coord.Y() - Track2_Coord.Z()*Track2_Slope.Y();

  double xp1 = Track1_Slope.X();
  double yp1 = Track1_Slope.Y();
  double xp2 = Track2_Slope.X();
  double yp2 = Track2_Slope.Y();

  TMatrixD Mclose(2,2);
  TVectorD bclose(2);

  Mclose(0,0) = 1.0 + pow(xp1,2) + pow(yp1,2);
  Mclose(0,1) = -(1.0 + xp1*xp2 + yp1*yp2);
  Mclose(1,0) = Mclose(0,1);
  Mclose(1,1) = 1.0 + pow(xp2,2) + pow(yp2,2);

  bclose(0) = xp1*(x2-x1) + yp1*(y2-y1);
  bclose(1) = xp2*(x1-x2) + yp2*(y1-y2);

  TVectorD zClose = Mclose.Invert() * bclose;

  double z1 = zClose(0);
  double z2 = zClose(1);

  double sClose2 = pow( x1 + xp1*z1 - (x2 + xp2*z2), 2 ) + pow( y1 + yp1*z1 - (y2 + yp2*z2), 2 ) + pow( z1-z2, 2 );

  sclose = sqrt(sClose2);
  zclose = 0.5*(z1 + z2 );
}

//-----------------------------------------------------------------------------------------------------------------------------
// conetest
//-----------------------------------------------------------------------------------------------------------------------------

bool conetest( TVector3 Track1_Coord, TVector3 Track1_Slope, double theta, double zclose, double zback, double Lx=2.0, double Ly=0.6, double xcenter=0.0, double ycenter=0.0 ){
  double xfp, yfp, xpfp, ypfp;

  xfp = Track1_Coord.X() - Track1_Coord.Z() * Track1_Slope.X();
  yfp = Track1_Coord.Y() - Track1_Coord.Z() * Track1_Slope.Y();
  xpfp = Track1_Slope.X();
  ypfp = Track1_Slope.Y();

  double xclose = xfp + xpfp*zclose;
  double yclose = yfp + ypfp*zclose;

  double xpplus = (xpfp + tan(theta))/(1.-xpfp*tan(theta));
  double xpminus = (xpfp - tan(theta))/(1.+xpfp*tan(theta));
  double ypplus = (ypfp + tan(theta))/(1.-ypfp*tan(theta));
  double ypminus = (ypfp - tan(theta))/(1.+ypfp*tan(theta));

  double xmax = xclose + xpplus * (zback - zclose);
  double xmin = xclose + xpminus * (zback - zclose);
  double ymax = yclose + ypplus * (zback - zclose);
  double ymin = yclose + ypminus * (zback - zclose);

  return ( fabs( xmax - xcenter ) <= Lx/2.0 && fabs( xmin - xcenter ) <= Lx/2.0 && fabs( ymax - ycenter ) <= Ly/2.0 && fabs( ymin - ycenter ) <= Ly/2.0 );
  
}

//-----------------------------------------------------------------------------------------------------------------------------
// Get Mandelstam t
//-----------------------------------------------------------------------------------------------------------------------------

Double_t Get_t(Double_t p1, Double_t th3, Double_t xm1, Double_t xm2)
{
  // Calculate Mandlestam t for elastic scattering,
  // particle 1 (mass xm1) incident on particle 2 (mass xm2)
  // given incident lab momentum p1  and scattering angle th3
  Double_t e1 = sqrt(p1*p1 + xm1*xm1);
  Double_t p = p1;
  Double_t E = e1 + xm2;
  Double_t M = sqrt(E*E - p*p);
  Double_t e3cm = (M*M + xm1*xm1 - xm2*xm2)/(2.0*M);
  Double_t p3cm = sqrt( e3cm*e3cm - xm1*xm1 );
  Double_t sinth3 = sin(th3);
  Double_t costh3 = cos(th3);
  Double_t a = (M*M + xm1*xm1 - xm2*xm2)*p*costh3;
  Double_t b = 2*E*sqrt(M*M*p3cm*p3cm - xm1*xm1*p*p*sinth3*sinth3);
  Double_t c = 2*(M*M + p*p*sinth3*sinth3);
  Double_t p3 = (a + b)/c;
  Double_t e3 = sqrt(p3*p3 + xm1*xm1);
  Double_t e4 = E - e3;
  Double_t t = -2*xm2*(e4 - xm2);
  return t;
}

//-----------------------------------------------------------------------------------------------------------------------------

void PolAnalysisGEnRP( const char *setupfilename="in.txt", Long64_t nentries = 50000, Bool_t writeOut = false ) { 
  
  ifstream setupfile(setupfilename);
  TString filename;
  TChain *C = new TChain("T");
  while( setupfile >> filename )
    C->Add(filename);
  
  Long64_t ntree = C->GetEntries();
  cout << "Total number of generated events = " << ntree << endl;
  genrp_tree *T = new genrp_tree(C);
  if( nentries == -1 ) nentries = ntree;

  Double_t Lumi       = 1e38;
  Double_t Ngen_total = 50000; 
  
  Double_t Eb             = 4.4;  
  Double_t th_bb          = 41.9; 
  Double_t p_bb           = 2.0;  
  Double_t th_sbs         = 24.7;
  Double_t p_sbs          = 3.2;  
  Double_t hcal_dist      = 8.5;
  Double_t hcal_offset    = 0.45;
  Double_t hcal_sampfrac  = 1/18.2;  // ??
  Double_t polana_dist    = 4.15;    // ??
  Double_t polfgem_dist   = 3.9;     // ??
  Double_t polrgem_dist   = 4.4;     // ??
  Double_t polactana_dist = 4.8;     // ??
  Double_t pbbres         = 0.01;    // ??

  Double_t sh_min    = 0.60;
  Double_t sh_max    = 0.80;
  Double_t sh_e      = 0.70;
  Double_t ps_min    = 0.02; 
  Double_t W_min     = 0.0; 
  Double_t W_max     = 2.0;
  Double_t pdiff_cut = 0.2;
  Double_t hcalE_min = 0.07; 
  
  Double_t Mp   = 0.93827;  
  Double_t Mn   = 0.93957;  
  Double_t mu_p = 2.79284734462;

  TFile* fout{NULL};
  TDirectory *dir{NULL};

  if( writeOut) {
    fout = new TFile("Sim_GEnRP.root", "recreate");  
    dir = gDirectory; 
  }

  //-----------------------------------------------------------------------------------------------------------------------------

  const Int_t nphibins = 36;

  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadBottomMargin(.15);

  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleSize(0.055,"X");
  gStyle->SetTitleSize(0.055,"Y");

  gStyle->SetLabelOffset(0.01, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");

  gStyle->SetNdivisions(105,"X");
  gStyle->SetNdivisions(105,"Y");

  gStyle->SetStripDecimals(kFALSE);

  TH1D* hkin_p         = new TH1D("hkin_p","",100,0.5*p_bb,1.25*p_bb);
  TH1D* hkin_th        = new TH1D("hkin_th","",100,-0.3,0.3);
  TH1D* hkin_ph        = new TH1D("hkin_ph","",100,-0.1,0.1);
  TH1D* hkin_x         = new TH1D("hkin_x","",100,-1.0,1.0);
  TH1D* hkin_y         = new TH1D("hkin_y","",100,-0.4,0.4);
  TH1D* hkin_yt        = new TH1D("hkin_yt","",100,-0.15,0.15);
  TH1D* hkin_pdiff     = new TH1D("hkin_pdiff","",50,-0.5,0.5);
  TH1D* hkin_W         = new TH1D("hkin_W","",50,W_min,W_max);
  TH2D* hkin2d_thp     = new TH2D("hkin2d_thp","",100, th_bb-6,th_bb+6.,100,0.5*p_bb,1.25*p_bb);

  TH1D *hbbcal_psE     = new TH1D("hbbcal_psE","",100,0,1.25*p_bb/2.); 
  TH1D *hbbcal_shE     = new TH1D("hbbcal_shE","",100,0,1.25*p_bb); 
  TH1D* hbbcal_cale    = new TH1D("hbbcal_cale","",100,0.5*p_bb,1.25*p_bb);
  TH1D* hbbcal_edivp   = new TH1D("hbbcal_edivp","",100,-1.0,1.0);
  TH2D *hbbcal2d_pss   = new TH2D("hbbcal2d_pssh","",100,0.,1.2, 100,0.,1.2); 

  TH1D* hhcal_e        = new TH1D("hhcal_e","",100,0.0,2.5*p_sbs);
  TH1D* hhcal_x        = new TH1D("hhcal_x","",100,-2.5,2.5);
  TH1D* hhcal_y        = new TH1D("hhcal_y","",100,-2.,2.);
  TH2D* hhcal_xy       = new TH2D("hhcal_xy","",50,-2.,2.,100,-2.5,2.5);
  TH1D* hhcal_prede    = new TH1D("hhcal_prede","",100,0.5*p_sbs,1.5*p_sbs);
  TH1D* hhcal_predx    = new TH1D("hhcal_predx","",100,-2.5,2.5);
  TH1D* hhcal_predy    = new TH1D("hhcal_predy","",100,-2.,2.);
  TH2D* hhcal_predxy   = new TH2D("hhcal_predxy","",50,-2.,2.,50,-2.,2.);
  TH1D* hhcal_deltae   = new TH1D("hhcal_deltae","",100,-4.5,4.5);
  TH1D* hhcal_deltax   = new TH1D("hhcal_deltax","",100,-2.5,2.5);
  TH1D* hhcal_deltay   = new TH1D("hhcal_deltay","",100,-2,2.);
  TH2D* hhcal_deltaxy  = new TH2D("hhcal_deltaxy","",50,-2.,2.,50,-2.5,2.5);

  //-----------------------------------------------------------------------------------------------------------------------------

  TH1D* hkin_pc        = new TH1D("hkin_pc","",100,0.5*p_bb,1.25*p_bb);
  TH1D* hkin_thc       = new TH1D("hkin_thc","",100,-0.3,0.3);
  TH1D* hkin_phc       = new TH1D("hkin_phc","",100,-0.1,0.1);
  TH1D* hkin_xc        = new TH1D("hkin_xc","",100,-1.0,1.0);
  TH1D* hkin_yc        = new TH1D("hkin_yc","",100,-0.4,0.4);
  TH1D* hkin_ytc       = new TH1D("hkin_ytc","",100,-0.15,0.15);
  TH1D* hkin_pdiffc    = new TH1D("hkin_pdiffc","",50,-0.5,0.5);
  TH1D* hkin_Wc        = new TH1D("hkin_Wc","",50,W_min,W_max);
  TH2D* hkin2d_thpc    = new TH2D("hkin2d_thpc","",100, th_bb-6,th_bb+6.,100,0.5*p_bb,1.25*p_bb);

  TH1D *hbbcal_psEc    = new TH1D("hbbcal_psEc","",100,0,1.25*p_bb/2.); 
  TH1D *hbbcal_shEc    = new TH1D("hbbcal_shEc","",100,0,1.25*p_bb); 
  TH1D* hbbcal_calec   = new TH1D("hbbcal_calec","",100,0.5*p_bb,1.25*p_bb);
  TH1D* hbbcal_edivpc  = new TH1D("hbbcal_edivpc","",100,-1.0,1.0);
  TH2D *hbbcal2d_pssc  = new TH2D("hbbcal2d_psshc","",100,0.,1.2, 100,0.,1.2); 

  TH1D* hhcal_ec       = new TH1D("hhcal_ec","",100,0.0,2.5*p_sbs);
  TH1D* hhcal_xc       = new TH1D("hhcal_xc","",100,-2.5,2.5);
  TH1D* hhcal_yc       = new TH1D("hhcal_yc","",100,-2.,2.);
  TH2D* hhcal_xyc      = new TH2D("hhcal_xyc","",50,-2.,2.,100,-2.5,2.5);
  TH1D* hhcal_predec   = new TH1D("hhcal_predec","",100,0.5*p_sbs,1.5*p_sbs);
  TH1D* hhcal_predxc   = new TH1D("hhcal_predxc","",100,-2.5,2.5);
  TH1D* hhcal_predyc   = new TH1D("hhcal_predyc","",100,-2.,2.);
  TH2D* hhcal_predxyc  = new TH2D("hhcal_predxyc","",50,-2.,2.,50,-2.,2.);
  TH1D* hhcal_deltaec  = new TH1D("hhcal_deltaec","",100,-4.5,4.5);
  TH1D* hhcal_deltaxc  = new TH1D("hhcal_deltaxc","",100,-2.5,2.5);
  TH1D* hhcal_deltayc  = new TH1D("hhcal_deltayc","",100,-2.,2.);
  TH2D* hhcal_deltaxyc = new TH2D("hhcal_deltaxyc","",50,-2.5,2.5,50,-2.5,2.5);

  TH1F* hpol_thpp      = new TH1F("hpol_thpp","",50,0,30);
  TH1F* hpol_thppc     = new TH1F("hpol_thppc","",50,0,30);
  TH1F* hpol_phppp     = new TH1F("hpol_phppp","",nphibins,-180,180);
  TH1F* hpol_phppm     = new TH1F("hpol_phppm","",nphibins,-180,180);
  TH1F* hpol_tpp       = new TH1F("hpol_tpp","",50,0,1);
  TH1F* hpol_tppc      = new TH1F("hpol_tppc","",50,0,1);
  TH1F* hpol_sclpp     = new TH1F("hpol_sclpp","",50,0, 0.04);
  TH2F* hpol_zclpp     = new TH2F("hpol_zclpp","",50,3.15, 5.15,50,0,30);

  TH1F* hpol_thnp_cx   = new TH1F("hpol_thnp_cx","",50,0,30);
  TH1F* hpol_thnp_cxc  = new TH1F("hpol_thnp_cxc","",50,0,30);
  TH1F* hpol_phnp_cxp  = new TH1F("hpol_phnp_cxp","",nphibins,-180,180);
  TH1F* hpol_phnp_cxm  = new TH1F("hpol_phnp_cxm","",nphibins,-180,180);
  TH1F* hpol_tnp_cx    = new TH1F("hpol_tnp_cx","",50,0,1);
  TH1F* hpol_tnp_cxc   = new TH1F("hpol_tnp_cxc","",50,0,1);
  TH1F* hpol_sclnp     = new TH1F("hpol_sclnp","",50,0, 0.01);
  TH2F* hpol_zclnp     = new TH2F("hpol_zclnp","",50,3.15, 5.15,50,0,30);

  TH1F* hpol_Ptp       = new TH1F("hpol_Ptp","",50,-0.5,0.5);
  TH1F* hpol_Plp       = new TH1F("hpol_Plp","",50,0.2,1.2);
  TH1F* hpol_Ptn       = new TH1F("hpol_Ptn","",50,0.2,1.2);
  TH1F* hpol_Pln       = new TH1F("hpol_Pln","",50,-0.5,0.5);

  TH1F* hpol_Szp       = new TH1F("hpol_Szp","",50,0.2,1.2);
  TH1F* hpol_Syp       = new TH1F("hpol_Syp","",50,-0.5,0.5);
  TH1F* hpol_Szn       = new TH1F("hpol_Szn","",50,0.2,1.2);
  TH1F* hpol_Syn       = new TH1F("hpol_Syn","",50,-0.5,0.5);
  TH1F* hpol_SxpFP     = new TH1F("hpol_SxpFP","",50,-1.2,1.2);
  TH1F* hpol_SypFP     = new TH1F("hpol_SypFP","",50,-0.5,0.5);

  TH1F* hpol_Chip      = new TH1F("hpol_Chip","",50,0.,90.);

  //-----------------------------------------------------------------------------------------------------------------------------

  TRandom3* fRand = new TRandom3();
  
  TLorentzVector Tp4(0,0,0,Mp); 
  TLorentzVector kp4(0,0,Eb,Eb); 
  TLorentzVector Qp4, kpp4, Rp4; 

  th_bb  *= M_PI/180.;
  th_sbs *= M_PI/180.;

  Int_t helicity{0}, nQEp{0}, nQEn{0}, npp{0}, nnp_cx{0}, nppc{0}, nnp_cxc{0}, npr{0};
  Int_t nprLm{0}, nprRm{0}, nprTm{0}, nprBm{0}, nprLp{0}, nprRp{0}, nprTp{0}, nprBp{0};

  TVector3 SBS_zaxis( -sin(th_sbs), 0, cos(th_sbs) );
  TVector3 SBS_xaxis(0, -1, 0 );
  TVector3 SBS_yaxis = SBS_zaxis.Cross(SBS_xaxis).Unit();
  

  //-----------------------------------------------------------------------------------------------------------------------------

  for(Long64_t ev=0; ev<nentries;ev++) {

    T->GetEntry(ev);
    
    if( ev%100 == 0 )
      cout << ev << endl;

    double weight = T->ev_sigma * T->ev_solang * Lumi / Ngen_total;
    weight = 1;
    
    helicity = T->BeamPol;
    
    //-----------------------------------------------------------------------------------------------------------------------------
    // Pre-cuts -- these are normally in the replay cdef
    //-----------------------------------------------------------------------------------------------------------------------------

    if( T->Earm_BBGEM_Track_ntracks <= 0 ) continue; 
    if( T->Earm_BBPSTF1_det_esum < 0.05 ) continue;

    //-----------------------------------------------------------------------------------------------------------------------------
    // HCAL clustering -- from Andrew's GEN_ERR.C macro
    //-----------------------------------------------------------------------------------------------------------------------------

    TVector3 hcal_hitpos;
    double   hcal_Emax_clust{0};
    
    //Need to do nucleon reconstruction:
    //Find hit with max. energy deposit:
    set<int> unused_hits; //list of hits not yet assigned to a cluster by position in good hit array
    vector<int> hitlist_hcal; //list of hits above threshold by position in tree hit array
    vector<double> xhit_hcal,yhit_hcal,Ehit_hcal,thit_hcal;
    vector<int> row_hcal,col_hcal;
    
    //For island algorithm, we need to be able to look up cells by row and column. So we need
    //a map indexed by unique cell number:
    map<int,int> hitcell_hcal; //key = cell, val = hit index in good hit array:
    double thresh_hcal = 0.007;
    int ngood=0;
    int imax_hcal=-1;
    double emax_hcal = 0.0;
    
    //If we want here we can "cheat" and only include hits attributable to the primary nucleon:
    TVector3 SDTrack_pos; //global position of SD track boundary crossing:
    for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
      double edep = (*(T->Harm_HCalScint_hit_sumedep))[ihit];
      int MID = (*(T->SDTrack_MID))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit]];
      int PrTID = (*(T->PTrack_TID))[(*(T->Harm_HCalScint_hit_ptridx))[ihit]];
      int PID = (*(T->SDTrack_PID))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit]];
      int PrPID = (*(T->PTrack_PID))[(*(T->Harm_HCalScint_hit_ptridx))[ihit]];
            
      if( edep >= thresh_hcal ){
	hitlist_hcal.push_back( ihit );
	xhit_hcal.push_back( (*(T->Harm_HCalScint_hit_xcell))[ihit] );
	yhit_hcal.push_back( (*(T->Harm_HCalScint_hit_ycell))[ihit] );
	Ehit_hcal.push_back( edep );
	thit_hcal.push_back( (*(T->Harm_HCalScint_hit_tavg))[ihit] );
	row_hcal.push_back( (*(T->Harm_HCalScint_hit_row))[ihit] );
	col_hcal.push_back( (*(T->Harm_HCalScint_hit_col))[ihit] );
	hitcell_hcal[(*(T->Harm_HCalScint_hit_cell))[ihit]] = ngood;
	if( edep > emax_hcal ){
	  imax_hcal = ngood;
	  emax_hcal = edep;
	  SDTrack_pos.SetXYZ( (*(T->SDTrack_posx))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit] ],
			      (*(T->SDTrack_posy))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit] ],
			      (*(T->SDTrack_posz))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit] ] );
	}
	unused_hits.insert(ngood);
	ngood++;
      }
    }
    
    int nclust=0;
    vector<int> nhitclust;
    vector<double> xclust, yclust, Eclust,tclust;
    
    if( imax_hcal >= 0 ){
      bool foundclust = true;

      while( foundclust ){
	foundclust=false;
	if( imax_hcal >= 0 ){
	  vector<int> hitsincluster;
	  double xsum = 0.0, ysum = 0.0, esum = 0.0, tsum = 0.0;
	  int nhittemp = 0;
	  
	  tsum += Ehit_hcal[imax_hcal]*thit_hcal[imax_hcal];
	  xsum += Ehit_hcal[imax_hcal]*xhit_hcal[imax_hcal];
	  ysum += Ehit_hcal[imax_hcal]*yhit_hcal[imax_hcal];
	  esum += Ehit_hcal[imax_hcal];
	  
	  nhittemp++;	    
	  unused_hits.erase( imax_hcal );
	  
	  int rowmax = row_hcal[imax_hcal], colmax=col_hcal[imax_hcal];
	  int cellmax = colmax + 12*rowmax;
	  hitsincluster.push_back(imax_hcal);
	    
	  int ihittemp=0;
	  while( ihittemp < nhittemp ){ //island algorithm for clustering hits:
	    //search left, right, up, down until we run out of unused hits
	    //contiguous with any hit already added to the cluster:
	    int rowtemp = row_hcal[hitsincluster[ihittemp]];
	    int coltemp = col_hcal[hitsincluster[ihittemp]];
	    
	    int rowtest=rowtemp;
	    int coltest=coltemp-1;
	    int celltest = coltest + 12*rowtest;
	    
	    std::map<int,int>::iterator foundhit = hitcell_hcal.find(celltest);
	    
	    if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24){ //new hit
	      int hitidx = foundhit->second;
	      if( unused_hits.find(hitidx) != unused_hits.end() ){
		hitsincluster.push_back( hitidx );
		nhittemp++;
		xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		esum += Ehit_hcal[hitidx];
		unused_hits.erase(hitidx);
	      }
	    }
	    
	    coltest=coltemp+1;
	    celltest=coltest + 12*rowtest;
	    foundhit = hitcell_hcal.find(celltest);
	    
	    if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24){ //new hit
	      int hitidx = foundhit->second;
	      if( unused_hits.find(hitidx) != unused_hits.end() ){
		hitsincluster.push_back( hitidx );
		nhittemp++;
		xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		esum += Ehit_hcal[hitidx];
		unused_hits.erase(hitidx);
	      }
	    }
	    
	    coltest=coltemp;
	    rowtest=rowtemp-1;
	    celltest=coltest + 12*rowtest;
	    foundhit = hitcell_hcal.find(celltest);
	    
	    if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24){ //new hit
	      int hitidx = foundhit->second;
	      if( unused_hits.find(hitidx) != unused_hits.end() ){
		hitsincluster.push_back( hitidx );
		nhittemp++;
		xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		esum += Ehit_hcal[hitidx];
		unused_hits.erase(hitidx);
	      }
	    }
	    
	    coltest=coltemp;
	    rowtest=rowtemp+1;
	    celltest=coltest + 12*rowtest;
	    foundhit = hitcell_hcal.find(celltest);
	    
	    if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24 ){ //new hit
	      int hitidx = foundhit->second;
	      if( unused_hits.find(hitidx) != unused_hits.end() ){
		hitsincluster.push_back( hitidx );
		nhittemp++;
		xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		esum += Ehit_hcal[hitidx];
		unused_hits.erase(hitidx);
	      }
	    }
	    ihittemp++; //go to next available hit.
	  } //while (ihittemp < nhittemp)
	  
	  foundclust = true;
	  nhitclust.push_back( nhittemp );
	  xclust.push_back( xsum/esum );
	  yclust.push_back( ysum/esum );
	  tclust.push_back( tsum/esum );
	  Eclust.push_back( esum );
	  nclust++;
	  
	} //if (imax_hcal >= 0 ) subsequent maxima
	
	imax_hcal = -1;
	emax_hcal = 0.0;
	//Loop over unused hits and see if we can find another maximum:
	for(set<int>::iterator ihit=unused_hits.begin(); ihit!=unused_hits.end(); ++ihit ){
	  int hitidx = *ihit;
	  if( Ehit_hcal[hitidx] > emax_hcal ){
	    emax_hcal = Ehit_hcal[hitidx];
	    imax_hcal = hitidx;
	  }
	}
      } //while (foundclust)
    } //if(imax_hcal >= 0) first maximum
    
    if( nclust > 0 ){ //choose cluster with maximum energy:
      int imax_clust=0;
      hcal_Emax_clust = 0.0;
      for( int iclust=0; iclust<nclust; iclust++ ){
	if( Eclust[iclust] > hcal_Emax_clust ){
	  imax_clust = iclust;
	  hcal_Emax_clust = Eclust[iclust];
	}
      }
      hcal_hitpos.SetXYZ(xclust[imax_clust],yclust[imax_clust],hcal_dist);
    }
    else
      hcal_hitpos.SetXYZ(-9999., -9999., -9999.);
    
    //-----------------------------------------------------------------------------------------------------------------------------
    // Kinematic reconstruction
    //-----------------------------------------------------------------------------------------------------------------------------
  
    Double_t p = (*(T->Earm_BBGEM_Track_P))[0];
    p += fRand->Gaus(0,pbbres*p);

    TVector3 pvect( p*sin(T->ev_th)*cos(T->ev_ph), p*sin(T->ev_th)*sin(T->ev_ph), p*cos(T->ev_th) ); 

    TVector3 BB_zaxis( sin(th_bb), 0.0, cos(th_bb) ); //BB is on beam left, global x axis points to beam left 
    TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down: 
    TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit(); 

    double tracker_pitch_angle = 10. * M_PI/180.;
    
    TVector3 BBGEM_zaxis = BB_zaxis; 
    TVector3 BBGEM_yaxis = BB_yaxis; 
    BBGEM_zaxis.Rotate( -tracker_pitch_angle, BBGEM_yaxis ); 
    TVector3 BBGEM_xaxis = BBGEM_yaxis.Cross(BBGEM_zaxis).Unit(); 
    
    TVector3 pvect_BB( pvect.Dot(BB_xaxis), pvect.Dot(BB_yaxis), pvect.Dot(BB_zaxis) ); 
    
    double xptar = pvect_BB.X()/pvect_BB.Z(); 
    double yptar = pvect_BB.Y()/pvect_BB.Z(); 
       
    double vx = T->ev_vx; 
    double vy = T->ev_vy; 
    double vz = T->ev_vz; 
    
    TVector3 vertex(vx,vy,vz); 
       
    TVector3 vertex_BB(vertex.Dot(BB_xaxis), vertex.Dot(BB_yaxis), vertex.Dot(BB_zaxis) ); 
 
    double ytar = vertex_BB.Y() - yptar * vertex_BB.Z(); 
    double xtar = vertex_BB.X() - xptar * vertex_BB.Z(); 

    Double_t th      = acos( pvect.Z() / (*(T->Earm_BBGEM_Track_P))[0] );
    Double_t pexp_th = Eb / (1.0 + Eb/Mp*(1.0-cos(th)));
    Double_t pdiff   = p - pexp_th; 
    
    Double_t px = pvect.X();
    Double_t py = pvect.Y();
    Double_t pz = pvect.Z();
    
    kpp4.SetPxPyPzE(px,py,pz,p);
    Qp4 = kp4 - kpp4;
    Rp4 = Tp4 + Qp4;
    Double_t W2 = Rp4.M2(); 
    
    hkin_p->Fill(p,weight);
    hkin_x->Fill((*(T->Earm_BBGEM_Track_X))[0],weight);
    hkin_y->Fill((*(T->Earm_BBGEM_Track_Y))[0],weight);
    hkin_th->Fill( xptar ,weight);
    hkin_ph->Fill( yptar ,weight);
    hkin_yt->Fill(ytar,weight);
    hkin_pdiff->Fill( pdiff ,weight);
    hkin_W->Fill(W2,weight);
    hkin2d_thp->Fill((180/M_PI)*th, p,weight);

    //-----------------------------------------------------------------------------------------------------------------------------

    hbbcal_psE->Fill( T->Earm_BBPSTF1_det_esum  , weight);
    hbbcal_shE->Fill( T->Earm_BBSHTF1_det_esum  , weight);
    hbbcal_cale->Fill(T->Earm_BBPSTF1_det_esum  + T->Earm_BBSHTF1_det_esum , weight);
    hbbcal_edivp->Fill( T->Earm_BBPSTF1_det_esum  +  T->Earm_BBSHTF1_det_esum  - p , weight);
    hbbcal2d_pss->Fill( T->Earm_BBSHTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0],  T->Earm_BBPSTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0] , weight);

    //-----------------------------------------------------------------------------------------------------------------------------
    // Project reconstructed recoil nucleon p4 to SBS detectors
    //-----------------------------------------------------------------------------------------------------------------------------
    
    Rp4.RotateY(th_sbs);
    
    Double_t recoil_th    = TMath::ATan(Rp4.Px()/Rp4.Pz());
    Double_t recoil_ph    = TMath::ATan(Rp4.Py()/Rp4.Pz());
    
    Double_t hcal_x       = ( hcal_hitpos.Y() );
    Double_t hcalpred_x   = (hcal_dist * TMath::Sin( recoil_ph )) + hcal_offset;
    Double_t polanapred_x = (polana_dist * TMath::Sin( recoil_th ));
    
    Double_t hcal_y       = ( hcal_hitpos.X() );
    Double_t hcalpred_y   = -hcal_dist * TMath::Sin( recoil_th );
    Double_t polanapred_y = polana_dist * TMath::Sin( recoil_ph );
    
    Double_t hcaldelta_x  = hcal_x - hcalpred_x;
    Double_t hcaldelta_y  = hcal_y + hcalpred_y;
    
    hhcal_x->Fill(hcal_x, weight);
    hhcal_y->Fill(hcal_y, weight);
    hhcal_xy->Fill(hcal_y,hcal_x, weight);
    hhcal_e->Fill( hcal_Emax_clust * 1./hcal_sampfrac , weight);

    hhcal_predx->Fill(hcalpred_x, weight);
    hhcal_predy->Fill(hcalpred_y, weight);
    hhcal_predxy->Fill(hcalpred_y,hcalpred_x, weight);
    hhcal_prede->Fill( Rp4.T() , weight);

    hhcal_deltax->Fill(hcaldelta_x, weight);
    hhcal_deltay->Fill(hcaldelta_y, weight);
    hhcal_deltaxy->Fill(hcaldelta_y,hcaldelta_x, weight);
    hhcal_deltae->Fill( hcal_Emax_clust * 1./hcal_sampfrac - Rp4.T() , weight);

    //-----------------------------------------------------------------------------------------------------------------------------
    // Apply cuts
    //-----------------------------------------------------------------------------------------------------------------------------
    
    if( ApplyElec ) { 
      if(  T->Earm_BBPSTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0]  < ps_min ) continue;
      if( (T->Earm_BBPSTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0]  + sh_e*T->Earm_BBSHTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0]  )< sh_min ) continue;
      if( (T->Earm_BBPSTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0] + sh_e* T->Earm_BBSHTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0] )> sh_max) continue;
      if( (*(T->Earm_BBGEM_Track_P))[0] < 0.5*p_bb ) continue;
      if( W2 < W_min ) continue; 
      if( W2 > W_max ) continue;   
    }
    if( ApplyElas ) {
      if( fabs(pdiff) > pdiff_cut ) continue;
      if( hcal_Emax_clust < hcalE_min ) continue;
    }

    //-----------------------------------------------------------------------------------------------------------------------------
    
    hkin_pc->Fill(p, weight);
    hkin_xc->Fill((*(T->Earm_BBGEM_Track_X))[0], weight);
    hkin_yc->Fill((*(T->Earm_BBGEM_Track_Y))[0], weight);
    hkin_thc->Fill( xptar , weight);
    hkin_phc->Fill( yptar , weight);
    hkin_ytc->Fill(ytar, weight);
    hkin_pdiffc->Fill( pdiff , weight);
    hkin_Wc->Fill(W2, weight);
    hkin2d_thpc->Fill((180./M_PI)*th, p, weight);

    //-----------------------------------------------------------------------------------------------------------------------------

    hbbcal_psEc->Fill( T->Earm_BBPSTF1_det_esum  , weight);
    hbbcal_shEc->Fill( T->Earm_BBSHTF1_det_esum  , weight);
    hbbcal_calec->Fill(T->Earm_BBPSTF1_det_esum  + T->Earm_BBSHTF1_det_esum , weight);
    hbbcal_edivpc->Fill( T->Earm_BBPSTF1_det_esum  +  T->Earm_BBSHTF1_det_esum  - p , weight);
    hbbcal2d_pssc->Fill( T->Earm_BBSHTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0],  T->Earm_BBPSTF1_det_esum/(*(T->Earm_BBGEM_Track_P))[0] , weight);

    //-----------------------------------------------------------------------------------------------------------------------------

    hhcal_ec->Fill( hcal_Emax_clust * 1./hcal_sampfrac , weight);
    hhcal_xc->Fill(hcal_x, weight);
    hhcal_yc->Fill(hcal_y, weight);
    hhcal_xyc->Fill(-hcal_y,hcal_x, weight);

    hhcal_predec->Fill( Rp4.T() , weight);
    hhcal_predxc->Fill(hcalpred_x, weight);
    hhcal_predyc->Fill(hcalpred_y, weight);
    hhcal_predxyc->Fill(hcalpred_y,hcalpred_x, weight);

    hhcal_deltaec->Fill( hcal_Emax_clust *  1./hcal_sampfrac  - Rp4.T() , weight);
    hhcal_deltaxc->Fill(hcaldelta_x, weight);
    hhcal_deltayc->Fill(hcaldelta_y, weight);
    hhcal_deltaxyc->Fill(hcaldelta_y,hcaldelta_x, weight);
    
    //-----------------------------------------------------------------------------------------------------------------------------
    // Polarimeter Analysis
    //-----------------------------------------------------------------------------------------------------------------------------
    
    if( T->ev_fnucl == 0 ) 
      nQEn++;
    else if( T->ev_fnucl == 1 ) 
      nQEp++;
    
    Double_t anaf_x{0}, anaf_y{0}, anar_x{0}, anar_y{0}, polhcal_x{0}, polhcal_y{0};
    Double_t angz{0}, thsc{0}, phsc{0};
    
    TVector3 front_vtx(0.,0.,0);
    TVector3 anaf_vtx(0.,0.,0);
    TVector3 anar_vtx(0.,0.,0);
    TVector3 rear_vtx(0.,0.,0);
    TVector3 hcal_vtx(0.,0.,0);
    TVector3 in_track(0.,0.,0);
    TVector3 out_track(0.,0.,0);

    //-----------------------------------------------------------------------------------------------------------------------------
    // pN -> pN in passive (steel) analyser
    //-----------------------------------------------------------------------------------------------------------------------------
    
    // Still to be done:
    // Choose front track based on BB reconstruction
    // Choose rear track based on matching HCAL max cluster

    if( T->Harm_CEPolFront_Track_ntracks == 1  && T->Harm_CEPolRear_Track_ntracks == 1 ) {
      
      if( helicity == 1 && T->ev_fnucl == 1 ) {
	
	hpol_Ptp->Fill( T->ev_Pt );
	hpol_Plp->Fill( T->ev_Pl );

	hpol_Szp->Fill( T->ev_Sz );
	hpol_Syp->Fill( T->ev_Sy );

 	hpol_SxpFP->Fill( (*(T->Harm_CEPolFront_Track_Sx))[0] );
 	hpol_SypFP->Fill( (*(T->Harm_CEPolFront_Track_Sy))[0] );
      }
      
      anaf_x =  (*(T->Harm_CEPolFront_Track_X))[0] + (polana_dist - polfgem_dist) * (*(T->Harm_CEPolFront_Track_Xp))[0];
      anaf_y =  (*(T->Harm_CEPolFront_Track_Y))[0] + (polana_dist - polfgem_dist) * (*(T->Harm_CEPolFront_Track_Yp))[0];
      
      anar_x =  (*(T->Harm_CEPolRear_Track_X))[0] + (polana_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Xp))[0];
      anar_y =  (*(T->Harm_CEPolRear_Track_Y))[0] + (polana_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Yp))[0];
      
      front_vtx.SetXYZ( (*(T->Harm_CEPolFront_Track_X))[0], (*(T->Harm_CEPolFront_Track_Y))[0], polfgem_dist );
      anaf_vtx.SetXYZ( anaf_x, anaf_y, polana_dist );
      anar_vtx.SetXYZ( anar_x, anar_y, polana_dist );
      rear_vtx.SetXYZ( (*(T->Harm_CEPolRear_Track_X))[0], (*(T->Harm_CEPolRear_Track_Y))[0], polrgem_dist );
      
      polhcal_x =  (*(T->Harm_CEPolRear_Track_X))[0] + (hcal_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Xp))[0];
      polhcal_y =  (*(T->Harm_CEPolRear_Track_Y))[0] + (hcal_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Yp))[0];
      
      hcal_vtx.SetXYZ( polhcal_x, polhcal_y, hcal_dist );
            
      in_track.SetXYZ( tan ((anaf_vtx.X() - front_vtx.X())/(anaf_vtx.Z() - front_vtx.Z())) ,
		       tan ((anaf_vtx.Y() - front_vtx.Y())/(anaf_vtx.Z() - front_vtx.Z())) ,
		       1. );
      
      out_track.SetXYZ( tan( (rear_vtx.X() - anar_vtx.X())/(rear_vtx.Z() - anar_vtx.Z())),
			tan ((rear_vtx.Y() - anar_vtx.Y())/(rear_vtx.Z() - anar_vtx.Z())),
			1. );
      
      in_track  = in_track.Unit();
      out_track = out_track.Unit();
      
      TVector3 yaxistemp(0,1,0);
      TVector3 xaxistemp = yaxistemp.Cross(in_track).Unit();
      yaxistemp = in_track.Cross(xaxistemp).Unit();
      
      thsc = 180./M_PI * acos( out_track.Dot(in_track) ); 
      phsc = 180./M_PI * TMath::ATan2( out_track.Dot(yaxistemp), out_track.Dot(xaxistemp) );
      
      double sclose,zclose;
      
      calc_sclose_zclose( front_vtx, rear_vtx, in_track, out_track, sclose, zclose );
      
      hpol_sclpp->Fill( sclose , weight); 
      hpol_zclpp->Fill( zclose , thsc, weight); 

      hpol_thpp->Fill ( thsc , weight);

      Double_t t = fabs(Get_t(Rp4.Mag(), thsc * M_PI/180, Mp, Mp));
      hpol_tpp->Fill( t , weight);
      
      bool conetestpp = conetest( front_vtx, in_track, thsc, zclose, polrgem_dist );

      TVector3 ppvect_global(T->ev_npx, T->ev_npy, T->ev_npz );
      TVector3 ppunit_global = ppvect_global.Unit();
      
      TVector3 ppunit_sbs( ppunit_global.Dot( SBS_xaxis ),
			   ppunit_global.Dot( SBS_yaxis ),
			   ppunit_global.Dot( SBS_zaxis ) );

      double thetabend = acos( in_track.Dot( ppunit_sbs ) );
      double pp        = T->ev_np;
      double gamma     = sqrt(1.+pow(pp/Mp,2));
      double chi       = gamma*(mu_p - 1.0)*thetabend;
      double sinchi    = sin(chi);

      hpol_Chip->Fill( chi * 180/M_PI );

      if( thsc >= 5 && sclose < 0.1 &&  fabs(phsc)<180 ) { 
     	nppc++;
	hpol_thppc->Fill ( thsc , weight);
	hpol_tppc->Fill( t , weight);
	if( helicity == -1 ) hpol_phppm->Fill ( phsc , weight);
	else if( helicity == 1 ) hpol_phppp->Fill ( phsc , weight);
      }
      //	cout << endl << "New Event " << endl;
      // 	cout << front_vtx->X() << "\t" << front_vtx->Y() << "\t" << front_vtx->Z() << endl;
      // 	cout << anaf_vtx->X() << "\t" << anaf_vtx->Y() << "\t" << anaf_vtx->Z() << endl;
      // 	cout << anar_vtx->X() << "\t" << anar_vtx->Y() << "\t" << anar_vtx->Z() << endl;
      // 	cout << rear_vtx->X() << "\t" << rear_vtx->Y() << "\t" << rear_vtx->Z() << endl;
      // 	cout << hcal_vtx->X() << "\t" << hcal_vtx->Y() << "\t" << hcal_vtx->Z() << endl;
      // 	cout << -hcal_hitpos.Y()-hcal_offset << "\t" << hcal_hitpos.X() << "\t" << hcal_hitpos.Z() << endl;
      
      // 	in_track.SetXYZ( (*(T->Harm_CEPolFront_Track_Xp))[0],
      // 			 (*(T->Harm_CEPolFront_Track_Yp))[0],
      // 			 1.0 );
      
      // 	out_track.SetXYZ( (*(T->Harm_CEPolRear_Track_Xp))[0],
      // 			  (*(T->Harm_CEPolRear_Track_Yp))[0],
      // 			  1.0 );
      
    }

    //-----------------------------------------------------------------------------------------------------------------------------
    // np -> pn ChEx in passive (steel) analyser, coordinate system = SBS TRANSPORT
    //-----------------------------------------------------------------------------------------------------------------------------

    // Still to be done:
    // Compare BB projected front track to analyzer hit position
    // Choose rear track based on matching HCAL max cluster
    
    if( T->Harm_CEPolFront_Track_ntracks == 0 && T->Harm_CEPolRear_Track_ntracks == 1 ) {
      
      if( helicity == 1  && T->ev_fnucl == 0 ) {

	hpol_Ptn->Fill( T->ev_Pt );
	hpol_Pln->Fill( T->ev_Pl );

	hpol_Szn->Fill( T->ev_Sz );
	hpol_Syn->Fill( T->ev_Sy );
      }
      
      anar_x    =  (*(T->Harm_CEPolRear_Track_X))[0] + (polana_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Xp))[0];
      anar_y    =  (*(T->Harm_CEPolRear_Track_Y))[0] + (polana_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Yp))[0];
      polhcal_x =  (*(T->Harm_CEPolRear_Track_X))[0] + (hcal_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Xp))[0];
      polhcal_y =  (*(T->Harm_CEPolRear_Track_Y))[0] + (hcal_dist - polrgem_dist) * (*(T->Harm_CEPolRear_Track_Yp))[0];
      
      front_vtx = vertex;
      front_vtx.RotateY( th_sbs );
      anaf_vtx.SetXYZ( polanapred_x, polanapred_y, polana_dist ); 
      anar_vtx.SetXYZ( anar_x, anar_y, polana_dist );
      rear_vtx.SetXYZ( (*(T->Harm_CEPolRear_Track_X))[0], (*(T->Harm_CEPolRear_Track_Y))[0], polrgem_dist );
      hcal_vtx.SetXYZ( -hcal_hitpos.Y()-hcal_offset, hcal_hitpos.X(), hcal_hitpos.Z() );
      
      in_track.SetXYZ( tan ((anar_vtx.X() - front_vtx.X())/(anar_vtx.Z() - front_vtx.Z())) ,
		       tan ((anar_vtx.Y() - front_vtx.Y())/(anar_vtx.Z() - front_vtx.Z())) ,
		       1. );
      
      out_track.SetXYZ( tan( (rear_vtx.X() - anar_vtx.X())/(rear_vtx.Z() - anar_vtx.Z())),
			tan ((rear_vtx.Y() - anar_vtx.Y())/(rear_vtx.Z() - anar_vtx.Z())),
			1. );
      
      in_track  = in_track.Unit();
      out_track = out_track.Unit();
      
      TVector3 yaxistemp(0,1,0);
      TVector3 xaxistemp = yaxistemp.Cross(in_track).Unit();
      yaxistemp = in_track.Cross(xaxistemp).Unit();
      
      thsc = 180./M_PI * acos( out_track.Dot(in_track) ); 
      phsc = 180./M_PI * TMath::ATan2( out_track.Dot(yaxistemp), out_track.Dot(xaxistemp) );
      
      double sclose,zclose;
      
      calc_sclose_zclose( front_vtx, rear_vtx, in_track, out_track, sclose, zclose );
      
      hpol_sclnp->Fill( sclose , weight); 
      hpol_zclnp->Fill( zclose , thsc, weight); 

      Double_t t = fabs(Get_t(Rp4.Mag(), thsc * M_PI/180, Mp, Mp));

      hpol_thnp_cx->Fill ( thsc , weight);
      hpol_tnp_cx->Fill( t , weight);
      
      bool conetestnp = conetest( anaf_vtx, in_track, thsc, zclose, polrgem_dist );

      if( thsc >= 0 && sclose < 0.2 && fabs(phsc)<140 ) { 
	nnp_cxc++;
	
	hpol_thnp_cxc->Fill ( thsc , weight);
	hpol_tnp_cxc->Fill( t, weight );
      
	if( helicity == -1 ) hpol_phnp_cxm->Fill ( phsc , weight);
	else if( helicity == 1 ) hpol_phnp_cxp->Fill ( phsc , weight);
      }

      //	cout << endl << "New Event " << endl;
      //        cout << front_vtx.X() << "\t" << front_vtx.Y() << "\t" << front_vtx.Z() << endl;
      // 	cout << anaf_vtx.X() << "\t" << anaf_vtx.Y() << "\t" << anaf_vtx.Z() << endl;
      // 	cout << anar_vtx.X() << "\t" << anar_vtx.Y() << "\t" << anar_vtx.Z() << endl;
      // 	cout << rear_vtx.X() << "\t" << rear_vtx.Y() << "\t" << rear_vtx.Z() << endl;
      // 	cout << hcal_vtx.X() << "\t" << hcal_vtx.Y() << "\t" << hcal_vtx.Z() << endl;
      // 	cout << polhcal_x << "\t" << polhcal_y << "\t" << hcal_dist << endl;

    }

  }

 
  //-----------------------------------------------------------------------------------------------------------------------------
  // Kinematic plots
  //-----------------------------------------------------------------------------------------------------------------------------

  TLatex* tex;

  if( PlotKine ) {

    TCanvas* ckin1 = new TCanvas("ckin1","",1200,800);
    ckin1->Divide(4,2);

    ckin1->cd(1);
    hkin_p->Draw();
    hkin_p->GetXaxis()->SetTitle("p_{bbtrack} [GeV/c]");
    hkin_pc->Draw("same");
    hkin_pc->SetLineColor(2);

    tex = new TLatex( 0.25, 0.8, "All tracks");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(4);
    tex->SetTextSize(0.05);
    tex->Draw();

    tex = new TLatex( 0.25, 0.75, "After cuts");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.05);
    tex->Draw();

    ckin1->cd(2);
    hkin_th->Draw();
    hkin_th->GetXaxis()->SetTitle("#theta_tgt_{bbtrack} [rad]");
    hkin_thc->Draw("same");
    hkin_thc->SetLineColor(2);

    ckin1->cd(3);
    hkin_ph->Draw();
    hkin_ph->GetXaxis()->SetTitle("#phi_tgt_{bbtrack} [rad]");
    hkin_phc->Draw("same");
    hkin_phc->SetLineColor(2);

    ckin1->cd(4);
    hkin_yt->Draw();
    hkin_yt->GetXaxis()->SetTitle("y_tgt_{bbtrack} [m]");
    hkin_ytc->Draw("same");
    hkin_ytc->SetLineColor(2);

    ckin1->cd(5);
    hkin_W->Draw();
    hkin_W->GetXaxis()->SetTitle("W^{2} [GeV^{2}]");
    hkin_Wc->Draw("same");
    hkin_Wc->SetLineColor(2);
    
    ckin1->cd(6);
    hkin_pdiff->Draw();
    hkin_pdiff->GetXaxis()->SetTitle("p_{bbtrack} - p_{#theta exp} [GeV/c]");
    hkin_pdiffc->Draw("same");
    hkin_pdiffc->SetLineColor(2);

    ckin1->cd(7);
    hkin2d_thpc->Draw("colz");
    hkin2d_thpc->GetXaxis()->SetTitle("#theta_{bblab} [degrees]");
    hkin2d_thpc->GetYaxis()->SetTitle("p_{bbtrack} [GeV/c]");

    ckin1->cd(8);
    hhcal_deltaxyc->Draw("colz");
    hhcal_deltaxyc->GetXaxis()->SetTitle("HCal (Meas - Pred) y[m]");
    hhcal_deltaxyc->GetYaxis()->SetTitle("HCal (Meas - Pred) x[m]");
    
    ckin1->Print("temp_aalltracks.pdf");
  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  // Pre-shower and Shower correlation plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotBBCal ) {
    TCanvas* cbbcal = new TCanvas("cbbcal","",1200,800);
    cbbcal->Divide(3,2);
    cbbcal->cd(1);
    hbbcal_psE->Draw("");
    hbbcal_psEc->SetLineColor(2);
    hbbcal_psEc->Draw("same");
    hbbcal_psE->GetXaxis()->SetTitle("BBCal PS Energy [GeV]");

    cbbcal->cd(2);
    hbbcal_shE->Draw("");
    hbbcal_shEc->SetLineColor(2);
    hbbcal_shEc->Draw("same");
    hbbcal_shE->GetXaxis()->SetTitle("BBCal Shower Energy [GeV]");

    cbbcal->cd(3);
    hbbcal_cale->Draw();
    hbbcal_calec->SetLineColor(2);
    hbbcal_calec->Draw("same");
    hbbcal_cale->GetXaxis()->SetTitle("E_{bbcal} [GeV]");

    cbbcal->cd(4);
    hbbcal_edivp->Draw();
    hbbcal_edivpc->SetLineColor(2);
    hbbcal_edivpc->Draw("same");
    hbbcal_edivp->GetXaxis()->SetTitle("(E_{bbcal} - p_{bbtrack}) [GeV]");

    cbbcal->cd(5);
    hbbcal2d_pss->Draw("colz");
    hbbcal2d_pss->GetXaxis()->SetTitle("E_{SH}/p_{bbtrack}");
    hbbcal2d_pss->GetYaxis()->SetTitle("E_{PS}/p_{bbtrack}");
    TLine *line1 = new TLine(0,sh_min,sh_min/sh_e,0); 
    line1->SetLineColor(2); 
    line1->SetLineWidth(2); 
    line1->Draw(); 
    TLine *line3 = new TLine(0,sh_max,sh_max/sh_e,0); 
    line3->SetLineColor(2); 
    line3->SetLineWidth(2); 
    line3->Draw(); 
    TLine *line2 = new TLine(0,ps_min,1.2,ps_min); 
    line2->SetLineColor(2); 
    line2->SetLineWidth(2); 
    line2->Draw(); 

    cbbcal->cd(6);
    hbbcal2d_pssc->Draw("colz");
    hbbcal2d_pssc->GetXaxis()->SetTitle("E_{SH}/p_{bbtrack}");
    hbbcal2d_pssc->GetYaxis()->SetTitle("E_{PS}/p_{bbtrack}");

    cbbcal->Print("temp_bbcal.pdf");
  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  // HCal plots
  //-----------------------------------------------------------------------------------------------------------------------------

  if( PlotHCal ) {

    TCanvas* chcal = new TCanvas("chcal","",1200,800);
    chcal->Divide(4,3);
    chcal->cd(1);
    hhcal_x->Draw("");
    hhcal_x->GetXaxis()->SetTitle("HCal x [m]");
    hhcal_xc->SetLineColor(2);
    hhcal_xc->Draw("same");
    
    chcal->cd(2);
    hhcal_y->Draw("");
    hhcal_y->GetXaxis()->SetTitle("HCal y [m]");
    hhcal_yc->SetLineColor(2);
    hhcal_yc->Draw("same");
    
    chcal->cd(3)->SetLogy(1);
    hhcal_e->Draw("");
    hhcal_e->GetXaxis()->SetTitle("HCal E [GeV]");
    hhcal_ec->SetLineColor(2);
    hhcal_ec->Draw("same");

    chcal->cd(4);
    hhcal_xyc->Draw("colz");
    hhcal_xyc->GetXaxis()->SetTitle("HCal y [m]");
    hhcal_xyc->GetYaxis()->SetTitle("HCal x [m]");
    
    chcal->cd(5);
    hhcal_predx->Draw("");
    hhcal_predx->GetXaxis()->SetTitle("HCal predicted x [m]");
    hhcal_predxc->SetLineColor(2);
    hhcal_predxc->Draw("same");
    
    chcal->cd(6);
    hhcal_predy->Draw("");
    hhcal_predy->GetXaxis()->SetTitle("HCal predicted y [m]");
    hhcal_predyc->SetLineColor(2);
    hhcal_predyc->Draw("same");
    
    chcal->cd(7);
    hhcal_prede->Draw("");
    hhcal_prede->GetXaxis()->SetTitle("HCal predicted E [GeV]");
    hhcal_predec->SetLineColor(2);
    hhcal_predec->Draw("same");

    chcal->cd(8);
    hhcal_predxyc->Draw("colz");
    hhcal_predxyc->GetXaxis()->SetTitle("HCal Pred y [m]");
    hhcal_predxyc->GetYaxis()->SetTitle("HCal Pred x [m]");
   
    chcal->cd(9);
    hhcal_deltax->Draw("");
    hhcal_deltax->GetXaxis()->SetTitle("HCal (meas - pred) x [m]");
    hhcal_deltaxc->SetLineColor(2);
    hhcal_deltaxc->Draw("same");
    
    chcal->cd(10);
    hhcal_deltay->Draw("");
    hhcal_deltay->GetXaxis()->SetTitle("HCal (meas - pred) y [m]");
    hhcal_deltayc->SetLineColor(2);
    hhcal_deltayc->Draw("same");
    
    chcal->cd(11);
    hhcal_deltae->Draw("");
    hhcal_deltae->GetXaxis()->SetTitle("HCal (meas -pred) E [GeV]");
    hhcal_deltaec->SetLineColor(2);
    hhcal_deltaec->Draw("same");
    
    chcal->cd(12);
    hhcal_deltaxyc->Draw("colz");
    hhcal_deltaxyc->GetXaxis()->SetTitle("HCal (meas - pred) y [m]");
    hhcal_deltaxyc->GetYaxis()->SetTitle("HCal (meas - pred) x [m]");
    
    chcal->Print("temp_hcal.pdf");
  }

  //-----------------------------------------------------------------------------------------------------------------------------
  // Polarimetry plots
  //-----------------------------------------------------------------------------------------------------------------------------
  
  if( PlotPol ) {

    TCanvas* cpolpp = new TCanvas("cpolpp","",1200,800);
    cpolpp->Divide(4,2);

    cpolpp->cd(1);

    tex = new TLatex( 0.08, 0.85, "pN -> pN (steel)");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();

    Double_t effpp = (Float_t)nppc/nQEp;
    tex = new TLatex( 0.08, 0.7, Form("Efficiency = %3.2f %%", 100.*effpp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    Double_t Ptp = hpol_Syp->GetMean();
    Double_t Plp = hpol_Szp->GetMean();
    tex = new TLatex( 0.08, 0.6, Form("Pt = %3.2f (#rightarrow Py^{FP})", Ptp)); 
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    tex = new TLatex( 0.08, 0.5, Form("Pt = %3.2f (#rightarrow Px^{FP})", Plp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    cpolpp->cd(2);
    hpol_sclpp->Draw("");
    hpol_sclpp->GetXaxis()->SetTitle("sclose [m]");

    cpolpp->cd(3);
    hpol_zclpp->Draw("colz");
    hpol_zclpp->GetXaxis()->SetTitle("zclose [m]");

    cpolpp->cd(4)->SetLogy(1);
    hpol_tppc->Draw("");
    hpol_tppc->GetXaxis()->SetTitle("|t| [GeV]");

    cpolpp->cd(5)->SetLogy(1);
    hpol_thpp->Draw("");
    hpol_thpp->GetXaxis()->SetTitle("#theta [deg]");

    hpol_thppc->Draw("same");
    hpol_thppc->SetLineColor(2);

    cpolpp->cd(6);
    hpol_phppm->Sumw2();
    hpol_phppm->Draw("");
    hpol_phppm->SetMarkerStyle(21);
    hpol_phppm->SetMarkerColor(1);
    hpol_phppm->SetLineColor(1);
    hpol_phppm->SetMarkerSize(0.5);
    hpol_phppm->GetXaxis()->SetTitle("#phi [deg]");

    hpol_phppp->Sumw2();
    hpol_phppp->SetMarkerStyle(21);
    hpol_phppp->SetMarkerColor(6);
    hpol_phppp->SetLineColor(6);
    hpol_phppp->SetMarkerSize(0.5);
    hpol_phppp->Draw("same");

    cpolpp->cd(7);
    TH1F* hsum  = new TH1F("hsum","",nphibins,-180,180);
    TH1F* hdiff = new TH1F("hdiff","",nphibins,-180,180);
    TH1F* hasym = new TH1F("hasym","",nphibins,-180,180);
    hsum->Sumw2();
    hdiff->Sumw2();
    hasym->Sumw2();
    hsum->Add( hpol_phppp, hpol_phppm );
    hdiff->Add( hpol_phppp, hpol_phppm, 1, -1 );
    hasym->Divide( hdiff, hsum, 1, 1 );
//     Double_t n2plus  = hpol_phppp->Integral(0,nphibins);
//     n2plus = nphibins/(n2plus);
//     Double_t n2minus = hpol_phppm->Integral(0,nphibins);
//     n2minus = nphibins/(n2minus);
//     hasym->Add(hpol_phppp,hpol_phppm,n2plus,-n2minus);
    hasym->Draw("E1");
    hasym->SetMarkerStyle(21);
    hasym->SetMarkerColor(kBlue);
    hasym->SetLineColor(kBlue);
    hasym->SetMarkerSize(0.5);
    hasym->GetXaxis()->SetTitle("#varphi (degrees)");
    hasym->GetYaxis()->SetTitle("Asymmetry");
    hasym->GetYaxis()->SetRangeUser(-1.,1.);
    
    TF1 *fit = new TF1("fit", "[0]*sin(x/57.29578)-[1]*cos(x/57.29578)",-180,180);
    hasym->Fit(fit,"","",-180,180);

    Double_t Pypp  = fit->GetParameter(1);
    Double_t ePypp = fit->GetParError(1);
    Double_t Pxpp  = fit->GetParameter(0);
    Double_t ePxpp = fit->GetParError(0);

    Double_t Chipp = hpol_Chip->GetMean();

    cpolpp->cd(8);
    Double_t tpp[]   = { -1, hpol_tppc->GetMean(), 2 };
    Double_t etpp[]  = { 0, 0, 0, 0 };
    
    Double_t Aypp[]  = { -3, fabs(Pxpp/Plp/sin(Chipp*M_PI/180)), -3 };
    Double_t eAypp[] = { 0, fabs(ePxpp/Plp/sin(Chipp*M_PI/180)), 0 };

    TGraphErrors* gAypp = new TGraphErrors( 3, tpp, Aypp, etpp, eAypp );
    gAypp->SetMarkerColor( 4 );
    gAypp->SetLineColor( 4 );
    gAypp->SetMarkerStyle( 20 );
    gAypp->SetMarkerSize( 1.5 );
    gAypp->Draw("AP");
    gAypp->GetXaxis()->SetRangeUser( 0, 0.15 );
    gAypp->GetYaxis()->SetRangeUser( 0, 1.0 );
    gAypp->GetXaxis()->SetTitle( "|t| [GeV]");
    gAypp->GetYaxis()->SetTitle( "Ay");

    cpolpp->cd(1);
    tex = new TLatex( 0.08, 0.4, Form("Fit Py^{FP} = %3.2f +/- %3.2f", Pypp, ePypp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    tex = new TLatex( 0.08, 0.3, Form("Fit Px^{FP} = %3.2f +/- %3.2f", Pxpp, ePxpp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    tex = new TLatex( 0.08, 0.2, Form("A_{y}^{eff} = %3.2f +/- %3.2f", Aypp[1], eAypp[1]));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    Double_t FOMpp = effpp*Aypp[1]*Aypp[1];
    tex = new TLatex( 0.08, 0.1, Form("FOM = %3.2e", FOMpp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    Double_t PxpFP = hpol_SxpFP->GetMean();
    tex = new TLatex( 0.08, 0.0, Form("Chi = %3.1f degrees", Chipp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    cpolpp->Print("temp_polapp.pdf");
    
    //-----------------------------------------------------------------------------------------------------------------------------

    TCanvas* cpolnp = new TCanvas("cpolnp","",1200,800);
    cpolnp->Divide(4,2);

    cpolnp->cd(1);
    tex = new TLatex( 0.08, 0.85, "np -> pn ChEx (steel)");
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(1);
    tex->SetTextSize(0.075);
    tex->Draw();

    Double_t effnp = (Float_t)nnp_cxc/nQEn;
    tex = new TLatex( 0.08, 0.7, Form("efficiency = %3.2f %%", 100.*effnp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    cout << Form("np->pn (ChEx) efficiency = %3.2f %%", 100.*effnp) << endl;

    Double_t Ptn = hpol_Syn->GetMean();
    Double_t Pln = hpol_Szn->GetMean();
    tex = new TLatex( 0.08, 0.6, Form("Pt = %3.2f (#rightarrow Py^{FP})", Ptn)); 
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    tex = new TLatex( 0.08, 0.5, Form("Pl = %3.2f (#rightarrow Px^{FP})", Pln));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    cpolnp->cd(2)->SetLogy(1);
    hpol_sclnp->Draw("");
    hpol_sclnp->GetXaxis()->SetTitle("sclose [m]");

    cpolnp->cd(3);
    hpol_zclnp->Draw("colz");
    hpol_zclnp->GetXaxis()->SetTitle("zclose [m]");

    cpolnp->cd(4)->SetLogy(1);
    hpol_tnp_cxc->Draw("");
    hpol_tnp_cxc->GetXaxis()->SetTitle("|t| [GeV]");

    cpolnp->cd(5)->SetLogy(1);
    hpol_thnp_cx->Draw("");
    hpol_thnp_cx->GetXaxis()->SetTitle("#theta [deg]");

    hpol_thnp_cxc->Draw("SAME");
    hpol_thnp_cxc->SetLineColor(2);

    cpolnp->cd(6);
    hpol_phnp_cxm->Sumw2();
    hpol_phnp_cxm->Draw("");
    hpol_phnp_cxm->SetMarkerStyle(21);
    hpol_phnp_cxm->SetMarkerColor(1);
    hpol_phnp_cxm->SetLineColor(1);
    hpol_phnp_cxm->SetMarkerSize(0.5);
    hpol_phnp_cxm->GetXaxis()->SetTitle("#phi [deg]");

    hpol_phnp_cxp->Sumw2();
    hpol_phnp_cxp->SetMarkerStyle(21);
    hpol_phnp_cxp->SetMarkerColor(6);
    hpol_phnp_cxp->SetLineColor(6);
    hpol_phnp_cxp->SetMarkerSize(0.5);
    hpol_phnp_cxp->Draw("same");

    ofstream fOutFileTh, fOutFilePhP, fOutFilePhM;

    fOutFileTh.open("thpol.csv");
    fOutFilePhP.open("phpol_plus.csv");
    fOutFilePhM.open("phpol_minus.csv");
    
    for( int i=0; i < hpol_thnp_cxc->GetNbinsX(); i++ ) { 
      fOutFileTh << i << "," << hpol_thnp_cxc->GetBinCenter(i) << "," << hpol_thnp_cxc->GetBinContent(i) << endl; 
    }
    for( int i=0; i < hpol_phnp_cxp->GetNbinsX(); i++ ) { 
      fOutFilePhP << i << "," << hpol_phnp_cxp->GetBinCenter(i) << "," << hpol_phnp_cxp->GetBinContent(i) << endl; 
    }
    for( int i=0; i < hpol_phnp_cxm->GetNbinsX(); i++ ) {
      fOutFilePhM << i << "," << hpol_phnp_cxm->GetBinCenter(i) << "," << hpol_phnp_cxm->GetBinContent(i) << endl; 
    }
    
    fOutFileTh.close();
    fOutFilePhP.close();
    fOutFilePhM.close();

    cpolnp->cd(7);
    TH1F* hsumcx  = new TH1F("hsumcx","",nphibins,-180,180);
    TH1F* hdiffcx = new TH1F("hdiffcx","",nphibins,-180,180);
    TH1F* hasymcx = new TH1F("hasymcx","",nphibins,-180,180);
    hsumcx->Sumw2();
    hdiffcx->Sumw2();
    hasymcx->Sumw2();
    hsumcx->Add( hpol_phnp_cxp, hpol_phnp_cxm );
    hdiffcx->Add( hpol_phnp_cxp, hpol_phnp_cxm, 1, -1 );
    hasymcx->Divide( hdiffcx, hsumcx, 1, 1 );

//     TH1F* hasym1cx = new TH1F("hasym1cx","",nphibins,-180,180);
//     n2plus  = hpol_phnp_cxp->Integral(0,nphibins);
//     n2plus = nphibins/(n2plus);
//     n2minus = hpol_phnp_cxm->Integral(0,nphibins);
//     n2minus = nphibins/(n2minus);
//     hasym1cx->Add(hpol_phnp_cxp,hpol_phnp_cxm,n2plus,-n2minus);
    hasymcx->Draw("E1");
    hasymcx->SetMarkerStyle(21);
    hasymcx->SetMarkerColor(kBlue);
    hasymcx->SetLineColor(kBlue);
    hasymcx->SetMarkerSize(0.5);
    hasymcx->GetXaxis()->SetTitle("#varphi (degrees)");
    hasymcx->GetYaxis()->SetTitle("Asymmetry");
    hasymcx->GetYaxis()->SetRangeUser(-1.,1.);

    TF1 *fitcx = new TF1("fitcx", "[0]*sin(x/57.29578)-[1]*cos(x/57.29578)",-180,180);
    hasymcx->Fit(fitcx,"","",-180,180);

    Double_t Pynp  = fitcx->GetParameter(1);
    Double_t ePynp = fitcx->GetParError(1);
    Double_t Pxnp  = fitcx->GetParameter(0);
    Double_t ePxnp = fitcx->GetParError(0);

    Double_t Chinp = Chipp * 1.9/2.9;//85.6;

    cpolnp->cd(8);
    Double_t tnp[]   = { -1, hpol_tnp_cxc->GetMean(), 2 };
    Double_t etnp[]  = { 0, 0, 0 };
    
    Double_t Aynp[]  = { -3, fabs(Pxnp/Pln/sin(Chinp*M_PI/180)), -3 };
    Double_t eAynp[] = { 0, fabs(ePxnp/Pln/sin(Chinp*M_PI/180)), 0 };

    TGraphErrors* gAynp = new TGraphErrors( 3, tnp, Aynp, etnp, eAynp );
    gAynp->SetMarkerColor( 4 );
    gAynp->SetLineColor( 4 );
    gAynp->SetMarkerStyle( 20 );
    gAynp->SetMarkerSize( 1.5 );
    gAynp->Draw("AP");
    gAynp->GetXaxis()->SetRangeUser( 0, 0.15 );
    gAynp->GetYaxis()->SetRangeUser( 0, 1. );
    gAynp->GetXaxis()->SetTitle( "|t| [GeV]");
    gAynp->GetYaxis()->SetTitle( "Ay");

    cpolnp->cd(1);

    tex = new TLatex( 0.08, 0.4, Form("Fit Py^{FP} = %3.2f +/- %3.2f", Pynp, ePynp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    tex = new TLatex( 0.08, 0.3, Form("Fit Px^{FP} = %3.2f +/- %3.2f", Pxnp, ePxnp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    tex = new TLatex( 0.08, 0.2, Form("A_{y}^{eff} = %3.2f +/- %3.2f", Aynp[1], eAynp[1]));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    Double_t FOMnp = effnp*Aynp[1]*Aynp[1];
    tex = new TLatex( 0.08, 0.1, Form("FOM = %3.2e", FOMnp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();
    
    tex = new TLatex( 0.08, 0.0, Form("Chi = %3.1f degrees", Chinp));
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextColor(2);
    tex->SetTextSize(0.065);
    tex->Draw();

    cpolnp->Print("temp_polnp.pdf");


  }

  //-----------------------------------------------------------------------------------------------------------------------------
  
  gSystem->Exec("pdfunite  temp*.pdf Sim_GEnRP.pdf");  
  gSystem->Exec("rm  temp*.pdf");  
  
  if( writeOut ) {
    dir->GetList()->Write();   
    fout->Close(); 
  }
  
}

//-----------------------------------------------------------------------------------------------------------------------------
