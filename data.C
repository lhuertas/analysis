//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TAxis.h>
#include <TChainElement.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>
#include <TMath.h>
#include <TVectorT.h>
#include <THStack.h>
#include <TLatex.h>
#include <TRandom.h>
#include <TFractionFitter.h>
#include "Math/MinimizerOptions.h"
#include "TMinuit.h"


//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
using namespace std;

class BHRdmDataEvent {
   public:
      BHRdmDataEvent() {}
      ~BHRdmDataEvent() {}

      double xi_totem_right_rdm;
      double xi_totem_left_rdm;
      double t_totem_right_rdm;
      double t_totem_left_rdm;
      bool rp_right_rdm;
      bool rp_left_rdm;
      bool valid_vtx_rdm;
      bool valid_proton_right_rdm;
      bool valid_proton_left_rdm;
      double xi_cms_right;
      double xi_cms_left;
      double x_right;
      double x_left;
      double xi_totem_xicmscut_right;
      double t_totem_xicmscut_right;
      double xi_totem_xicmscut_left;
      double t_totem_xicmscut_left;
};      
///ROOUNFOLD
#include "/Users/lina/cernbox/doctorado/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/lina/cernbox/doctorado/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/lina/cernbox/doctorado/RooUnfold/src/RooUnfoldBinByBin.h"

// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldInvert.h"


double weight_right_bh, weight_right, weight_right_full, area_data_right;
double weight_left_bh, weight_left, weight_left_full, area_data_left;
int bmin_right, bmax_right, bmin_left, bmax_left;

double luminosity = 37.5;//nb^{-1} 0.91->merging efficiency rereco
//double luminosity = 38.53;//nb^{-1} 0.91->merging efficiency promptreco
TFile* data_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/data_ntuple.root","READ");
// TFile* data_file_rereco = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/rereco/data_ntuple_rereco.root","READ");
TFile* data_file_rereco = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_jecwinter.root","READ");
TString id;

double corr_t_vs_xi [8][5]; //(t, xi) 
void xi_corr(double t_proton, double xi_proton, double &corr){
    

    for (int i = 1; i <= 8; ++i){
        for (int j = 2; j <= 6; ++j){
                if (t_proton<tbins[i] && t_proton>tbins[i-1] && xi_proton<xi_bins[j] && xi_proton>xi_bins[j-1]) corr = corr_t_vs_xi[i][j];
                if (t_proton<tbins[i] && t_proton>tbins[i-1] && xi_proton>xi_bins[6]) corr = corr_t_vs_xi[i][6];
                if (t_proton>tbins[8] && xi_proton<xi_bins[j] && xi_proton>xi_bins[j-1]) corr = corr_t_vs_xi[8][j];
            }
        }  
}        


void data(TH1F* &xi_cms_minus_totem_right, TH1F* &xi_cms_minus_totem_right_bin, TH1F* &xi_cms_minus_totem_left, TH1F* &xi_cms_minus_totem_left_bin, TH1F* &xi_cms_minus_totem_right_bh, 
          TH1F* &xi_cms_minus_totem_left_bh, TH1F* &xi_right, TH1F* &xi_left, TH1F* &xi_right_cut, TH1F* &xi_right_bh, TH1F* &xi_right_bh_cut, TH1F* &xi_right_cut_noeff, TH1F* &xi_right_nojet,
          TH1F* &xi_left_cut, TH1F* &xi_left_bh, TH1F* &xi_left_bh_cut, TH1F* &xi_left_cut_noeff, TH1F* &t_right_cut, TH1F* &t_right_bh_cut, TH1F* &t_right_cut_zb, TH1F* &t_right_cut_full, 
          TH1F* &t_right_cut_bh, TH1F* &t_right_cut_noeff, TH1F* &t_right_nojet, TH1F* &t_left_cut, TH1F* &t_left_bh_cut, TH1F* &t_left_cut_zb, TH1F* &t_left_cut_full, TH1F* &t_left_cut_bh,
          TH1F* &t_left_cut_noeff, TH1F* &beta_right_cut, TH1F* &beta_left_cut, TH1F* &pt_jet1, TH1F* &pt_jet1_right_cut, TH1F* &pt_jet1_left_cut, TH1F* &pt_jet2, TH1F* &pt_jet2_right_cut, 
          TH1F* &pt_jet2_left_cut, TH1F* &eta_jet1, TH1F* &eta_jet1_right_cut, TH1F* &eta_jet1_left_cut, TH1F* &eta_jet2, TH1F* &eta_jet2_right_cut, TH1F* &eta_jet2_left_cut,
          TH1F* &delta_eta_jets_right_cut, TH1F* &delta_eta_jets_left_cut, TH1F* &delta_phi_jets_right_cut, TH1F* &delta_phi_jets_left_cut,
          TH1F* &log_x_right, TH1F* &log_x_right_noeff, TH1F* &log_x_right_cut,  TH1F* &log_x_right_bh_cut, TH1F* &log_x_right_cut_noeff, TH1F* &log_x_left, TH1F* &log_x_left_noeff, 
          TH1F* &log_x_left_cut, TH1F* &log_x_left_bh_cut, TH1F* &log_x_left_cut_noeff, 
          TH1F* &xi_right_cut_sasha, TH1F* &xi_left_cut_sasha, TH1F* &xi_cms_right_cut_sasha, TH1F* &xi_cms_left_cut_sasha, TH1F* &sigma_xi_cms_right_sasha, TH1F* &sigma_xi_cms_left_sasha, TH1F* &vtx_z, TH1F* &th_x_right, 
          TH1F* &th_y_right, TH1F* &th_x_left, TH1F* &th_y_left, TH1F* &mass_jj_right, TH1F* &mass_jj_left, TH1F* &mass_x_right, TH1F* &mass_x_left, TH1F* &r_jj_right, TH1F* &r_jj_left, 
          bool data_unc_rp = false, bool rereco=true, bool rp_xmin = false, bool rp_xmax = false, bool sys_xi_up = false, bool sys_xi_dw = false, bool single_vertex = false){

    id = "xx";
    if (data_unc_rp && !rp_xmin && !sys_xi_up && !sys_xi_dw) id = "unc_rp";
    else if (!data_unc_rp && rp_xmin && !sys_xi_up && !sys_xi_dw) id = "unc_rp_x";
    else if (!data_unc_rp && !rp_xmin && sys_xi_up && !sys_xi_dw) id = "unc_xi_up";
    else if (!data_unc_rp && !rp_xmin && !sys_xi_up && sys_xi_dw) id = "unc_xi_dw";
    else if (!data_unc_rp && !rp_xmin && !sys_xi_up && !sys_xi_dw) id = "nominal";
    cout<<id<<endl;

 
   // xi correction
    TFile* accep_t_vs_xi = TFile::Open("~/cernbox/doctorado/note/scripts/accep_t_vs_xi_nojet_right.root","READ");
    TH2F* accep_t_vs_xi_nojet_right = (TH2F*)accep_t_vs_xi->Get("t_rec_vs_xi_rec_minus_pomwig");
    for (int i = 1; i <= accep_t_vs_xi_nojet_right->GetNbinsX(); ++i)
    { 
        for (int j = 2; j <= 6; ++j)
        {
            cout<<i<<" "<<j<<" "<<accep_t_vs_xi_nojet_right->GetBinContent(i,j)<<endl;
            corr_t_vs_xi[i][j] = accep_t_vs_xi_nojet_right->GetBinContent(i,j);
        }   
    }
    accep_t_vs_xi->Close();

 

    xi_cms_minus_totem_right = new TH1F("xi_cms_minus_totem_right_"+id,"",50,-0.4,0.4);
    xi_cms_minus_totem_right_bin = new TH1F("xi_cms_minus_totem_right_bin_"+id,"", 15, bin);
    xi_cms_minus_totem_left = new TH1F("xi_cms_minus_totem_left_"+id,"",50,-0.4,0.4);
    xi_cms_minus_totem_left_bin = new TH1F("xi_cms_minus_totem_left_bin_"+id,"", 15, bin);
    xi_cms_minus_totem_right_bh = new TH1F("xi_cms_minus_totem_right_bh_"+id,"",50,-0.4,0.4);
    xi_cms_minus_totem_left_bh = new TH1F("xi_cms_minus_totem_left_bh_"+id,"",50,-0.4,0.4);
    xi_right = new TH1F("xi_right_"+id,"",50,-0.04,0.2);
    xi_left = new TH1F("xi_left_"+id,"",50,-0.04,0.2);
    xi_right_cut = new TH1F("xi_right_cut_"+id,"",11,xi_bins);
    xi_right_bh = new TH1F("xi_right_bh_"+id,"",50,-0.04,0.2);
    xi_right_bh_cut = new TH1F("xi_right_bh_cut_"+id,"",11,xi_bins);
    xi_right_cut_noeff = new TH1F("xi_right_cut_noeff_"+id,"",11,xi_bins);
    xi_right_nojet = new TH1F("xi_right_nojet_"+id,"",11,xi_bins);
    xi_left_cut = new TH1F("xi_left_cut_"+id,"",11,xi_bins);
    xi_left_bh = new TH1F("xi_left_bh_"+id,"",50,-0.04,0.2);
    xi_left_bh_cut = new TH1F("xi_left_bh_cut_"+id,"",11,xi_bins);
    xi_left_cut_noeff = new TH1F("xi_left_cut_noeff_"+id,"",11,xi_bins);
    t_right_cut = new TH1F("t_right_cut_"+id,"", 8, tbins);
    t_right_bh_cut = new TH1F("t_right_bh_cut_"+id,"", 8, tbins);
    t_right_cut_zb = new TH1F("t_right_cut_zb_"+id,"", 8, tbins);
    t_right_cut_full = new TH1F("t_right_cut_full_"+id,"", 8, tbins);
    t_right_cut_bh = new TH1F("t_right_cut_bh_"+id,"", 8, tbins);
    t_right_cut_noeff = new TH1F("t_right_cut_noeff_"+id,"", 8, tbins);
    t_right_nojet = new TH1F("t_right_nojet_"+id,"", 8, tbins);
    t_left_cut = new TH1F("t_left_cut_"+id,"", 8, tbins);
    t_left_bh_cut = new TH1F("t_left_bh_cut_"+id,"", 8, tbins);
    t_left_cut_zb = new TH1F("t_left_cut_full_"+id,"", 8, tbins);
    t_left_cut_full = new TH1F("t_left_cut_"+id,"", 8, tbins);
    t_left_cut_bh = new TH1F("t_left_cut_bh_"+id,"", 8, tbins);
    t_left_cut_noeff = new TH1F("t_left_cut_noeff_"+id,"", 8, tbins);
    beta_right_cut = new TH1F("beta_right_cut_"+id,"",15,0,1);
    beta_left_cut = new TH1F("beta_left_cut_"+id,"",15,0,1);
    pt_jet1 = new TH1F("pt_jet1_"+id,"", 15, 0, 200);
    pt_jet1_right_cut = new TH1F("pt_jet1_right_cut_"+id,"", 15, 0, 200);
    pt_jet1_left_cut = new TH1F("pt_jet1_left_cut_"+id,"", 15, 0, 200);
    pt_jet2 = new TH1F("pt_jet2_"+id,"", 15, 0, 200);
    pt_jet2_right_cut = new TH1F("pt_jet2_right_cut_"+id,"", 15, 0, 200);
    pt_jet2_left_cut = new TH1F("pt_jet2_left_cut_"+id,"", 15, 0, 200);
    eta_jet1 = new TH1F("eta_jet1_"+id,"", 20, -5.2, 5.2);
    eta_jet1_right_cut = new TH1F("eta_jet1_right_cut_"+id,"", 20, -5.2, 5.2);
    eta_jet1_left_cut = new TH1F("eta_jet1_left_cut_"+id,"", 20, -5.2, 5.2);
    eta_jet2 = new TH1F("eta_jet2_"+id,"", 20, -5.2, 5.2);
    eta_jet2_right_cut = new TH1F("eta_jet2_right_cut_"+id,"", 20, -5.2, 5.2);
    eta_jet2_left_cut = new TH1F("eta_jet2_left_cut_"+id,"", 20, -5.2, 5.2);
    delta_eta_jets_right_cut = new TH1F("delta_eta_jets_right_cut_"+id,"", 40, -5.2, 5.2);
    delta_phi_jets_right_cut = new TH1F("delta_phi_jets_right_cut_"+id,"", 40, -5.2, 5.2);
    delta_eta_jets_left_cut = new TH1F("delta_eta_jets_left_cut_"+id,"", 40, -5.2, 5.2);
    delta_phi_jets_left_cut = new TH1F("delta_phi_jets_left_cut_"+id,"", 40, -5.2, 5.2);
    log_x_right = new TH1F("log_x_right_"+id,"",15, -4, 0);
    log_x_right_noeff = new TH1F("log_x_right_noeff_"+id,"",15, -4, 0);
    log_x_right_cut = new TH1F("log_x_right_cut_"+id,"",15, -4, 0);
    log_x_right_bh_cut = new TH1F("log_x_right_bh_cut_"+id,"",15, -4, 0);
    log_x_right_cut_noeff = new TH1F("log_x_right_cut_noeff_"+id,"",15, -4, 0);
    log_x_left = new TH1F("log_x_left_"+id,"",15, -4, 0);
    log_x_left_noeff = new TH1F("log_x_left_noeff_"+id,"",15, -4, 0);
    log_x_left_cut = new TH1F("log_x_left_cut_"+id,"",15, -4, 0);
    log_x_left_bh_cut = new TH1F("log_x_left_bh_cut_"+id,"",15, -4, 0);
    log_x_left_cut_noeff = new TH1F("log_x_left_cut_noeff_"+id,"",15, -4, 0);
    xi_right_cut_sasha = new TH1F("xi_right_cut_sasha_"+id,"",11,xi_bins);
    xi_left_cut_sasha = new TH1F("xi_left_cut_sasha_"+id,"",8, bin_sasha);
    xi_cms_right_cut_sasha = new TH1F("xi_cms_right_sasha_data_"+id,"",8, bin_sasha);
    xi_cms_left_cut_sasha = new TH1F("xi_cms_left_sasha_data_"+id,"",8, bin_sasha);
    sigma_xi_cms_right_sasha = new TH1F("sigma_xi_cms_right_sasha_"+id,"",8, bin_sasha);
    sigma_xi_cms_left_sasha = new TH1F("sigma_xi_cms_left_sasha_"+id,"",8, bin_sasha);
    vtx_z = new TH1F("vtx_z_"+id,"",150,-30, 30);
    TH2F* proton_y_vs_x_rp_024_025_accept = new TH2F("proton_y_vs_x_rp_024_025_accept","proton_y_vs_x_rp_024_025",100,-10,10,100,-40,40);
    TH2F* proton_y_vs_x_rp_120_121_accept = new TH2F("proton_y_vs_x_rp_120_121_accept","proton_y_vs_x_rp_120_121",100,-10,10,100,-40,40);
    TH2F* proton_y_vs_x_rp_124_125_accept = new TH2F("proton_y_vs_x_rp_124_125_accept","proton_y_vs_x_rp_124_125",100,-10,10,100,-40,40);
    th_x_right = new TH1F("th_x_right_"+id,"",20, -0.4e-3, 0.4e-3);
    th_y_right = new TH1F("th_y_right_"+id,"",20, -0.4e-3, 0.4e-3);
    th_x_left = new TH1F("th_x_left_"+id,"",20, -0.4e-3, 0.4e-3);
    th_y_left = new TH1F("th_y_left_"+id,"",20, -0.4e-3, 0.4e-3);
    mass_x_right = new TH1F("mass_x_right","",20, 0, 1800);
    mass_x_left = new TH1F("mass_x_left","",20, 0, 1800);
    mass_jj_right = new TH1F("mass_jj_right","",20, 0, 1000);
    mass_jj_left = new TH1F("mass_jj_left","",20, 0, 1000);
    r_jj_right = new TH1F("r_jj_right","",20, 0, 1);
    r_jj_left = new TH1F("r_jj_left","",20, 0, 1);
     TH1F* x_pos_top_right = new TH1F("x_pos_top_right","", 20, -1, 10);
    TH1F* x_pos_bottom_right = new TH1F("x_pos_bottom_right","", 20, -1, 10);
    TH1F* y_pos_top_right = new TH1F("y_pos_top_right","", 20, 0, 35);
    TH1F* y_pos_bottom_right = new TH1F("y_pos_bottom_right","", 20, -35, 0);
    TH1F* x_pos_top_left = new TH1F("x_pos_top_left","", 20, -1, 10);
    TH1F* x_pos_bottom_left = new TH1F("x_pos_bottom_left","", 20, -1, 10);
    TH1F* y_pos_top_left = new TH1F("y_pos_top_left","", 20, 0, 35);
    TH1F* y_pos_bottom_left = new TH1F("y_pos_bottom_left","", 20, -35, 0);

    x_pos_top_right->Sumw2();
    y_pos_top_right->Sumw2();
    x_pos_bottom_right->Sumw2();
    y_pos_bottom_right->Sumw2();
    x_pos_top_left->Sumw2();
    y_pos_top_left->Sumw2();
    x_pos_bottom_left->Sumw2();
    y_pos_bottom_left->Sumw2();

   TH2F* xi_cms_vs_xi_proton = new TH2F("","", 11,xi_bins, 11,xi_bins);

    TH1F* eta_jet1_xi1_right = new TH1F("eta_jet1_xi1_right","", 20, -5.2, 5.2);
    TH1F* eta_jet1_xi2_right = new TH1F("eta_jet1_xi2_right","", 20, -5.2, 5.2);
    TH1F* eta_jet1_xi3_right = new TH1F("eta_jet1_xi3_right","", 20, -5.2, 5.2);
    TH1F* eta_jet1_xi4_right = new TH1F("eta_jet1_xi4_right","", 20, -5.2, 5.2);
    TH1F* eta_jet1_xi1_left = new TH1F("eta_jet1_xi1_right","", 20, -5.2, 5.2);
    TH1F* eta_jet1_xi2_left = new TH1F("eta_jet1_xi2_right","", 20, -5.2, 5.2);
    TH1F* eta_jet1_xi3_left = new TH1F("eta_jet1_xi3_right","", 20, -5.2, 5.2);
    TH1F* eta_jet1_xi4_left = new TH1F("eta_jet1_xi4_right","", 20, -5.2, 5.2);

    TTree* tree_data;
    if(rereco == false) tree_data = (TTree*) data_file->Get( treeName.c_str() );
    else tree_data = (TTree*) data_file_rereco->Get( treeName.c_str() );
    int nev_data = int(tree_data->GetEntriesFast());
    cout <<"The data file has " << nev_data << " entries : " << endl;
 
    double jet1_pt, jet1_eta, jet1_phi, jet2_pt, jet2_eta, jet2_phi, mjj2;
    double xi_cms_minus_data, xi_proton_right_data, xi_cms_plus_data, xi_proton_left_data, x_right_data, x_left_data;
    double t_proton_right_data, t_proton_left_data, beta_proton_right_data, beta_proton_left_data, eff_trigger;
    bool valid_proton_right_data, valid_proton_left_data;
    bool rp_right_top_data, rp_right_bottom_data, rp_left_top_data, rp_left_bottom_data; 
    bool primVtx = false;
    int nVtx;
    double x_pos_024, x_pos_025, x_pos_124, x_pos_125, y_pos_024, y_pos_025, y_pos_124, y_pos_125, vtx_x, vtx_y, vtx_z_pos, x_pos_120, x_pos_121, y_pos_120, y_pos_121, lumisec;
    double thx_proton_right, thx_proton_left, thy_proton_right, thy_proton_left;
    tree_data->SetBranchAddress("eff_trigger",&eff_trigger);
    tree_data->SetBranchAddress("xi_cms_minus",&xi_cms_minus_data);
    tree_data->SetBranchAddress("xi_cms_plus",&xi_cms_plus_data);
    tree_data->SetBranchAddress("xi_totem_right",&xi_proton_right_data);
    tree_data->SetBranchAddress("xi_totem_left",&xi_proton_left_data);
    tree_data->SetBranchAddress("t_totem_right",&t_proton_right_data);
    tree_data->SetBranchAddress("t_totem_left",&t_proton_left_data);
    tree_data->SetBranchAddress("beta_proton_right",&beta_proton_right_data);
    tree_data->SetBranchAddress("beta_proton_left",&beta_proton_left_data);
    tree_data->SetBranchAddress("x_right",&x_right_data);
    tree_data->SetBranchAddress("x_left",&x_left_data);
    tree_data->SetBranchAddress("valid_proton_right",&valid_proton_right_data);
    tree_data->SetBranchAddress("valid_proton_left",&valid_proton_left_data);
    tree_data->SetBranchAddress("rp_right_top",&rp_right_top_data);
    tree_data->SetBranchAddress("rp_right_bottom",&rp_right_bottom_data);
    tree_data->SetBranchAddress("rp_left_top",&rp_left_top_data);
    tree_data->SetBranchAddress("rp_left_bottom",&rp_left_bottom_data);
    tree_data->SetBranchAddress("jet1_pt",&jet1_pt);
    tree_data->SetBranchAddress("jet1_eta",&jet1_eta);
    tree_data->SetBranchAddress("jet1_phi",&jet1_phi);
    tree_data->SetBranchAddress("jet2_pt",&jet2_pt);
    tree_data->SetBranchAddress("jet2_eta",&jet2_eta);
    tree_data->SetBranchAddress("jet2_phi",&jet2_phi);
    tree_data->SetBranchAddress("x_pos_024",&x_pos_024);
    tree_data->SetBranchAddress("y_pos_024",&y_pos_024);
    tree_data->SetBranchAddress("x_pos_025",&x_pos_025);
    tree_data->SetBranchAddress("y_pos_025",&y_pos_025);
    tree_data->SetBranchAddress("x_pos_120",&x_pos_120);
    tree_data->SetBranchAddress("y_pos_120",&y_pos_120);
    tree_data->SetBranchAddress("x_pos_121",&x_pos_121);
    tree_data->SetBranchAddress("y_pos_121",&y_pos_121);
    tree_data->SetBranchAddress("x_pos_124",&x_pos_124);
    tree_data->SetBranchAddress("y_pos_124",&y_pos_124);
    tree_data->SetBranchAddress("x_pos_125",&x_pos_125);
    tree_data->SetBranchAddress("y_pos_125",&y_pos_125);
    tree_data->SetBranchAddress("vtx_x",&vtx_x);
    tree_data->SetBranchAddress("vtx_y",&vtx_y);
    tree_data->SetBranchAddress("vtx_z",&vtx_z_pos);
    if (rereco) tree_data->SetBranchAddress("select_Vertex",&primVtx);
    if (rereco) tree_data->SetBranchAddress("nVtx",&nVtx);
    if (rereco) tree_data->SetBranchAddress("thetax_proton_right",&thx_proton_right);
    if (rereco) tree_data->SetBranchAddress("thetay_proton_right",&thy_proton_right);
    if (rereco) tree_data->SetBranchAddress("thetax_proton_left",&thx_proton_left);
    if (rereco) tree_data->SetBranchAddress("thetay_proton_left",&thy_proton_left);
    if (rereco) tree_data->SetBranchAddress("mjj2",&mjj2);


 TRandom* rdm = new TRandom;
  rdm->SetSeed(12345);

    //beam halo
    std::vector<BHRdmDataEvent> beamhalo;
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev_data; ++i_evt){
        tree_data->GetEntry(i_evt);
        bool rp_right_data = rp_right_bottom_data || rp_right_top_data;
        bool rp_left_data = rp_left_bottom_data || rp_left_top_data;

        BHRdmDataEvent bh;
        bh.xi_cms_right = xi_cms_minus_data;
        bh.xi_cms_left = xi_cms_plus_data;
        if (xi_cms_minus_data>0.12 && valid_proton_right_data && rp_right_data && xi_proton_right_data>0 && xi_proton_right_data<0.1 && fabs(t_proton_right_data)>0.03 && fabs(t_proton_right_data)<1){
            bh.xi_totem_xicmscut_right = xi_proton_right_data;
            bh.t_totem_xicmscut_right = t_proton_right_data;
            bh.x_right = x_right_data;
        }    
        if (xi_cms_plus_data>0.12 && valid_proton_left_data && rp_left_data && xi_proton_left_data>0 && xi_proton_left_data<0.1 && t_proton_left_data>0.03 && t_proton_left_data<1){
            bh.xi_totem_xicmscut_left = xi_proton_left_data;
            bh.t_totem_xicmscut_left = t_proton_left_data;
            bh.x_left = x_left_data;
        }    
        beamhalo.push_back(bh);

    }//beam halo



    double eff_proton = 0.94;
    double eff_vtx;
      if (single_vertex) eff_vtx = 0.95;
      else eff_vtx = 1;

    int nevents_hera_right=0;
    int nevents_right=0;
    int nevents_hera_left=0;
    int nevents_left=0;
  
    for(int i_evt = 0; i_evt < nev_data; ++i_evt){
        tree_data->GetEntry(i_evt);
    // float tbins[9] = {0.03, 0.07, 0.12, 0.21, 0.31, 0.42,  0.55, 0.75, 1.};
        // if (!primVtx) continue;
        if (single_vertex && nVtx!=1) continue;
        if (!single_vertex && nVtx<1) continue;

        if (sys_xi_up) {
           xi_proton_right_data = xi_proton_right_data + xi_proton_right_data*0.1;
           xi_proton_left_data = xi_proton_left_data + xi_proton_left_data*0.1;
        }   
        if (sys_xi_dw){
           xi_proton_right_data = xi_proton_right_data - xi_proton_right_data*0.1;
           xi_proton_left_data = xi_proton_left_data - xi_proton_left_data*0.1;
        }
        
        bool jet_sel = jet1_pt>pt_threshold && jet2_pt>pt_threshold && fabs(jet1_eta)<4.4 && fabs(jet2_eta)<4.4;
        bool proton_right_kin_sel =  xi_proton_right_data>0 && xi_proton_right_data<0.1 && t_proton_right_data>0.03 && t_proton_right_data<1;
        bool proton_left_kin_sel = xi_proton_left_data>0 && xi_proton_left_data<0.1 && t_proton_left_data>0.03 && t_proton_left_data<1;

      	bool fid_cut_right_top = x_pos_124>0 && x_pos_124<7 && y_pos_124 >8.4 && y_pos_124<29;
	    bool fid_cut_right_bottom = x_pos_125>0 && x_pos_125<7 && y_pos_125<-8.4 && y_pos_125>-29;
        bool fid_cut_unc_right_top = x_pos_124>0 && x_pos_124<7 && y_pos_124 >8.2 && y_pos_124<27 ;
        bool fid_cut_unc_right_bottom = x_pos_125>0 && x_pos_125<7 && y_pos_125 <-8.2 && y_pos_125>-27 ;
        bool fid_cut_right_top_xmin = x_pos_124>0 && x_pos_124<6 && y_pos_124 >8.4 && y_pos_124<27;
        bool fid_cut_right_bottom_xmin = x_pos_125>0 && x_pos_125<6 && y_pos_125<-8.4 && y_pos_125>-27;
        bool fid_cut_right_top_xmax = x_pos_124>-0.5 && x_pos_124<8 && y_pos_124 >8.4 && y_pos_124<27;
        bool fid_cut_right_bottom_xmax = x_pos_125>-0.5 && x_pos_125<8 && y_pos_125<-8.4 && y_pos_125>-27;

        bool fid_cut_left_top = x_pos_024>0 && x_pos_024<7 && y_pos_024 >8.4 && y_pos_024<29;
        bool fid_cut_left_bottom = x_pos_025>0 && x_pos_025<7 && y_pos_025<-8.4 && y_pos_025>-29;
        bool fid_cut_unc_left_top = x_pos_024>0 && x_pos_024<7 && y_pos_024 >8.2 && y_pos_024<29 ;
        bool fid_cut_unc_left_bottom = x_pos_025>0 && x_pos_025<7 && y_pos_025 <-8.2 && y_pos_025>-29 ;
        bool fid_cut_left_top_xmin = x_pos_024>0 && x_pos_024<6 && y_pos_024 >8.4 && y_pos_024<29;
        bool fid_cut_left_bottom_xmin = x_pos_025>0 && x_pos_025<6 && y_pos_025<-8.4 && y_pos_025>-29;
        bool fid_cut_left_top_xmax = x_pos_024>-0.5 && x_pos_024<8 && y_pos_024 >8.4 && y_pos_024<29;
        bool fid_cut_left_bottom_xmax = x_pos_025>-0.5 && x_pos_025<8 && y_pos_025<-8.4 && y_pos_025>-29;

        bool rp_left_data_nom = (rp_left_top_data && fid_cut_left_top) || (rp_left_bottom_data && fid_cut_left_bottom);
        bool rp_left_data_unc = (rp_left_top_data && fid_cut_unc_left_top) || (rp_left_bottom_data && fid_cut_unc_left_bottom);
        // bool rp_left_data = (data_unc_rp==false) ? rp_left_data_nom : rp_left_data_unc;
        bool rp_left_data = false;
        if (data_unc_rp==false && rp_xmin == false && rp_xmax == false) rp_left_data = rp_left_data_nom;
        if (data_unc_rp==true && rp_xmin == false && rp_xmax == false) rp_left_data = rp_left_data_unc;
        if (data_unc_rp==false && rp_xmin == true && rp_xmax == false) rp_left_data = (rp_left_top_data && fid_cut_left_top_xmin) || (rp_left_bottom_data && fid_cut_left_bottom_xmin);
        if (data_unc_rp==false && rp_xmin == false && rp_xmax == true) rp_left_data = (rp_left_top_data && fid_cut_left_top_xmax) || (rp_left_bottom_data && fid_cut_left_bottom_xmax);

        bool rp_right_data_nom = (rp_right_top_data && fid_cut_right_top) || (rp_right_bottom_data && fid_cut_right_bottom);
        bool rp_right_data_unc = (rp_right_top_data && fid_cut_unc_right_top) || (rp_right_bottom_data && fid_cut_unc_right_bottom);
        // bool rp_right_data = (data_unc_rp==false) ? rp_right_data_nom : rp_right_data_unc;
        bool rp_right_data = false;
        if (data_unc_rp==false && rp_xmin == false && rp_xmax == false) rp_right_data = rp_right_data_nom;
        if (data_unc_rp==true && rp_xmin == false && rp_xmax == false) rp_right_data = rp_right_data_unc;
        if (data_unc_rp==false && rp_xmin == true && rp_xmax == false) rp_right_data = (rp_right_top_data && fid_cut_right_top_xmin) || (rp_right_bottom_data && fid_cut_right_bottom_xmin);
        if (data_unc_rp==false && rp_xmin == false && rp_xmax == true) rp_right_data = (rp_right_top_data && fid_cut_right_top_xmax) || (rp_right_bottom_data && fid_cut_right_bottom_xmax);

        int i_evt_backg_1 = 0 + rdm->Rndm()*(beamhalo.size());
        int i_evt_backg_2 = 0 + rdm->Rndm()*(beamhalo.size());
        BHRdmDataEvent const & beamhalo_data_1 = beamhalo.at(i_evt_backg_1);
        BHRdmDataEvent const & beamhalo_data_2 = beamhalo.at(i_evt_backg_2);
        double xi_right_cms_bh = beamhalo_data_1.xi_cms_right;
        double xi_right_totem_bh = beamhalo_data_2.xi_totem_xicmscut_right;
        double x_right_bh = beamhalo_data_2.x_right;
        double x_left_bh = beamhalo_data_2.x_left;
        double t_right_totem_bh = beamhalo_data_2.t_totem_xicmscut_right;
        double xi_left_cms_bh = beamhalo_data_1.xi_cms_left;
        double xi_left_totem_bh = beamhalo_data_2.xi_totem_xicmscut_left;
        double t_left_totem_bh = beamhalo_data_2.t_totem_xicmscut_left;
        vtx_z->Fill(vtx_z_pos, 1);
	
        double corr_xi_right = 1;//(xi_rec_proton_right_pom > 0.08) ? 1.32864 : 1.;
        xi_corr(t_proton_right_data, xi_proton_right_data, corr_xi_right);

        
       if (proton_right_kin_sel && valid_proton_right_data && rp_right_top_data && !rp_right_bottom_data){
           t_right_nojet->Fill(t_proton_right_data, 1);
           xi_right_nojet->Fill(xi_proton_right_data, 1);
        }
       if (proton_right_kin_sel && valid_proton_right_data && !rp_right_top_data && rp_right_bottom_data){
           t_right_nojet->Fill(t_proton_right_data, 1);
           xi_right_nojet->Fill(xi_proton_right_data, 1);
        }

        if (jet_sel){
            pt_jet1->Fill(jet1_pt, 1/(eff_trigger*eff_vtx));
            pt_jet2->Fill(jet2_pt, 1/(eff_trigger*eff_vtx));
            eta_jet1->Fill(jet1_eta, 1/(eff_trigger*eff_vtx));
            eta_jet2->Fill(jet2_eta, 1/(eff_trigger*eff_vtx));
            log_x_right->Fill(log10(x_right_data),1/(eff_trigger*eff_vtx));
            log_x_right_noeff->Fill(log10(x_right_data),1);
            log_x_left->Fill(log10(x_left_data),1/(eff_trigger*eff_vtx));
            log_x_left_noeff->Fill(log10(x_left_data),1);
            sigma_xi_cms_right_sasha->Fill(xi_cms_minus_data, 1/(eff_trigger*eff_vtx));
            sigma_xi_cms_left_sasha->Fill(xi_cms_plus_data, 1/(eff_trigger*eff_vtx));
        }
           
        // if (jet_sel && valid_proton_right_data && rp_right_data && xi_right_totem_bh<0.1 &&  t_right_totem_bh>0.03 && t_right_totem_bh<1 && xi_right_cms_bh - xi_right_totem_bh<0)t_right_bh_cut->Fill(t_right_totem_bh,0.8/(eff_trigger*eff_proton)); 
        if (jet_sel && proton_right_kin_sel && valid_proton_right_data && rp_right_data){
            xi_cms_minus_totem_right->Fill(xi_cms_minus_data - xi_proton_right_data, 1.);
            xi_cms_minus_totem_right_bin->Fill(xi_cms_minus_data - xi_proton_right_data, 1.);
            xi_cms_minus_totem_right_bh->Fill(xi_right_cms_bh - xi_right_totem_bh, 1.); 
            xi_right->Fill(xi_proton_right_data, 1./(eff_trigger*eff_proton*eff_vtx));
            xi_right_bh->Fill(xi_right_totem_bh, 1.); 
            if (xi_right_totem_bh>0 && xi_right_totem_bh<0.1 &&  t_right_totem_bh>0.03 && t_right_totem_bh<1 && xi_right_cms_bh - xi_right_totem_bh<0) t_right_cut_bh->Fill(t_right_totem_bh, 1);
            if (xi_cms_minus_data - xi_proton_right_data<0){++nevents_right;
                xi_right_cut->Fill(xi_proton_right_data, 1/(eff_trigger*eff_proton*eff_vtx));
                xi_cms_vs_xi_proton->Fill(xi_proton_right_data, xi_cms_minus_data, 1.);
                xi_right_cut_noeff->Fill(xi_proton_right_data, 1);
                xi_right_cut_sasha->Fill(xi_proton_right_data, 1/(eff_trigger*eff_proton*eff_vtx));
                t_right_cut->Fill(t_proton_right_data, 1/(eff_trigger*eff_proton*eff_vtx*corr_xi_right));
                t_right_cut_zb->Fill(t_proton_right_data, 1/(eff_trigger*eff_proton*eff_vtx));
                t_right_cut_full->Fill(t_proton_right_data, 1/(eff_trigger*eff_proton*eff_vtx));
                t_right_cut_noeff->Fill(t_proton_right_data, 1);
                beta_right_cut->Fill(beta_proton_right_data, 1/(eff_trigger*eff_proton*eff_vtx));
                log_x_right_cut->Fill(log10(x_right_data),1./(eff_trigger*eff_proton*eff_vtx));
                log_x_right_cut_noeff->Fill(log10(x_right_data),1.);
                pt_jet1_right_cut->Fill(jet1_pt, 1/(eff_trigger*eff_proton*eff_vtx));
                eta_jet1_right_cut->Fill(jet1_eta,1./(eff_trigger*eff_proton*eff_vtx));
                pt_jet2_right_cut->Fill(jet2_pt,1./(eff_trigger*eff_proton*eff_vtx));
                eta_jet2_right_cut->Fill(jet2_eta,1./(eff_trigger*eff_proton*eff_vtx));
                th_x_right->Fill(thx_proton_right, 1./(eff_trigger*eff_proton*eff_vtx));
                th_y_right->Fill(thy_proton_right, 1./(eff_trigger*eff_proton*eff_vtx));
                delta_eta_jets_right_cut->Fill(jet1_eta-jet2_eta, 1./(eff_trigger*eff_proton*eff_vtx)); 
                delta_phi_jets_right_cut->Fill(jet1_phi-jet2_phi, 1./(eff_trigger*eff_proton*eff_vtx));
                xi_cms_right_cut_sasha->Fill(xi_cms_minus_data, 1/(eff_trigger*eff_proton*eff_vtx));
                mass_jj_right->Fill(sqrt(mjj2), 1/(eff_trigger*eff_proton*eff_vtx));
                mass_x_right->Fill(4000*sqrt(xi_proton_right_data), 1/(eff_trigger*eff_proton*eff_vtx));
                r_jj_right->Fill(sqrt(mjj2)/(4000*sqrt(xi_proton_right_data)), 1/(eff_trigger*eff_proton*eff_vtx));

                if(xi_proton_right_data<0.048)eta_jet1_xi2_right->Fill(jet1_eta, 1/eff_trigger);
                if(xi_proton_right_data>0.048 && xi_proton_right_data<0.065) eta_jet1_xi3_right->Fill(jet1_eta, 1/eff_trigger);
                if(xi_proton_right_data>0.065 && xi_proton_right_data<0.1)eta_jet1_xi4_right->Fill(jet1_eta, 1/eff_trigger);

                if(xi_right_cms_bh - xi_right_totem_bh<0){
                    t_right_bh_cut->Fill(t_right_totem_bh,1); 
                    ++nevents_hera_right;
                    xi_right_bh_cut->Fill(xi_right_totem_bh, 1.); 
                    log_x_right_bh_cut->Fill(log10(x_right_bh), 1.); 
                }
            }

        }

        double corr_xi_left = 1;//(xi_rec_proton_right_pom > 0.08) ? 1.32864 : 1.;
        xi_corr(t_proton_left_data, xi_proton_left_data, corr_xi_left);

        if (jet_sel && proton_left_kin_sel && valid_proton_left_data && rp_left_data){
            xi_cms_minus_totem_left->Fill(xi_cms_plus_data - xi_proton_left_data, 1.); 
            xi_cms_minus_totem_left_bin->Fill(xi_cms_plus_data - xi_proton_left_data, 1.);
            xi_cms_minus_totem_left_bh->Fill(xi_left_cms_bh - xi_left_totem_bh, 1.); 
            xi_left->Fill(xi_proton_left_data, 1./(eff_trigger*eff_proton*eff_vtx));
            xi_left_bh->Fill(xi_left_totem_bh, 1.); 
            if (xi_left_totem_bh>0 && xi_left_totem_bh<0.1 &&  t_left_totem_bh>0.03 && t_left_totem_bh<1 && xi_left_cms_bh - xi_left_totem_bh<0) t_left_cut_bh->Fill(t_left_totem_bh, 1);
            if (xi_cms_plus_data - xi_proton_left_data<0){++nevents_left;
                xi_left_cut->Fill(xi_proton_left_data, 1/(eff_trigger*eff_proton*eff_vtx));
                xi_left_cut_noeff->Fill(xi_proton_left_data, 1);
                xi_left_cut_sasha->Fill(xi_proton_left_data, 1/(eff_trigger*eff_proton*eff_vtx));
                t_left_cut->Fill(t_proton_left_data, 1/(eff_trigger*eff_proton*eff_vtx*corr_xi_left));
                t_left_cut_zb->Fill(t_proton_left_data, 1/(eff_trigger*eff_proton*eff_vtx));
                t_left_cut_full->Fill(t_proton_left_data, 1/(eff_trigger*eff_proton*eff_vtx));
                t_left_cut_noeff->Fill(t_proton_left_data, 1);
                beta_left_cut->Fill(beta_proton_left_data, 1/(eff_trigger*eff_proton*eff_vtx));
                log_x_left_cut->Fill(log10(x_left_data),1./(eff_trigger*eff_proton*eff_vtx));
                log_x_left_cut_noeff->Fill(log10(x_left_data),1.);
                pt_jet1_left_cut->Fill(jet1_pt,1./(eff_trigger*eff_proton*eff_vtx));
                eta_jet1_left_cut->Fill(jet1_eta,1./(eff_trigger*eff_proton*eff_vtx));
                pt_jet2_left_cut->Fill(jet2_pt,1./(eff_trigger*eff_proton*eff_vtx));
                eta_jet2_left_cut->Fill(jet2_eta,1./(eff_trigger*eff_proton*eff_vtx));
                th_x_left->Fill(thx_proton_left, 1./(eff_trigger*eff_proton*eff_vtx)); 
                th_y_left->Fill(thy_proton_left, 1./(eff_trigger*eff_proton*eff_vtx)); 
                delta_eta_jets_left_cut->Fill(jet1_eta-jet2_eta, 1./(eff_trigger*eff_proton*eff_vtx)); 
                delta_phi_jets_left_cut->Fill(jet1_phi-jet2_phi, 1./(eff_trigger*eff_proton*eff_vtx));
                xi_cms_left_cut_sasha->Fill(xi_cms_plus_data, 1/(eff_trigger*eff_proton*eff_vtx));
                mass_jj_left->Fill(sqrt(mjj2), 1/(eff_trigger*eff_proton*eff_vtx));
                mass_x_left->Fill(4000*sqrt(xi_proton_left_data), 1/(eff_trigger*eff_proton*eff_vtx));
                r_jj_left->Fill(sqrt(mjj2)/(4000*sqrt(xi_proton_left_data)), 1/(eff_trigger*eff_proton*eff_vtx));

                if(xi_proton_left_data<0.048)eta_jet1_xi2_right->Fill(jet1_eta, 1/eff_trigger);
                if(xi_proton_left_data>0.048 && xi_proton_left_data<0.065) eta_jet1_xi3_right->Fill(jet1_eta, 1/eff_trigger);
                if(xi_proton_left_data>0.065 && xi_proton_left_data<0.1)eta_jet1_xi4_right->Fill(jet1_eta, 1/eff_trigger);

                if(xi_left_cms_bh - xi_left_totem_bh<0){
                  ++nevents_hera_left;
                  xi_left_bh_cut->Fill(xi_left_totem_bh, 1.); 
                  t_left_bh_cut->Fill(t_left_totem_bh, 1.); 
                    log_x_left_bh_cut->Fill(log10(x_left_bh), 1.); 
                }   
            }
        }
        //rp_pos
        if (jet_sel && proton_right_kin_sel && valid_proton_right_data && rp_right_top_data && !rp_right_bottom_data &&  fid_cut_right_top){
            proton_y_vs_x_rp_120_121_accept->Fill(x_pos_120, y_pos_120, 1.);
            proton_y_vs_x_rp_124_125_accept->Fill(x_pos_124, y_pos_124, 1.);
        }
        if (jet_sel && proton_right_kin_sel && valid_proton_right_data && !rp_right_top_data && rp_right_bottom_data && fid_cut_right_bottom){
            proton_y_vs_x_rp_120_121_accept->Fill(x_pos_121, y_pos_121, 1.);
            proton_y_vs_x_rp_124_125_accept->Fill(x_pos_125, y_pos_125, 1.);
        }
        if (jet_sel && proton_left_kin_sel && valid_proton_left_data && rp_left_top_data && !rp_left_bottom_data && fid_cut_left_top){
            proton_y_vs_x_rp_024_025_accept->Fill(x_pos_024, y_pos_024, 1.);
            // proton_y_vs_x_rp_024_025_accept->Fill(x_pos_025, y_pos_025, 1.);
          }
        if (jet_sel && proton_left_kin_sel && valid_proton_left_data && !rp_left_top_data && rp_left_bottom_data && fid_cut_left_bottom){
            // proton_y_vs_x_rp_024_025_accept->Fill(x_pos_024, y_pos_024, 1.);
            proton_y_vs_x_rp_024_025_accept->Fill(x_pos_025, y_pos_025, 1.);
        }

        if (valid_proton_right_data && proton_right_kin_sel){
            if (rp_right_top_data){
                x_pos_top_right->Fill(x_pos_124, 1);
                y_pos_top_right->Fill(y_pos_124, 1);
            }
            if (rp_right_bottom_data){
                x_pos_bottom_right->Fill(x_pos_125, 1);
                y_pos_bottom_right->Fill(y_pos_125, 1);
            }
        }
        if (valid_proton_left_data && proton_left_kin_sel){
            if (rp_left_top_data){
                x_pos_top_left->Fill(x_pos_024, 1);
                y_pos_top_left->Fill(y_pos_024, 1);
            }
            if (rp_left_bottom_data){
                x_pos_bottom_left->Fill(x_pos_025, 1);
                y_pos_bottom_left->Fill(y_pos_025, 1);
            }
        }

    

  
    }
    /*if (rereco) */cout<<"Rereco==    right: "<<nevents_right<<"  left: "<<nevents_left<<endl;

    // t_left_cut->SetXTitle("|t| (GeV^{2})");
    // t_left_cut->SetYTitle("dN/d|t|");
    // t_left_cut->Scale(1,"width");
    // t_left_cut->Draw("E1");
    // t_right_cut->Scale(1,"width");
    // t_right_cut->Draw("E1sames");
    // TLegend *leg1 = new TLegend(0.2,0.75,0.48,0.9);
    // leg1->AddEntry(t_left_cut,"Data - Left","lp");
    // leg1->AddEntry(t_right_cut,"Data - Right","lp");
    // leg1->SetFillColor(0);
    // leg1->SetLineColor(0);
    // leg1->Draw();     

    // proton_y_vs_x_rp_024_025_accept->Draw();
  // TAxis *axis_right = xi_cms_minus_totem_right->GetXaxis();
  // bmin_right = axis_right->FindBin(0.0); 
  // bmax_right = axis_right->FindBin(0.4); 
  // TAxis *axis_left = xi_cms_minus_totem_left->GetXaxis();
  // bmin_left = axis_left->FindBin(0.0); 
  // bmax_left = axis_left->FindBin(0.4); 

  // area_data_right = xi_cms_minus_totem_right->Integral(bmin_right,bmax_right);
  // double area_simulated_right_bh = xi_cms_minus_totem_right_bh->Integral(bmin_right,bmax_right);
  // // double area_simulated_right_full = xi_cms_minus_totem_backg_right_full->Integral(bmin_right,bmax_right);
  // weight_right_bh = area_data_right/area_simulated_right_bh;
  // // // double weight_right_full = area_data_right/area_simulated_right_full;
  // area_data_left = xi_cms_minus_totem_left->Integral(bmin_left,bmax_left);
  // double area_simulated_left_bh = xi_cms_minus_totem_left_bh->Integral(bmin_left,bmax_left);
  // // // double area_simulated_left_full = xi_cms_minus_totem_backg_left_full->Integral(bmin_left,bmax_left);
  // // // double weight_left_full = area_data_left/area_simulated_left_full;
  // weight_left_bh = area_data_left/area_simulated_left_bh;

// cout<<xi_right_cut->GetEntries()<<"  "<<xi_left_cut->GetEntries()<<endl;
  // if (rereco)  cout<<"ReReco=    hera right: "<<nevents_hera_right/t_right_bh_cut->GetEntries()<<"  hera left: "<<nevents_hera_left/t_left_bh_cut->GetEntries()<<endl;

//  xi_cms_vs_xi_proton->GetYaxis()->SetTitle("#xi_{CMS}");
// xi_cms_vs_xi_proton->GetXaxis()->SetTitle("#xi_{TOTEM}");
// xi_cms_vs_xi_proton->Draw(); 
//   t_left_cut->Scale(1,"width");
//   t_right_cut->Scale(1,"width");
//   t_left_cut->Fit("expo","I", "", 0.07, 0.45);
//   t_right_cut->Fit("expo","I", "", 0.07, 0.45);
//   t_left_cut->Draw("E1");
//   t_right_cut->Draw("E1same");

// eta_jet1_xi2_right->Draw("E1");
// eta_jet1_xi3_right->Draw("E1same");
// eta_jet1_xi4_right->Draw("E1same");

    // TFile * outfile = new TFile("data_histos_nojet.root", "RECREATE");
    // x_pos_top_right->Write();
    // y_pos_top_right->Write();
    // x_pos_bottom_right->Write();
    // y_pos_bottom_right->Write();
    // x_pos_top_left->Write();
    // y_pos_top_left->Write();
    // x_pos_bottom_left->Write();
    // y_pos_bottom_left->Write();
    // outfile->Close();

    // TFile * outfile = new TFile("data_histos_unc_up_y.root", "RECREATE");
    // t_right_cut->Write();
    // t_left_cut->Write();
    // xi_right_cut->Write();
    // xi_left_cut->Write();
    // outfile->Close();

//   TCanvas *c1_3 = new TCanvas("c1_3", "newpad");
// proton_y_vs_x_rp_124_125_accept->Draw();

//   TCanvas *c1_2 = new TCanvas("c1_2", "newpad");
// proton_y_vs_x_rp_024_025_accept->Draw();
    // mass_x_right->Draw();
    // mass_x_left->Draw("same");
}
