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
//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;

///ROOUNFOLD
#include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldBinByBin.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldInvert.h"
// Double_t beta_fit(Double_t *x, Double_t *par ){
//   Double_t result = 0;
//   // result = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0];
//   result = par[0]+par[1]*x[0]+ par[2]*pow(x[0],2) + par[3]*pow(x[0],3) +par[4]*pow(x[0],4);
//   return result;
// }
 

// TFile* pythia8 = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/pythia8_diff_ntuple.root","READ");
// TFile* pythia8 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia8_diff_ntuple_newtheta.root","READ");
TFile* pythia8 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pythia8_diff_ntuple_beam_smearing.root","READ");
// TFile* pythia8_CUETP8M1 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia8_diff_CUETP8M1_ntuple.root","READ");
TFile* pythia8_CUETP8M1 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pythia8_diff_CUETP8M1_ntuple_beam_smearing.root","READ");
// TFile* pythia8 = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pythia8_diff_test.root","READ");
TFile* mc_diff_file;
// double eff = 1.7e-4;//0.00133*survival_prob*1.22393;
double scale_rew_pythia_right;
double scale_rew_pythia_left;
double M_P = 0.938272046;
 /////////// reweight slope
double t_slope_right_pythia;
double t_slope_left_pythia; 
double scale_rew_slope_right;
double scale_rew_slope_left;

void pythia8_diff(TH1F* &xi_cms_minus_totem_rec_right, TH1F* &xi_cms_minus_totem_rec_left, TH1F* &xi_rec_right, TH1F* &xi_rec_left, TH1F* &xi_rec_right_cut, TH1F* &xi_gen_right_cut, TH1F* &xi_rec_left_cut, 
                  TH1F* &xi_gen_left_cut, TH1F* &t_rec_right_cut, TH1F* &t_gen_right_cut, TH1F* &t_rec_left_cut, TH1F* &t_gen_left_cut, TH1F* &beta_right_cut, TH1F* &beta_left_cut, TH1F* &log_x_rec_right_cut,
	                TH1F* &log_x_gen_right_cut, TH1F* &log_x_rec_left_cut, TH1F* &log_x_gen_left_cut, TH1F* &pt_jet1_right_cut, TH1F* &pt_jet1_left_cut, TH1F* &pt_jet2_right_cut, TH1F* &pt_jet2_left_cut, 
                  TH1F* &eta_jet1_right_cut, TH1F* &eta_jet1_left_cut, TH1F* &eta_jet2_right_cut, TH1F* &eta_jet2_left_cut,  
                  TH1F* &delta_eta_jets_right_cut, TH1F* &delta_eta_jets_left_cut, TH1F* &delta_phi_jets_right_cut, TH1F* &delta_phi_jets_left_cut,
                  TH1F* &th_x_rec_right, TH1F* &th_y_rec_right, TH1F* &th_x_rec_left, TH1F* &th_y_rec_left,
                  TH1F* &th_x_gen_right, TH1F* &th_y_gen_right, TH1F* &th_x_gen_left, TH1F* &th_y_gen_left,
                  TH1F* &mass_jj_rec_right, TH1F* &mass_jj_rec_left, TH1F* &mass_x_rec_right, TH1F* &mass_x_rec_left, TH1F* &r_jj_rec_right, TH1F* &r_jj_rec_left,
                  TH1F* &t_right_unfolded,  TH1F* &t_left_unfolded, TH1F* &xi_right_unfolded,  TH1F* &xi_left_unfolded,
                  TH1F* &logx_right_unfolded,  TH1F* &logx_left_unfolded, TH1F* &t_right_unfolded_noeff,  TH1F* &t_left_unfolded_noeff, TH1F* &xi_right_unfolded_noeff,  TH1F* &xi_left_unfolded_noeff,
                  TH1F* &logx_right_cut_unfolded_noeff,  TH1F* &logx_left_cut_unfolded_noeff, TH1F* &logx_right_unfolded_noeff,  TH1F* &logx_left_unfolded_noeff,  TH1F* &t_right_unfolded_binbybin,  
                  TH1F* &t_left_unfolded_binbybin, TH1F* &xi_right_unfolded_binbybin,  TH1F* &xi_left_unfolded_binbybin,
                  TH1F* &logx_right_unfolded_binbybin,  TH1F* &logx_left_unfolded_binbybin,  TH1F* &t_right_unfolded_up,  TH1F* &t_right_unfolded_dw, TH1F* &t_left_unfolded_up,  TH1F* &t_left_unfolded_dw, TH1F* &xi_right_unfolded_up,  TH1F* &xi_right_unfolded_dw,
                  TH1F* &xi_left_unfolded_up,  TH1F* &xi_left_unfolded_dw, TH1F* &log_x_right_unfolded_up,  TH1F* &log_x_right_unfolded_dw, TH1F* &log_x_left_unfolded_up,  TH1F* &log_x_left_unfolded_dw,  
                  TH1F* &t_right_unfolded_up_pf,  TH1F* &t_right_unfolded_dw_pf, TH1F* &t_left_unfolded_up_pf,  TH1F* &t_left_unfolded_dw_pf, TH1F* &xi_right_unfolded_up_pf,  TH1F* &xi_right_unfolded_dw_pf,
                  TH1F* &xi_left_unfolded_up_pf,  TH1F* &xi_left_unfolded_dw_pf, TH1F* &log_x_right_unfolded_up_pf,  TH1F* &log_x_right_unfolded_dw_pf, TH1F* &log_x_left_unfolded_up_pf,  TH1F* &log_x_left_unfolded_dw_pf, 
                  TH1F* &t_right_unfolded_herabackg, TH1F* &t_left_unfolded_herabackg, TH1F* &xi_right_unfolded_herabackg, TH1F* &xi_left_unfolded_herabackg, TH1F* &logx_right_unfolded_herabackg, TH1F* &logx_left_unfolded_herabackg,
                  TH1F* &t_right_pythia_unfolded, TH1F* &t_left_pythia_unfolded, TH1F* &xi_right_pythia_unfolded, TH1F* &xi_left_pythia_unfolded, TH1F* &logx_right_pythia_unfolded, TH1F* &logx_left_pythia_unfolded,
                  bool reweight = false, bool reweight_slope = false, bool unc_rp = false, bool unc_gauss = false, bool sys_xi_up = false, 
                  bool sys_xi_dw = false, bool rp_x_min = false, bool rp_x_max = false, bool single_vertex = false, bool iter_variation = false, 
                  string const& tune = "Tune4C"){

 id = "xx";
 if (!reweight && !reweight_slope && tune == "Tune4C") id = "norew_tune4C";
 if (reweight && !reweight_slope && tune == "Tune4C") id = "rew_tune4C";
 if (reweight && reweight_slope && tune == "Tune4C") id = "rew_slope_tune4C";
 if (!reweight && !reweight_slope && tune == "CUETP8M1") id = "norew_CUETP8M1";
 if (reweight && !reweight_slope && tune == "CUETP8M1") id = "rew_CUETP8M1";
 if (reweight && reweight_slope && tune == "CUETP8M1") id = "rew_slope_CUETP8M1";
 cout<<id<<endl;

 TF1* func = new TF1("func", beta_fit, 0., 1., 3);
 
 if (pt_threshold==35){
  func->SetParameter(0, 1.7614);
  func->SetParameter(1, -4.17946);
  func->SetParameter(2, 5.39813);
//   func_right->SetParameter(3, -107.509);
//   func_right->SetParameter(4, 68.6142);
 }
 if (pt_threshold==35){
  func->SetParameter(0, 1.57392);
  func->SetParameter(1, -4.14697);
  func->SetParameter(2, 4.78234);
//   func_left->SetParameter(3, -239.39);
//   func_left->SetParameter(4, 148.212);
 }

 if (pt_threshold==40){
  // func->SetParameter(0, 1.65312);
  // func->SetParameter(1, -4.34300);
  // func->SetParameter(2, 4.87932);
  //rereco

  if (tune == "Tune4C"){  
    func->SetParameter(0, 1.39814);//1.39861);//1.6199);
    func->SetParameter(1,  -3.28756);//-3.46778);//-4.38166);
    func->SetParameter(2, 3.39368);//3.70415);//5.45443);
    scale_rew_pythia_right = 1.0351;
    scale_rew_pythia_left = 1.0339;
    t_slope_right_pythia = 6.12568;
    t_slope_left_pythia = 6.12568; 
    scale_rew_slope_right = 0.8857;
    scale_rew_slope_left = 0.8855;
  } 
  if (tune == "CUETP8M1"){  
    func->SetParameter(0, 1.3652);
    func->SetParameter(1, -3.00332);
    func->SetParameter(2, 3.10908);
    scale_rew_pythia_right = 1.0333;
    scale_rew_pythia_left = 1.0353;
    t_slope_right_pythia = 6.04272;
    t_slope_left_pythia = 6.04272;
    scale_rew_slope_right = 0.8968;
    scale_rew_slope_left = 0.89731;
  }
 }

   //xi correction
    // TFile* accep_t_vs_xi = TFile::Open("~/cernbox/doctorado/note/scripts/accep_t_vs_xi_nojet_right.root","READ");
    // TH2F* accep_t_vs_xi_nojet_right = (TH2F*)accep_t_vs_xi->Get("t_rec_vs_xi_rec_minus_pomwig");
    // double corr_t_vs_xi [8][5]; //(t, xi) 
    // bool region_t_vs_xi [8][5];
    // for (int i = 1; i <= accep_t_vs_xi_nojet_right->GetNbinsX(); ++i)
    // { 
    //     for (int j = 2; j <= 6; ++j)
    //     {
    //         cout<<i<<" "<<j<<" "<<accep_t_vs_xi_nojet_right->GetBinContent(i,j)<<endl;
    //         corr_t_vs_xi[i][j] = accep_t_vs_xi_nojet_right->GetBinContent(i,j);
    //     }   
    // }

 ////////////// difff mc
     // float tbins[9] = {0.03, 0.08, 0.13, 0.21, 0.31, 0.41,  0.55, 0.75, 1.};
    float tbins[9] = {0.03, 0.07, 0.11, 0.21, 0.31, 0.45,  0.58, 0.72, 1.};
    float xi_bins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};
    float bin_sasha[9] ={0.0003, 0.002, 0.0045, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1};
     
    TH1F* xi_gen_right_cut_bin_sasha = new TH1F("xi_gen_right_cut_bin_sasha","", 8, bin_sasha);
    TH1F* xi_gen_right_cut_bin_sasha_pt40 = new TH1F("xi_gen_right_cut_bin_sasha_pt40","", 8, bin_sasha);
    TH1F* xi_gen_left_cut_bin_sasha = new TH1F("","", 8, bin_sasha);
    TH1F* xi_cms_gen_right_cut_bin_sasha = new TH1F("xi_cms_gen_right_cut_bin_sasha","", 8, bin_sasha);
    TH1F* xi_cms_gen_right_cut_bin_sasha_pt40 = new TH1F("xi_cms_gen_right_cut_bin_sasha_pt40","", 8, bin_sasha);
    TH1F* xi_cms_gen_left_cut_bin_sasha = new TH1F("","", 8, bin_sasha);
    xi_cms_minus_totem_rec_right = new TH1F("xi_cms_minus_totem_rec_right_pythia8diff_"+id,"",100,-0.4,0.4);
    xi_cms_minus_totem_rec_left = new TH1F("xi_cms_minus_totem_rec_left_pythia8diff_"+id,"",100,-0.4,0.4);
    xi_rec_right = new TH1F("xi_rec_right_pythia8diff_"+id,"",100,-0.04,0.2);
    xi_rec_left = new TH1F("xi_rec_left_pythia8diff_"+id,"",100,-0.04,0.2);
    xi_rec_right_cut = new TH1F("xi_rec_right_cut_pythia8diff_"+id,"",11,xi_bins);
    xi_gen_right_cut = new TH1F("xi_gen_right_cut_pythia8diff_"+id,"",11,xi_bins);
    xi_rec_left_cut = new TH1F("xi_rec_left_cut_pythia8diff_"+id,"",11,xi_bins);
    xi_gen_left_cut = new TH1F("xi_gen_left_cut_pythia8diff_"+id,"",11,xi_bins);
    t_rec_right_cut = new TH1F("t_rec_right_cut_pythia8diff_"+id,"", 8, tbins);
    t_gen_right_cut = new TH1F("t_gen_right_cut_pythia8diff_"+id,"", 8, tbins);
    t_rec_left_cut = new TH1F("t_rec_left_cut_pythia8diff_"+id,"", 8, tbins);
    t_gen_left_cut = new TH1F("t_gen_left_cut_pythia8diff_"+id,"", 8, tbins);
    beta_right_cut = new TH1F("beta_right_cut_pythia8diff_"+id,"",15,0,1);
    beta_left_cut = new TH1F("beta_left_cut_pythia8diff_"+id,"",15,0,1);
    log_x_rec_right_cut = new TH1F("log_x_rec_right_cut_pythia8diff_"+id,"",15, -4, 0);
    log_x_gen_right_cut = new TH1F("log_x_gen_right_cut_pythia8diff_"+id,"",15, -4, 0);
    log_x_rec_left_cut = new TH1F("log_x_rec_left_cut_pythia8diff_"+id,"",15, -4, 0);
    log_x_gen_left_cut = new TH1F("log_x_gen_left_cut_pythia8diff_"+id,"",15, -4, 0);
    pt_jet1_right_cut = new TH1F("pt_jet1_right_cut_pythia8diff_"+id,"", 15, 0, 200);
    pt_jet1_left_cut = new TH1F("pt_jet1_left_cut_pythia8diff_"+id,"", 15, 0, 200);
    pt_jet2_right_cut = new TH1F("pt_jet2_right_cut_pythia8diff_"+id,"", 15, 0, 200);
    pt_jet2_left_cut = new TH1F("pt_jet2_left_cut_pythia8diff_"+id,"", 15, 0, 200);
    eta_jet1_right_cut = new TH1F("eta_jet1_right_cut_pythia8diff_"+id,"", 20, -5.2, 5.2);
    eta_jet1_left_cut = new TH1F("eta_jet1_left_cut_pythia8diff_"+id,"", 20, -5.2, 5.2);
    eta_jet2_right_cut = new TH1F("eta_jet2_right_cut_pythia8diff_"+id,"", 20, -5.2, 5.2);
    eta_jet2_left_cut = new TH1F("eta_jet2_left_cut_pythia8diff_"+id,"", 20, -5.2, 5.2);
    TH2F* pt2_xi_gen_minus = new TH2F("pt2_xi_gen_minus_pythia8diff_"+id,"",20,0,1,20,0,0.1);
    TH2F* pt2_xi_gen_minus_rp = new TH2F("pt2_xi_gen_minus_rp_pythia8diff_"+id,"",20,0,1,20,0,0.1);
    TH2F* rp_pos_left = new TH2F("rp_pos_left","",100,-1,10,100,-40,40);
    TH2F* rp_pos_right = new TH2F("rp_pos_left","",100,-1,10,100,-40,40);
    TH1F* t_gen_left_bin = new TH1F("t_gen_left_cut_pythia8diff_"+id,"", 100,0,1);
    th_x_rec_right = new TH1F("th_x_rec_right_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    th_y_rec_right = new TH1F("th_y_rec_right_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    th_x_rec_left = new TH1F("th_x_rec_left_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    th_y_rec_left = new TH1F("th_y_rec_left_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    th_x_gen_right = new TH1F("th_x_gen_right_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    th_y_gen_right = new TH1F("th_y_gen_right_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    th_x_gen_left = new TH1F("th_x_gen_left_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    th_y_gen_left = new TH1F("th_y_gen_left_pythia8diff_"+id,"",20, -0.4e-3, 0.4e-3);
    delta_eta_jets_right_cut = new TH1F("delta_eta_jets_right_cut_pythia8diff_"+id,"", 40, -5.2, 5.2);
    delta_phi_jets_right_cut = new TH1F("delta_phi_jets_right_cut_pythia8diff_"+id,"", 40, -5.2, 5.2);
    delta_eta_jets_left_cut = new TH1F("delta_eta_jets_left_cut_pythia8diff_"+id,"", 40, -5.2, 5.2);
    delta_phi_jets_left_cut = new TH1F("delta_phi_jets_left_cut_pythia8diff_"+id,"", 40, -5.2, 5.2);
    TH2F* t_minus_th2 = new TH2F("t_minus_th2","", 8, tbins,8,tbins);
    TH2F* xi_minus_th2 = new TH2F("xi_minus_th2","",11,xi_bins,11,xi_bins);
    TH2F* logx_minus_th2 = new TH2F("logx_minus_th2","", 15, -4, 0,15, -4, 0);
    TH2F* t_plus_th2 = new TH2F("t_plus_th2","", 8, tbins,8,tbins);
    TH2F* xi_plus_th2 = new TH2F("xi_plus_th2","",11,xi_bins,11,xi_bins);
    TH2F* logx_plus_th2 = new TH2F("logx_plus_th2","", 15, -4, 0,15, -4, 0);
    TH2F* log_t_vs_log_xi_right = new TH2F("log_t_vs_log_xi_right","",100,-4,1,100,-4,-0.5);
    TH2F* log_t_vs_log_xi_left = new TH2F("log_t_vs_log_xi_left","",100,-4,1,100,-4,-0.5);
    mass_x_rec_right = new TH1F("mass_x_right","",20, 0, 1800);
    mass_x_rec_left = new TH1F("mass_x_left","",20, 0, 1800);
    mass_jj_rec_right = new TH1F("mass_jj_right","",20, 0, 1000);
    mass_jj_rec_left = new TH1F("mass_jj_left","",20, 0, 1000);
    r_jj_rec_right = new TH1F("r_jj_right","",20, 0, 1);
    r_jj_rec_left = new TH1F("r_jj_right","",20, 0, 1);
    TH1F* xi_rec_right_cut_nojet = new TH1F("xi_rec_right_cut_nojet","",11,xi_bins);
    TH1F* xi_gen_right_cut_nojet = new TH1F("xi_gen_right_cut_nojet","",11,xi_bins);
    TH1F* t_rec_right_cut_nojet = new TH1F("t_rec_right_cut_nojet","",8, tbins);
    TH1F* t_gen_right_cut_nojet = new TH1F("t_gen_right_cut_nojet","",8, tbins);
    TH2F* t_gen_vs_xi_gen_right_cut = new TH2F("t_gen_vs_xi_gen_right_cut", "", 8, tbins, 11, xi_bins);
    TH2F* t_rec_vs_xi_rec_right_cut = new TH2F("t_rec_vs_xi_rec_right_cut", "", 8, tbins, 11, xi_bins);
    TH1F* pt_jet1_gen_right_cut = new TH1F("", "", 20, 10, 200);
    TH2F* theta_y_vs_x_rec_nojet = new TH2F("theta_y_vs_x_rec_nojet","",20, -0.4e-3, 0.4e-3,20, -0.4e-3, 0.4e-3);
    TH2F* theta_y_vs_x_gen_nojet = new TH2F("theta_y_vs_x_gen_nojet","",20, -0.4e-3, 0.4e-3,20, -0.4e-3, 0.4e-3);
    TH1F* theta_x_rec_nojet = new TH1F("theta_x_rec_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F* theta_x_gen_nojet = new TH1F("theta_x_gen_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F* theta_y_rec_nojet = new TH1F("theta_y_rec_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F* theta_y_gen_nojet = new TH1F("theta_y_gen_nojet","",20, -0.4e-3, 0.4e-3);
    TH1F* x_pos_top_right = new TH1F("x_pos_top_right","", 20, -1, 10);
    TH1F* x_pos_bottom_right = new TH1F("x_pos_bottom_right","", 20, -1, 10);
    TH1F* y_pos_top_right = new TH1F("y_pos_top_right","", 20, 0, 35);
    TH1F* y_pos_bottom_right = new TH1F("y_pos_bottom_right","", 20, -35, 0);
    TH1F* x_pos_top_left = new TH1F("x_pos_top_left","", 20, -1, 10);
    TH1F* x_pos_bottom_left = new TH1F("x_pos_bottom_left","", 20, -1, 10);
    TH1F* y_pos_top_left = new TH1F("y_pos_top_left","", 20, 0, 35);
    TH1F* y_pos_bottom_left = new TH1F("y_pos_bottom_left","", 20, -35, 0);
  
    xi_rec_right_cut_nojet->Sumw2();
    xi_gen_right_cut_nojet->Sumw2();
    t_rec_right_cut_nojet->Sumw2();
    t_gen_right_cut_nojet->Sumw2();
    xi_cms_gen_right_cut_bin_sasha->Sumw2();
    xi_cms_gen_right_cut_bin_sasha_pt40->Sumw2();
    theta_x_rec_nojet->Sumw2();
    theta_y_rec_nojet->Sumw2();
    theta_x_gen_nojet->Sumw2();
    theta_y_gen_nojet->Sumw2();
    x_pos_top_right->Sumw2();
    y_pos_top_right->Sumw2();
    x_pos_bottom_right->Sumw2();
    y_pos_bottom_right->Sumw2();
    x_pos_top_left->Sumw2();
    y_pos_top_left->Sumw2();
    x_pos_bottom_left->Sumw2();
    y_pos_bottom_left->Sumw2();

    RooUnfoldResponse t_minus_response (t_rec_right_cut, t_gen_right_cut, "t_minus_unfolded", "t_minus_unfolded");
    RooUnfoldResponse t_plus_response (t_rec_left_cut, t_gen_left_cut, "t_plus_unfolded", "t_plus_unfolded");
    RooUnfoldResponse xi_minus_response (xi_rec_right_cut, xi_gen_right_cut, "xi_minus_unfolded", "xi_minus_unfolded");
    RooUnfoldResponse xi_plus_response (xi_rec_left_cut, xi_gen_left_cut, "xi_plus_unfolded", "xi_plus_unfolded");
    RooUnfoldResponse logx_minus_response (log_x_rec_right_cut, log_x_gen_right_cut, "logx_minus_unfolded", "logx_minus_unfolded");
    RooUnfoldResponse logx_plus_response (log_x_rec_left_cut, log_x_gen_left_cut, "logx_plus_unfolded", "logx_plus_unfolded");
    RooUnfoldResponse t_minus_response_back (t_gen_right_cut, t_rec_right_cut, "t_minus_unfolded", "t_minus_unfolded");
    RooUnfoldResponse xi_minus_response_back (xi_gen_right_cut, xi_rec_right_cut, "xi_minus_unfolded", "xi_minus_unfolded");
    RooUnfoldResponse logx_minus_response_back (log_x_gen_right_cut, log_x_rec_right_cut, "logx_minus_unfolded", "logx_minus_unfolded");

    double norm;
    if (tune == "Tune4C"){ norm = 2109.82*1.7e-4; mc_diff_file = pythia8;}
    if (tune == "CUETP8M1"){ norm = 15790.5*4.3e-5; mc_diff_file = pythia8_CUETP8M1;}
    string treeName = "small_tree"; 
    TTree* tree= (TTree*) mc_diff_file->Get( treeName.c_str() );
    int nev = int(tree->GetEntriesFast());
    cout <<"The diffractive file has " << nev << " entries : " << endl;
 
    double jet1_rec_pt, jet1_rec_eta, jet1_rec_phi, jet2_rec_pt, jet2_rec_eta, jet2_rec_phi, mjj2_rec;
    double jet1_gen_pt, jet1_gen_eta, jet1_gen_phi, jet2_gen_pt, jet2_gen_eta, jet2_gen_phi;
    double xi_rec_cms_minus, xi_rec_proton_right, xi_rec_cms_plus, xi_rec_proton_left;
    double xi_gen_cms_minus, xi_gen_proton_right, xi_gen_cms_plus, xi_gen_proton_left;
    double t_rec_proton_right, beta_rec_proton_right, t_rec_proton_left, beta_rec_proton_left;
    double t_gen_proton_right, beta_gen_proton_right, t_gen_proton_left, beta_gen_proton_left;
    double weight, x_rec_left, x_rec_right, x_gen_left, x_gen_right;
    bool rp_right, rp_left, rp_right_accep_bottom, rp_right_accep_top, rp_left_accep_bottom, rp_left_accep_top;
    double rp_xpos_125, rp_xpos_124, rp_ypos_125, rp_ypos_124, rp_xpos_24, rp_xpos_25, rp_ypos_24, rp_ypos_25;
    double px_plus, py_plus, pz_plus, xi_rec_proton_right_gauss, xi_rec_proton_left_gauss, t_rec_proton_right_gauss, t_rec_proton_left_gauss;
    int nVtx;
    double theta_x_minus_smear, theta_y_minus_smear, theta_x_plus_smear, theta_y_plus_smear, pz_minus_smear, pz_plus_smear, e_minus_smear, e_plus, e_plus_smear, pz_minus, px_minus, py_minus, e_minus, mass_plus, mass_minus;
    double theta_x_minus, theta_y_minus, theta_x_plus, theta_y_plus;
    tree->SetBranchAddress("xi_rec_cms_right",&xi_rec_cms_minus);
    tree->SetBranchAddress("xi_rec_cms_left",&xi_rec_cms_plus);
    tree->SetBranchAddress("xi_gen_cms_right",&xi_gen_cms_minus);
    tree->SetBranchAddress("xi_gen_cms_left",&xi_gen_cms_plus);
    tree->SetBranchAddress("xi_rec_proton_right",&xi_rec_proton_right);
    tree->SetBranchAddress("xi_rec_proton_left",&xi_rec_proton_left);
    tree->SetBranchAddress("xi_gen_proton_right",&xi_gen_proton_right);
    tree->SetBranchAddress("xi_gen_proton_left",&xi_gen_proton_left);
    tree->SetBranchAddress("t_rec_proton_right",&t_rec_proton_right);
    tree->SetBranchAddress("t_rec_proton_left",&t_rec_proton_left);
    tree->SetBranchAddress("t_gen_proton_right",&t_gen_proton_right);
    tree->SetBranchAddress("t_gen_proton_left",&t_gen_proton_left);
    tree->SetBranchAddress("beta_rec_right",&beta_rec_proton_right);
    tree->SetBranchAddress("beta_rec_left",&beta_rec_proton_left);
    tree->SetBranchAddress("beta_gen_right",&beta_gen_proton_right);
    tree->SetBranchAddress("beta_gen_left",&beta_gen_proton_left);
    tree->SetBranchAddress("rp_right",&rp_right);
    tree->SetBranchAddress("rp_left",&rp_left);
    tree->SetBranchAddress("jet1_rec_pt",&jet1_rec_pt);
    tree->SetBranchAddress("jet1_gen_pt",&jet1_gen_pt);
    tree->SetBranchAddress("jet1_rec_eta",&jet1_rec_eta);
    tree->SetBranchAddress("jet1_gen_eta",&jet1_gen_eta);
    tree->SetBranchAddress("jet1_rec_phi",&jet1_rec_phi);
    tree->SetBranchAddress("jet2_rec_pt",&jet2_rec_pt);
    tree->SetBranchAddress("jet2_gen_pt",&jet2_gen_pt);
    tree->SetBranchAddress("jet2_rec_eta",&jet2_rec_eta);
    tree->SetBranchAddress("jet2_gen_eta",&jet2_gen_eta);
    tree->SetBranchAddress("jet2_rec_phi",&jet2_rec_phi);
    tree->SetBranchAddress("x_rec_right",&x_rec_right);
    tree->SetBranchAddress("x_gen_right",&x_gen_right);
    tree->SetBranchAddress("x_rec_left",&x_rec_left);
    tree->SetBranchAddress("x_gen_left",&x_gen_left);
    tree->SetBranchAddress("weight_mc",&weight);
    tree->SetBranchAddress("rp_right_accep_bottom",&rp_right_accep_bottom);
    tree->SetBranchAddress("rp_right_accep_top",&rp_right_accep_top);
    tree->SetBranchAddress("rp_xpos_124",&rp_xpos_124);
    tree->SetBranchAddress("rp_ypos_124",&rp_ypos_124);
    tree->SetBranchAddress("rp_xpos_125",&rp_xpos_125);
    tree->SetBranchAddress("rp_ypos_125",&rp_ypos_125);
    tree->SetBranchAddress("rp_left_accep_bottom",&rp_left_accep_bottom);
    tree->SetBranchAddress("rp_left_accep_top",&rp_left_accep_top);
    tree->SetBranchAddress("rp_xpos_24",&rp_xpos_24);
    tree->SetBranchAddress("rp_ypos_24",&rp_ypos_24);
    tree->SetBranchAddress("rp_xpos_25",&rp_xpos_25);
    tree->SetBranchAddress("rp_ypos_25",&rp_ypos_25);
    tree->SetBranchAddress("px_proton_right",&px_minus);
    tree->SetBranchAddress("py_proton_right",&py_minus);
    tree->SetBranchAddress("pz_proton_right",&pz_minus);
    tree->SetBranchAddress("px_proton_left",&px_plus);
    tree->SetBranchAddress("py_proton_left",&py_plus);
    tree->SetBranchAddress("pz_proton_left",&pz_plus);
    // tree->SetBranchAddress("pz_proton_right_smear",&pz_minus_smear);
    tree->SetBranchAddress("e_proton_right",&e_minus);
    tree->SetBranchAddress("e_proton_left",&e_plus);
    // tree->SetBranchAddress("e_proton_right_smear",&e_minus_smear);
    tree->SetBranchAddress("theta_x_plus_smear",&theta_x_plus_smear);
    tree->SetBranchAddress("theta_y_plus_smear",&theta_y_plus_smear);
    tree->SetBranchAddress("theta_x_minus_smear",&theta_x_minus_smear);
    tree->SetBranchAddress("theta_y_minus_smear",&theta_y_minus_smear);
    tree->SetBranchAddress("theta_x_plus",&theta_x_plus);
    tree->SetBranchAddress("theta_y_plus",&theta_y_plus);
    tree->SetBranchAddress("theta_x_minus",&theta_x_minus);
    tree->SetBranchAddress("theta_y_minus",&theta_y_minus);
    tree->SetBranchAddress("mass_proton_left",&mass_plus);
    tree->SetBranchAddress("mass_proton_right",&mass_minus);
    tree->SetBranchAddress("nVtx",&nVtx);
    tree->SetBranchAddress("xi_rec_proton_right_gauss",&xi_rec_proton_right_gauss);
    tree->SetBranchAddress("t_rec_proton_right_gauss",&t_rec_proton_right_gauss);
    tree->SetBranchAddress("xi_rec_proton_left_gauss",&xi_rec_proton_left_gauss);
    tree->SetBranchAddress("t_rec_proton_left_gauss",&t_rec_proton_left_gauss);
    tree->SetBranchAddress("mjj2_rec",&mjj2_rec);

     TH1F* e_minus_histo = new TH1F("","",100,0,5000);
    double pi = 4000;

    for(int i_evt = 0; i_evt < nev; ++i_evt){
        tree->GetEntry(i_evt);

        if (single_vertex && nVtx!=1) continue;
        if (!single_vertex && nVtx<1) continue;

        xi_rec_proton_right = (unc_gauss == false) ? xi_rec_proton_right : xi_rec_proton_right_gauss;
        t_rec_proton_right = (unc_gauss == false) ? t_rec_proton_right : t_rec_proton_right_gauss; 
        xi_rec_proton_left = (unc_gauss == false) ? xi_rec_proton_left : xi_rec_proton_left_gauss;
        t_rec_proton_left = (unc_gauss == false) ? t_rec_proton_left : t_rec_proton_left_gauss; 

        double reweigth_beta_minus = (beta_gen_proton_right<=0.7) ? func->Eval(beta_gen_proton_right) : 1;
        double reweigth_beta_plus = (beta_gen_proton_left<=0.7) ? func->Eval(beta_gen_proton_left) : 1;
	      double reweight_slope_minus = (fabs(t_gen_proton_right)<=0.45) ? scale_rew_slope_right*exp(-t_slope_right_data*fabs(t_gen_proton_right))/exp(-t_slope_right_pythia*fabs(t_gen_proton_right)) : 1;
	      double reweight_slope_plus = (fabs(t_gen_proton_left)<=0.45) ? scale_rew_slope_left*exp(-t_slope_left_data*fabs(t_gen_proton_left))/exp(-t_slope_left_pythia*fabs(t_gen_proton_left)) : 1;
        double event_weight_minus, event_weight_plus;
	
	    if (reweight && !reweight_slope){
	        event_weight_minus = scale_rew_pythia_right*reweigth_beta_minus*weight;
	        event_weight_plus = scale_rew_pythia_left*reweigth_beta_plus*weight;}
	    if (reweight && reweight_slope){
	       event_weight_minus = scale_rew_pythia_right*reweigth_beta_minus*weight*reweight_slope_minus;
	       event_weight_plus = scale_rew_pythia_left*reweigth_beta_plus*weight*reweight_slope_plus;}
	    if (!reweight && !reweight_slope){
	       event_weight_minus = weight;
	       event_weight_plus = weight;}

        // double px_minus_smear = -pi*tan(theta_x_minus_smear);
        // double py_minus_smear = pi*tan(theta_y_minus_smear);
        // TLorentzVector p_beam_minus (0, 0, -pi, pi);
        // TLorentzVector p_scatt_minus (px_minus_smear, py_minus_smear, pz_minus_smear, e_minus_smear);
        // TLorentzVector t_vec_minus = (p_beam_minus - p_scatt_minus);
        // double t_rec_proton_right_theta = t_vec_minus.Mag2(); 

        // double px_plus_smear = -pi*tan(theta_x_plus_smear);
        // double py_plus_smear = pi*tan(theta_y_plus_smear);
        // double pi_mass = sqrt(pi*pi - M_P*M_P);
        // TLorentzVector p_beam_plus (0, 0, pi_mass, pi);
        // double e_gen_plus_from_p = sqrt(px_plus*px_plus + py_plus*py_plus + pz_plus*pz_plus + M_P*M_P);
        // TLorentzVector p_scatt_plus (px_plus, py_plus, pz_plus, e_gen_plus_from_p);
        // TLorentzVector t_vec_plus = (p_beam_plus - p_scatt_plus);
        // double t_rec_proton_left_theta = t_vec_plus.Mag2(); 

        bool jet_rec_sel = jet1_rec_pt>pt_threshold && jet2_rec_pt>pt_threshold && fabs(jet1_rec_eta)<4.4 && fabs(jet2_rec_eta)<4.4;
        bool jet_gen_sel = jet1_gen_pt>pt_threshold && jet2_gen_pt>pt_threshold && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4;
        bool proton_right_rec_sel =  xi_rec_proton_right>0 && xi_rec_proton_right<0.1 && fabs(t_rec_proton_right)>0.03 && fabs(t_rec_proton_right)<1;
        bool proton_right_gen_sel = xi_gen_proton_right>0 && xi_gen_proton_right<0.1 && fabs(t_gen_proton_right)>0.03 && fabs(t_gen_proton_right)<1;
        bool proton_left_rec_sel = xi_rec_proton_left>0 && xi_rec_proton_left<0.1 && fabs(t_rec_proton_left)>0.03 && fabs(t_rec_proton_left)<1;
        bool proton_left_gen_sel =  xi_gen_proton_left>0 && xi_gen_proton_left<0.1 && fabs(t_gen_proton_left)>0.03 && fabs(t_gen_proton_left)<1;

        bool fid_cuts_nom_right_top = rp_xpos_124>0 && rp_xpos_124<0.007 && rp_ypos_124 >0.0084 && rp_ypos_124<0.029 ;
        bool fid_cuts_nom_right_bottom = rp_xpos_125>0 && rp_xpos_125<0.007 && rp_ypos_125 <-0.0084 && rp_ypos_125>-0.029 ;
        bool fid_cuts_unc_right_top = rp_xpos_124>0 && rp_xpos_124<0.007 && rp_ypos_124 >0.0082 && rp_ypos_124<0.027 ;
        bool fid_cuts_unc_right_bottom = rp_xpos_125>0 && rp_xpos_125<0.007 && rp_ypos_125 <-0.0082 && rp_ypos_125>-0.027 ;
        bool rp_right_unc = (rp_right_accep_top && fid_cuts_unc_right_top) || (rp_right_accep_bottom && fid_cuts_unc_right_bottom);

        bool fid_cuts_unc_right_top_xmax = rp_xpos_124>-0.0005 && rp_xpos_124<0.008 && rp_ypos_124 >0.0084 && rp_ypos_124<0.027 ;
        bool fid_cuts_unc_right_bottom_xmax = rp_xpos_125>-0.0005 && rp_xpos_125<0.008 && rp_ypos_125 <-0.0084 && rp_ypos_125>-0.027 ;
        bool fid_cuts_unc_right_top_xmin = rp_xpos_124>0.00 && rp_xpos_124<0.006 && rp_ypos_124 >0.0084 && rp_ypos_124<0.027 ;
        bool fid_cuts_unc_right_bottom_xmin = rp_xpos_125>0.00 && rp_xpos_125<0.006 && rp_ypos_125 <-0.0084 && rp_ypos_125>-0.027 ;

        bool rp_right_sel = false; 
        if(unc_rp == false && rp_x_min == false && rp_x_max == false) rp_right_sel = (rp_right_accep_top && fid_cuts_nom_right_top) || (rp_right_accep_bottom && fid_cuts_nom_right_bottom);
        if(unc_rp == true && rp_x_min == false && rp_x_max == false) rp_right_sel = rp_right_unc ;
        if(unc_rp == false && rp_x_min == true && rp_x_max == false) rp_right_sel = (rp_right_accep_top && fid_cuts_unc_right_top_xmin) || (rp_right_accep_bottom && fid_cuts_unc_right_bottom_xmin);
        if(unc_rp == false && rp_x_min == false && rp_x_max == true) rp_right_sel = (rp_right_accep_top && fid_cuts_unc_right_top_xmax) || (rp_right_accep_bottom && fid_cuts_unc_right_bottom_xmax);

        bool fid_cuts_nom_left_top = rp_xpos_24>0 && rp_xpos_24<0.007 && rp_ypos_24 >0.0084 && rp_ypos_24<0.029 ;
        bool fid_cuts_nom_left_bottom = rp_xpos_25>0 && rp_xpos_25<0.007 && rp_ypos_25 <-0.0084 && rp_ypos_25>-0.029 ;
        bool fid_cuts_unc_left_top = rp_xpos_24>0 && rp_xpos_24<0.007 && rp_ypos_24 >0.0082 && rp_ypos_24<0.029 ;
        bool fid_cuts_unc_left_bottom = rp_xpos_25>0 && rp_xpos_25<0.007 && rp_ypos_25 <-0.0082 && rp_ypos_25>-0.029 ;
        bool fid_cuts_unc_left_top_xmin = rp_xpos_24>0.000 && rp_xpos_24<0.006 && rp_ypos_24 >0.0084 && rp_ypos_24<0.029 ;
        bool fid_cuts_unc_left_top_xmax = rp_xpos_24>-0.0005 && rp_xpos_24<0.008 && rp_ypos_24 >0.0084 && rp_ypos_24<0.029 ;
        bool fid_cuts_unc_left_bottom_xmin = rp_xpos_25>0.000 && rp_xpos_25<0.006 && rp_ypos_25 <-0.0084 && rp_ypos_25>-0.029 ;
        bool fid_cuts_unc_left_bottom_xmax = rp_xpos_25>-0.0005 && rp_xpos_25<0.008 && rp_ypos_25 <-0.0084 && rp_ypos_25>-0.029 ;

        bool rp_left_unc = (rp_left_accep_top && fid_cuts_unc_left_top) || (rp_left_accep_bottom && fid_cuts_unc_left_bottom);
        // bool rp_left_sel = (unc_rp == false) ? rp_left : rp_left_unc;
        bool rp_left_sel = false; 
        if(unc_rp == false && rp_x_min == false && rp_x_max == false) rp_left_sel = (rp_left_accep_top && fid_cuts_nom_left_top) || (rp_left_accep_bottom && fid_cuts_nom_left_bottom);
        if(unc_rp == true && rp_x_min == false && rp_x_max == false) rp_left_sel = rp_left_unc ;
        if(unc_rp == false && rp_x_min == true && rp_x_max == false) rp_left_sel = (rp_left_accep_top && fid_cuts_unc_left_top_xmin) || (rp_left_accep_bottom && fid_cuts_unc_left_bottom_xmin);
        if(unc_rp == false && rp_x_min == false && rp_x_max == true) rp_left_sel = (rp_left_accep_top && fid_cuts_unc_left_top_xmax) || (rp_left_accep_bottom && fid_cuts_unc_left_bottom_xmax);

        double corr_xi = 1;//(xi_rec_proton_right_pom > 0.08) ? 1.32864 : 1.;
        // for (int i = 1; i <= 8; ++i){
        //     for (int j = 2; j <= 6; ++j){
        //         if (fabs(t_rec_proton_right)<tbins[i] && fabs(t_rec_proton_right)>tbins[i-1] && xi_rec_proton_right<xi_bins[j]
        //             && xi_rec_proton_right>xi_bins[j-1]) corr_xi = corr_t_vs_xi[i][j];
        //         if (fabs(t_rec_proton_right)<tbins[i] && fabs(t_rec_proton_right)>tbins[i-1] && xi_rec_proton_right>xi_bins[6]) 
        //             corr_xi = corr_t_vs_xi[i][6];
        //         if (fabs(t_rec_proton_right)>tbins[8] && xi_rec_proton_right<xi_bins[j] && xi_rec_proton_right>xi_bins[j-1])
        //             corr_xi = corr_t_vs_xi[8][j];
        //     }
        // }        
           
        log_t_vs_log_xi_right->Fill(log10(fabs(t_gen_proton_right)), log10(xi_gen_proton_right));
        log_t_vs_log_xi_left->Fill(log10(fabs(t_gen_proton_left)), log10(xi_gen_proton_left));

        if (proton_right_rec_sel){
           if (rp_right_accep_top){
               x_pos_top_right->Fill(rp_xpos_124*1000, 1);
               y_pos_top_right->Fill(rp_ypos_124*1000, 1);
           }    
           if (rp_right_accep_bottom){
               x_pos_bottom_right->Fill(rp_xpos_125*1000, 1);
               y_pos_bottom_right->Fill(rp_ypos_125*1000, 1);
           }    
        }
        if (proton_left_rec_sel){
           if (rp_left_accep_top){
               x_pos_top_left->Fill(rp_xpos_24*1000, 1);
               y_pos_top_left->Fill(rp_ypos_24*1000, 1);
           }    
           if (rp_left_accep_bottom){
               x_pos_bottom_left->Fill(rp_xpos_25*1000, 1);
               y_pos_bottom_left->Fill(rp_ypos_25*1000, 1);
           }    
        }

        if (proton_right_rec_sel && rp_right_sel){
            xi_rec_right_cut_nojet->Fill(xi_rec_proton_right, 0.94/corr_xi);
            t_rec_right_cut_nojet->Fill(fabs(t_rec_proton_right), 0.94/corr_xi);
            theta_x_rec_nojet->Fill(theta_x_minus_smear, 0.94/corr_xi);
            theta_y_rec_nojet->Fill(theta_y_minus_smear, 0.94/corr_xi);
            t_rec_vs_xi_rec_right_cut->Fill(fabs(t_rec_proton_right), xi_rec_proton_right, 0.94/corr_xi);
            theta_y_vs_x_rec_nojet->Fill(theta_x_minus_smear, theta_y_minus_smear, 0.94/corr_xi);
        }    

        if (proton_right_gen_sel){
            xi_gen_right_cut_nojet->Fill(xi_gen_proton_right, 1.);
            t_gen_right_cut_nojet->Fill(fabs(t_gen_proton_right), 1.);
            theta_x_gen_nojet->Fill(theta_x_minus, 1.);
            theta_y_gen_nojet->Fill(theta_y_minus, 1.);
            t_gen_vs_xi_gen_right_cut->Fill(fabs(t_gen_proton_right), xi_gen_proton_right, 1.);
            theta_y_vs_x_gen_nojet->Fill(theta_x_minus, theta_y_minus, 1.);
        }    
        
        //right  
        if (jet_rec_sel && proton_right_rec_sel && rp_right_sel){
            xi_cms_minus_totem_rec_right->Fill(xi_rec_cms_minus - xi_rec_proton_right, event_weight_minus);
            xi_rec_right->Fill(xi_rec_proton_right, event_weight_minus);
            if (xi_rec_cms_minus - xi_rec_proton_right<0){
                xi_rec_right_cut->Fill(xi_rec_proton_right, event_weight_minus);
                t_rec_right_cut->Fill(fabs(t_rec_proton_right), event_weight_minus);
                beta_right_cut->Fill(beta_rec_proton_right, event_weight_minus);
                log_x_rec_right_cut->Fill(log10(x_rec_right), event_weight_minus);
		            pt_jet1_right_cut->Fill(jet1_rec_pt, event_weight_minus);
		            pt_jet2_right_cut->Fill(jet2_rec_pt, event_weight_minus);
		            eta_jet1_right_cut->Fill(jet1_rec_eta, event_weight_minus);
		            eta_jet2_right_cut->Fill(jet2_rec_eta, event_weight_minus);
                th_x_rec_right->Fill(theta_x_minus_smear, event_weight_minus);
                th_y_rec_right->Fill(theta_y_minus_smear, event_weight_minus);
                delta_eta_jets_right_cut->Fill(jet1_rec_eta-jet2_rec_eta, event_weight_minus);
                delta_phi_jets_right_cut->Fill(jet1_rec_phi-jet2_rec_phi, event_weight_minus);
                mass_jj_rec_right->Fill(sqrt(mjj2_rec), event_weight_minus);
                mass_x_rec_right->Fill(4000*sqrt(xi_rec_proton_right), event_weight_minus);
                r_jj_rec_right->Fill(sqrt(mjj2_rec)/(4000*sqrt(xi_rec_proton_right)), event_weight_minus);

                if(!(jet_gen_sel && proton_right_gen_sel)){
                    t_minus_response.Fake(fabs(t_rec_proton_right),event_weight_minus*norm);
                    xi_minus_response.Fake(xi_rec_proton_right,event_weight_minus*norm);
                    logx_minus_response.Fake(log10(x_rec_right),event_weight_minus*norm);
                    t_minus_response_back.Miss(fabs(t_rec_proton_right),event_weight_minus*norm);
                    xi_minus_response_back.Miss(xi_rec_proton_right,event_weight_minus*norm);
                    logx_minus_response_back.Miss(log10(x_rec_right),event_weight_minus*norm);

                }
                if (jet_gen_sel && proton_right_gen_sel){
                    t_minus_response.Fill(fabs(t_rec_proton_right),fabs(t_gen_proton_right),event_weight_minus*norm);
                    xi_minus_response.Fill(xi_rec_proton_right,xi_gen_proton_right,event_weight_minus*norm);
                    logx_minus_response.Fill(log10(x_rec_right),log10(x_gen_right),event_weight_minus*norm);
                    t_minus_response_back.Fill(fabs(t_gen_proton_right),fabs(t_rec_proton_right),event_weight_minus*norm);
                    xi_minus_response_back.Fill(xi_gen_proton_right,xi_rec_proton_right,event_weight_minus*norm);
                    logx_minus_response_back.Fill(log10(x_gen_right),log10(x_rec_right),event_weight_minus*norm);
                    xi_minus_th2->Fill(xi_rec_proton_right, xi_gen_proton_right, event_weight_minus*norm);
                    logx_minus_th2->Fill(log10(x_rec_right),log10(x_gen_right), event_weight_minus*norm);
                    t_minus_th2->Fill(fabs(t_rec_proton_right), fabs(t_gen_proton_right), event_weight_minus*norm);
                }   
            }

        }
        if (jet_gen_sel && proton_right_gen_sel){
           t_gen_right_cut->Fill(fabs(t_gen_proton_right), event_weight_minus); 
           xi_gen_right_cut->Fill(xi_gen_proton_right, event_weight_minus); 
           log_x_gen_right_cut->Fill(log10(x_gen_right), event_weight_minus);
           th_x_gen_right->Fill(theta_x_minus, event_weight_minus);
           th_y_gen_right->Fill(theta_y_minus, event_weight_minus);
           if(!(jet_rec_sel && proton_right_rec_sel && rp_right_sel && xi_rec_cms_minus - xi_rec_proton_right<0)){
                    t_minus_response.Miss(fabs(t_gen_proton_right),event_weight_minus*norm);
                    xi_minus_response.Miss(xi_gen_proton_right,event_weight_minus*norm);
                    logx_minus_response.Miss(log10(x_gen_right),event_weight_minus*norm);
                    t_minus_response_back.Fake(fabs(t_gen_proton_right),event_weight_minus*norm);
                    xi_minus_response_back.Fake(xi_gen_proton_right,event_weight_minus*norm);
                    logx_minus_response_back.Fake(log10(x_gen_right),event_weight_minus*norm);
           
           }
        } 

        if (jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4){
           xi_gen_right_cut_bin_sasha->Fill(xi_gen_proton_right, 1); 
           xi_cms_gen_right_cut_bin_sasha->Fill(xi_gen_cms_minus, 1); 
           if (proton_right_gen_sel) pt_jet1_gen_right_cut->Fill(jet1_gen_pt, event_weight_minus);
        }   

        if (jet1_gen_pt>40 && jet2_gen_pt>40 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4){
           xi_gen_right_cut_bin_sasha_pt40->Fill(xi_gen_proton_right, 1); 
           xi_cms_gen_right_cut_bin_sasha_pt40->Fill(xi_gen_cms_minus, 1); 
        } 

        pt2_xi_gen_minus->Fill( fabs(t_gen_proton_right)/*(1 + xi_gen_proton_right)*/, xi_gen_proton_right, event_weight_minus);
        if (rp_right_sel) pt2_xi_gen_minus_rp->Fill( fabs(t_gen_proton_right)/**(1 + xi_gen_proton_right)*/, xi_gen_proton_right, event_weight_minus);

        //left
        if (jet_rec_sel && proton_left_rec_sel && rp_left_sel){
            // xi_cms_plus_totem_rec_left->Fill(xi_rec_cms_plus - xi_rec_proton_left, event_weight);
            xi_rec_left->Fill(xi_rec_proton_left, event_weight_plus);
            if (xi_rec_cms_plus - xi_rec_proton_left<0){
                t_rec_left_cut->Fill(fabs(t_rec_proton_left), event_weight_plus);
                xi_rec_left_cut->Fill(xi_rec_proton_left, event_weight_plus);
                beta_left_cut->Fill(beta_rec_proton_left, event_weight_plus);
                log_x_rec_left_cut->Fill(log10(x_rec_left), event_weight_plus);
		            pt_jet1_left_cut->Fill(jet1_rec_pt, event_weight_plus);
		            pt_jet2_left_cut->Fill(jet2_rec_pt, event_weight_plus);
		            eta_jet1_left_cut->Fill(jet1_rec_eta, event_weight_plus);
		            eta_jet2_left_cut->Fill(jet2_rec_eta, event_weight_plus);
                th_x_rec_left->Fill(theta_x_plus_smear, event_weight_plus);
                th_y_rec_left->Fill(theta_y_plus_smear, event_weight_plus);
                delta_eta_jets_left_cut->Fill(jet1_rec_eta-jet2_rec_eta, event_weight_plus);
                delta_phi_jets_left_cut->Fill(jet1_rec_phi-jet2_rec_phi, event_weight_plus);
                mass_jj_rec_left->Fill(sqrt(mjj2_rec), event_weight_plus);
                mass_x_rec_left->Fill(4000*sqrt(xi_rec_proton_left), event_weight_plus);
                r_jj_rec_left->Fill(sqrt(mjj2_rec)/(4000*sqrt(xi_rec_proton_left)), event_weight_plus);

                if (!(jet_gen_sel && proton_left_gen_sel)){
                   t_plus_response.Fake(fabs(t_rec_proton_left),event_weight_plus*norm);
                   xi_plus_response.Fake(xi_rec_proton_left,event_weight_plus*norm);
                   logx_plus_response.Fake(log10(x_rec_left),event_weight_plus*norm);
                }
                if (jet_gen_sel && proton_left_gen_sel){
                   t_plus_response.Fill(fabs(t_rec_proton_left),fabs(t_gen_proton_left),event_weight_plus*norm);
                   xi_plus_response.Fill(xi_rec_proton_left,xi_gen_proton_left,event_weight_plus*norm);
                   logx_plus_response.Fill(log10(x_rec_left),log10(x_gen_left),event_weight_plus*norm);
                    xi_plus_th2->Fill(xi_rec_proton_left, xi_gen_proton_left, event_weight_plus*norm);
                    logx_plus_th2->Fill(log10(x_rec_left),log10(x_gen_left), event_weight_plus*norm);
                    t_plus_th2->Fill(fabs(t_rec_proton_left), fabs(t_gen_proton_left), event_weight_plus*norm);
                }
             }
	 }   
             t_gen_left_bin->Fill(fabs(t_rec_proton_left), event_weight_plus);

         if (jet_gen_sel && proton_left_gen_sel){
             t_gen_left_cut->Fill(fabs(t_gen_proton_left), event_weight_plus);
             xi_gen_left_cut->Fill(xi_gen_proton_left, event_weight_plus);
             log_x_gen_left_cut->Fill(log10(x_gen_left), event_weight_plus);
             th_x_gen_left->Fill(theta_x_plus, event_weight_plus);
             th_y_gen_left->Fill(theta_y_plus, event_weight_plus);
             if (!(jet_rec_sel && proton_left_rec_sel && rp_left_sel && xi_rec_cms_plus - xi_rec_proton_left<0)){
                   t_plus_response.Miss(fabs(t_gen_proton_left), event_weight_plus*norm);
                   xi_plus_response.Miss(xi_gen_proton_left,event_weight_plus*norm);
                   logx_plus_response.Miss(log10(x_gen_left),event_weight_plus*norm);
             }
	 
         }

        if (jet1_gen_pt>20 && jet2_gen_pt>20 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4){
             xi_cms_gen_left_cut_bin_sasha->Fill(xi_gen_cms_plus, 1); 
             xi_gen_left_cut_bin_sasha->Fill(xi_gen_proton_left, event_weight_plus);
        }     
    }


pt_jet1_gen_right_cut->Scale(norm);
    xi_cms_minus_totem_rec_right->Scale(norm);
    xi_rec_right->Scale(norm);
    xi_rec_right_cut->Scale(norm);
    xi_gen_right_cut->Scale(norm);
    xi_gen_right_cut_bin_sasha->Scale(norm);
    xi_cms_gen_right_cut_bin_sasha->Scale(norm);
    xi_cms_gen_right_cut_bin_sasha_pt40->Scale(norm);
    t_rec_right_cut->Scale(norm);
    t_gen_right_cut->Scale(norm);
    beta_right_cut->Scale(norm);  
    log_x_rec_right_cut->Scale(norm);  
    log_x_gen_right_cut->Scale(norm);  
    pt_jet1_right_cut->Scale(norm);  
    pt_jet2_right_cut->Scale(norm);  
    eta_jet1_right_cut->Scale(norm);  
    eta_jet1_left_cut->Scale(norm);  
    eta_jet2_right_cut->Scale(norm);  
    th_x_rec_right->Scale(norm);  
    th_y_rec_right->Scale(norm);  
    th_x_rec_left->Scale(norm);  
    th_y_rec_left->Scale(norm);  
    th_x_gen_right->Scale(norm);  
    th_y_gen_right->Scale(norm);  
    th_x_gen_left->Scale(norm);  
    th_y_gen_left->Scale(norm);  
    delta_eta_jets_right_cut->Scale(norm); 
    delta_eta_jets_left_cut->Scale(norm); 
    delta_phi_jets_right_cut->Scale(norm); 
    delta_phi_jets_left_cut->Scale(norm); 
    mass_jj_rec_right->Scale(norm); 
    mass_x_rec_right->Scale(norm); 
    r_jj_rec_right->Scale(norm); 
    mass_jj_rec_left->Scale(norm); 
    mass_x_rec_left->Scale(norm); 
    r_jj_rec_left->Scale(norm); 

    xi_cms_minus_totem_rec_left->Scale(norm);
    xi_rec_left->Scale(norm);
    xi_rec_left_cut->Scale(norm);
    xi_gen_left_cut->Scale(norm);
    xi_gen_left_cut_bin_sasha->Scale(norm);
    xi_cms_gen_left_cut_bin_sasha->Scale(norm);
    t_rec_left_cut->Scale(norm);
    t_gen_left_cut->Scale(norm);
    beta_left_cut->Scale(norm);  
    log_x_rec_left_cut->Scale(norm);  
    log_x_gen_left_cut->Scale(norm);  
    pt_jet1_left_cut->Scale(norm);  
    pt_jet2_left_cut->Scale(norm);  
    eta_jet2_left_cut->Scale(norm);  

    ///----------------------------------- Unfolding -----------------------------------

    TH1F* t_right_cut = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_right_cut_bh = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_right_cut_noeff = new TH1F("t_right_cut_noeff","", 8, tbins);
    TH1F* t_left_cut = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_left_cut_bh = new TH1F("t_right_cut","", 8, tbins);
    TH1F* t_left_cut_noeff = new TH1F("t_right_cut_noeff","", 8, tbins);
    TH1F* xi_right_cut = new TH1F("xi_right_cut","",11,xi_bins);
    TH1F* xi_right_cut_bh = new TH1F("xi_right_cut","",11,xi_bins);
    TH1F* xi_right_cut_noeff = new TH1F("xi_right_cut_noeff","",11,xi_bins);
    TH1F* xi_left_cut = new TH1F("xi_left_cut","",11,xi_bins);
    TH1F* xi_left_cut_bh = new TH1F("xi_left_cut","",11,xi_bins);
    TH1F* xi_left_cut_noeff = new TH1F("xi_left_cut_noeff","",11,xi_bins);
    TH1F* log_x_right_cut = new TH1F("log_x_right_cut","",15, -4, 0);
    TH1F* log_x_right_cut_bh = new TH1F("log_x_right_cut","",15, -4, 0);
    TH1F* log_x_right_cut_noeff = new TH1F("log_x_right_cut_noeff","",15, -4, 0);
    TH1F* log_x_right_noeff = new TH1F("log_x_right_noeff","",15, -4, 0);
    TH1F* log_x_left_cut = new TH1F("log_x_left_cut","",15, -4, 0);
    TH1F* log_x_left_cut_bh = new TH1F("log_x_left_cut","",15, -4, 0);
    TH1F* log_x_left_cut_noeff = new TH1F("log_x_left_cut_noeff","",15, -4, 0);
    TH1F* log_x_left_noeff = new TH1F("log_x_left_noeff","",15, -4, 0);
    TH1F* h1; TH1F* xi_cms_right_cut_sasha; TH1F* xi_cms_left_cut_sasha;
    TH1F* h2; 
    TH1F* h3; TH1F* h4; TH1F* h5; TH1F* h6; TH1F* h7; TH1F* h8; TH1F* h9; TH1F* h10; TH1F* h11; TH1F* h12; TH1F* h13; TH1F* h14; TH1F* h15; TH1F* h16; TH1F* h17; TH1F* h18; TH1F* h19; TH1F* h20; 
    TH1F* h21; TH1F* h22; TH1F* h23; TH1F* h24; TH1F* h25; TH1F* h26; TH1F* h27; TH1F* h28; TH1F* h29; TH1F* h30; TH1F* h31; TH1F* h32; TH1F* h33; TH1F* h34; TH1F* h35; TH1F* h36; TH1F* h37; TH1F* h38;
    TH1F* h39; TH1F* h40; TH1F* h41; 
    TH1F* h42; TH1F* h43; TH1F* h44; TH1F* h45; TH1F* h46; TH1F* h47; TH1F* h48; TH1F* h49; TH1F* h50; TH1F* h51; TH1F* h52; TH1F* h53; TH1F* h54; TH1F* h55; TH1F* h56;
    TH1F* h60; TH1F* h57; TH1F* h58; TH1F* h59; TH1F* h61; TH1F* h62;
    TH1F* th_x_right; TH1F* th_y_right; TH1F* th_x_left; TH1F* th_y_left;
    TH1F* sasha = new TH1F("sigma_xi_cms_right_sasha","",8, bin_sasha);

    if (unc_rp && !rp_x_min && !rp_x_max){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42, h57, h58, h59, h60,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, xi_cms_right_cut_sasha, xi_cms_left_cut_sasha, sasha, h55, h56, th_x_right, th_y_right, th_x_left, th_y_left, h57, h58, h59, h60, h61, h62, true, true,false,false, false, false, false);}

    if (unc_rp == false && rp_x_min == true && rp_x_max == false){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42, h57, h58, h59, h60,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, xi_cms_right_cut_sasha, xi_cms_left_cut_sasha, sasha, h55, h56, th_x_right, th_y_right, th_x_left, th_y_left, h57, h58, h59, h60, h61, h62, false, true,true,false, false, false, false);}

    if (unc_rp == false && rp_x_min == false && rp_x_max == true){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42, h57, h58, h59, h60,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, xi_cms_right_cut_sasha, xi_cms_left_cut_sasha, sasha, h55, h56, th_x_right, th_y_right, th_x_left, th_y_left, h57, h58, h59, h60, h61, h62, false, true,false,true, false, false, false);}

    if (!unc_rp && !rp_x_min && !rp_x_max){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42, h57, h58, h59, h60,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, xi_cms_right_cut_sasha, xi_cms_left_cut_sasha, sasha, h55, h56, th_x_right, th_y_right, th_x_left, th_y_left, h57, h58, h59, h60, h61, h62,  false, true, false, false, false, false, false);}

    if (sys_xi_up == true){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42, h57, h58, h59, h60,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, xi_cms_right_cut_sasha, xi_cms_left_cut_sasha, sasha, h55, h56, th_x_right, th_y_right, th_x_left, th_y_left, h57, h58, h59, h60, h61, h62, false, true, false, false, true, false, false);}

    if (sys_xi_dw == true){
        data(h1, h2, h3, h4, h5, 
          h6, h7, h8, xi_right_cut, xi_right_cut_bh, h11, xi_right_cut_noeff, h13,
          xi_left_cut, xi_left_cut_bh, h16, xi_left_cut_noeff, t_right_cut, t_right_cut_bh, h19, h20, 
          h21, t_right_cut_noeff, h23, t_left_cut, t_left_cut_bh, h25, h26, h27,
          t_left_cut_noeff, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,  h57, h58, h59, h60,
          h43, log_x_right_noeff, log_x_right_cut, log_x_right_cut_bh, log_x_right_cut_noeff, h47, log_x_left_noeff, log_x_left_cut, log_x_left_cut_bh, log_x_left_cut_noeff, 
          h51, h52, xi_cms_right_cut_sasha, xi_cms_left_cut_sasha, sasha, h55, h56, th_x_right, th_y_right, th_x_left, th_y_left, h57, h58, h59, h60, h61, h62, false, true, false, false, false, true, false);}
        

    TH1F* t_right_cut_herabackg = (TH1F*) t_right_cut->Clone();
    TH1F* t_left_cut_herabackg = (TH1F*) t_left_cut->Clone();
    TH1F* xi_right_cut_herabackg = (TH1F*) xi_right_cut->Clone();
    TH1F* xi_left_cut_herabackg = (TH1F*) xi_left_cut->Clone();
    TH1F* log_x_right_cut_herabackg = (TH1F*) log_x_right_cut->Clone();
    TH1F* log_x_left_cut_herabackg = (TH1F*) log_x_left_cut->Clone();
    TH1F* t_right_cut_backg;
    TH1F* t_left_cut_backg;
    TH1F* xi_right_cut_backg;
    TH1F* xi_left_cut_backg;
    TH1F* log_x_right_cut_backg;
    TH1F* log_x_left_cut_backg;

    pu_zb (h1, h2, h3, h4, h5, h6, xi_right_cut_backg, xi_left_cut_backg, h7, h8, t_right_cut_backg, t_left_cut_backg, log_x_right_cut_backg, log_x_left_cut_backg,h9, h10, false, true, false);

    t_right_cut_herabackg->Add(t_right_cut_bh, -1); 
    t_left_cut_herabackg->Add(t_left_cut_bh, -1); 
    xi_right_cut_herabackg->Add(xi_right_cut_bh, -1); 
    xi_left_cut_herabackg->Add(xi_right_cut_bh, -1); 
    log_x_right_cut_herabackg->Add(log_x_right_cut_bh, -1); 
    log_x_left_cut_herabackg->Add(log_x_left_cut_bh, -1); 

    t_right_cut->Add(t_right_cut_backg, -1); 
    t_left_cut->Add(t_left_cut_backg, -1); 
    xi_right_cut->Add(xi_right_cut_backg, -1); 
    xi_left_cut->Add(xi_left_cut_backg, -1); 
    log_x_right_cut->Add(log_x_right_cut_backg, -1); 
    log_x_left_cut->Add(log_x_left_cut_backg, -1); 

    int n_iter_t_right, n_iter_t_left, n_iter_xi_right, n_iter_xi_left, n_iter_logx_right, n_iter_logx_left;

    if (!iter_variation){
       if (tune == "Tune4C"){
        n_iter_t_right = 10;
        n_iter_t_left = 12; 
        n_iter_xi_right = 10;
        n_iter_xi_left = 6;
        n_iter_logx_right = 5;
        n_iter_logx_left = 6;
       }
       if (tune == "CUETP8M1"){
         n_iter_t_right = 12;
        n_iter_t_left = 10;
        n_iter_xi_right = 9;
        n_iter_xi_left = 9;
        n_iter_logx_right = 5;
        n_iter_logx_left = 6;
       } 
    } 
    else{
       if (tune == "Tune4C"){
        n_iter_t_right = 20;
        n_iter_t_left = 25;
        n_iter_xi_right = 15;
        n_iter_xi_left = 9;
        n_iter_logx_right = 7;
        n_iter_logx_left = 9; 
       }
       if (tune == "CUETP8M1"){
        n_iter_t_right = 22;
        n_iter_t_left = 18;
        n_iter_xi_right = 15;
        n_iter_xi_left = 15;
        n_iter_logx_right = 7;
        n_iter_logx_left = 8;
       }
    }


    RooUnfoldBayes unfold_right_t (&t_minus_response, t_right_cut, n_iter_t_right);
    t_right_unfolded = (TH1F*)unfold_right_t.Hreco();
    RooUnfoldBayes unfold_right_t_noeff (&t_minus_response, t_right_cut_noeff, n_iter_t_right);
    t_right_unfolded_noeff = (TH1F*)unfold_right_t_noeff.Hreco();
    RooUnfoldBinByBin unfold_right_t_binbybin (&t_minus_response, t_right_cut);
    t_right_unfolded_binbybin = (TH1F*)unfold_right_t_binbybin.Hreco();
    RooUnfoldBayes unfold_right_t_herabackg (&t_minus_response, t_right_cut_herabackg, n_iter_t_right);
    t_right_unfolded_herabackg = (TH1F*)unfold_right_t_herabackg.Hreco();
     
    RooUnfoldBayes unfold_left_t (&t_plus_response, t_left_cut, n_iter_t_left);
    t_left_unfolded = (TH1F*)unfold_left_t.Hreco();
    RooUnfoldBayes unfold_left_t_noeff (&t_plus_response, t_left_cut_noeff, n_iter_t_left);
    t_left_unfolded_noeff = (TH1F*)unfold_left_t_noeff.Hreco();
    RooUnfoldBinByBin unfold_left_t_binbybin (&t_plus_response, t_left_cut);
    t_left_unfolded_binbybin = (TH1F*)unfold_left_t_binbybin.Hreco();
    RooUnfoldBayes unfold_left_t_herabackg (&t_plus_response, t_left_cut_herabackg, n_iter_t_left);
    t_left_unfolded_herabackg = (TH1F*)unfold_left_t_herabackg.Hreco();

    RooUnfoldBayes unfold_right_xi (&xi_minus_response, xi_right_cut, n_iter_xi_right);
    xi_right_unfolded = (TH1F*)unfold_right_xi.Hreco();
    RooUnfoldBayes unfold_right_xi_noeff (&xi_minus_response, xi_right_cut_noeff, n_iter_xi_right);
    xi_right_unfolded_noeff = (TH1F*)unfold_right_xi_noeff.Hreco();
    RooUnfoldBinByBin unfold_right_xi_binbybin (&xi_minus_response, xi_right_cut);
    xi_right_unfolded_binbybin = (TH1F*)unfold_right_xi_binbybin.Hreco();
    RooUnfoldBayes unfold_right_xi_herabackg (&xi_minus_response, xi_right_cut_herabackg, n_iter_xi_right);
    xi_right_unfolded_herabackg = (TH1F*)unfold_right_xi_herabackg.Hreco();

    RooUnfoldBayes unfold_left_xi (&xi_plus_response, xi_left_cut, n_iter_xi_left);
    xi_left_unfolded = (TH1F*)unfold_left_xi.Hreco();
    RooUnfoldBayes unfold_left_xi_noeff (&xi_plus_response, xi_left_cut_noeff, n_iter_xi_left);
    xi_left_unfolded_noeff = (TH1F*)unfold_left_xi_noeff.Hreco();
    RooUnfoldBinByBin unfold_left_xi_binbybin (&xi_plus_response, xi_left_cut);
    xi_left_unfolded_binbybin = (TH1F*)unfold_left_xi_binbybin.Hreco();
    RooUnfoldBayes unfold_left_xi_herabackg (&xi_plus_response, xi_left_cut_herabackg, n_iter_xi_left);
    xi_left_unfolded_herabackg = (TH1F*)unfold_left_xi_herabackg.Hreco();

    RooUnfoldBayes unfold_right_logx (&logx_minus_response, log_x_right_cut, n_iter_logx_right);
    logx_right_unfolded = (TH1F*)unfold_right_logx.Hreco();
    RooUnfoldBayes unfold_right_logx_cut_noeff (&logx_minus_response, log_x_right_cut_noeff, n_iter_logx_right);
    logx_right_cut_unfolded_noeff = (TH1F*)unfold_right_logx_cut_noeff.Hreco();
    RooUnfoldBayes unfold_right_logx_noeff (&logx_minus_response, log_x_right_noeff, n_iter_logx_right);
    logx_right_unfolded_noeff = (TH1F*)unfold_right_logx_noeff.Hreco();
    RooUnfoldBinByBin unfold_right_logx_binbybin (&logx_minus_response, log_x_right_cut);
    logx_right_unfolded_binbybin = (TH1F*)unfold_right_logx_binbybin.Hreco();
    RooUnfoldBayes unfold_right_logx_herabackg (&logx_minus_response, log_x_right_cut_herabackg, n_iter_logx_right);
    logx_right_unfolded_herabackg = (TH1F*)unfold_right_logx_herabackg.Hreco();

    RooUnfoldBayes unfold_left_logx (&logx_plus_response, log_x_left_cut, n_iter_logx_left);
    logx_left_unfolded = (TH1F*)unfold_left_logx.Hreco();
    RooUnfoldBayes unfold_left_logx_cut_noeff (&logx_plus_response, log_x_left_cut_noeff, n_iter_logx_left);
    logx_left_cut_unfolded_noeff = (TH1F*)unfold_left_logx_cut_noeff.Hreco();
    RooUnfoldBayes unfold_left_logx_noeff (&logx_plus_response, log_x_left_noeff, n_iter_logx_left);
    logx_left_unfolded_noeff = (TH1F*)unfold_left_logx_noeff.Hreco();
    RooUnfoldBinByBin unfold_left_logx_binbybin (&logx_plus_response, log_x_left_cut);
    logx_left_unfolded_binbybin = (TH1F*)unfold_left_logx_binbybin.Hreco();
    RooUnfoldBayes unfold_left_logx_herabackg (&logx_plus_response, log_x_left_cut_herabackg, n_iter_logx_left);
    logx_left_unfolded_herabackg = (TH1F*)unfold_left_logx_herabackg.Hreco();

TH1F* t_right_cut_up; TH1F* t_left_cut_up; TH1F* xi_right_cut_up; TH1F* xi_left_cut_up; TH1F* log_x_minus_up; TH1F* log_x_plus_up; TH1F* log_x_right_cut_up; TH1F* log_x_left_cut_up;
TH1F* t_right_cut_dw; TH1F* t_left_cut_dw; TH1F* xi_right_cut_dw; TH1F* xi_left_cut_dw; TH1F* log_x_minus_dw; TH1F* log_x_plus_dw; TH1F* log_x_right_cut_dw; TH1F* log_x_left_cut_dw;
TH1F* t_right_cut_up_pf; TH1F* t_left_cut_up_pf; TH1F* xi_right_cut_up_pf; TH1F* xi_left_cut_up_pf; TH1F* log_x_minus_up_pf; TH1F* log_x_plus_up_pf; TH1F* log_x_right_cut_up_pf; TH1F* log_x_left_cut_up_pf;
TH1F* t_right_cut_dw_pf; TH1F* t_left_cut_dw_pf; TH1F* xi_right_cut_dw_pf; TH1F* xi_left_cut_dw_pf; TH1F* log_x_minus_dw_pf; TH1F* log_x_plus_dw_pf; TH1F* log_x_right_cut_dw_pf; TH1F* log_x_left_cut_dw_pf;
unc(t_right_cut_up, t_left_cut_up, xi_right_cut_up, xi_left_cut_up, log_x_minus_up, log_x_plus_up, log_x_right_cut_up, log_x_left_cut_up, true, true, true);
unc(t_right_cut_dw, t_left_cut_dw, xi_right_cut_dw, xi_left_cut_dw, log_x_minus_dw, log_x_plus_dw, log_x_right_cut_dw, log_x_left_cut_dw, true, false, true);
unc(t_right_cut_up_pf, t_left_cut_up_pf, xi_right_cut_up_pf, xi_left_cut_up_pf, log_x_minus_up_pf, log_x_plus_up_pf, log_x_right_cut_up_pf, log_x_left_cut_up_pf, false, true, true);
unc(t_right_cut_dw_pf, t_left_cut_dw_pf, xi_right_cut_dw_pf, xi_left_cut_dw_pf, log_x_minus_dw_pf, log_x_plus_dw_pf, log_x_right_cut_dw_pf, log_x_left_cut_dw_pf, false, false, true);

    RooUnfoldBayes unfold_right_t_up (&t_minus_response, t_right_cut_up, n_iter_t_right);
    t_right_unfolded_up = (TH1F*)unfold_right_t_up.Hreco();
    RooUnfoldBayes unfold_right_t_dw (&t_minus_response, t_right_cut_dw, n_iter_t_right);
    t_right_unfolded_dw = (TH1F*)unfold_right_t_dw.Hreco();

    RooUnfoldBayes unfold_left_t_up (&t_plus_response, t_left_cut_up, n_iter_t_left);
    t_left_unfolded_up = (TH1F*)unfold_left_t_up.Hreco();
    RooUnfoldBayes unfold_left_t_dw (&t_plus_response, t_left_cut_dw, n_iter_t_left);
    t_left_unfolded_dw = (TH1F*)unfold_left_t_dw.Hreco();

    RooUnfoldBayes unfold_right_xi_up (&xi_minus_response, xi_right_cut_up, n_iter_xi_right);
    xi_right_unfolded_up = (TH1F*)unfold_right_xi_up.Hreco();
    RooUnfoldBayes unfold_right_xi_dw (&xi_minus_response, xi_right_cut_dw, n_iter_xi_right);
    xi_right_unfolded_dw = (TH1F*)unfold_right_xi_dw.Hreco();

    RooUnfoldBayes unfold_left_xi_up (&xi_plus_response, xi_left_cut_up, n_iter_xi_left);
    xi_left_unfolded_up = (TH1F*)unfold_left_xi_up.Hreco();
    RooUnfoldBayes unfold_left_xi_dw (&xi_plus_response, xi_left_cut_dw, n_iter_xi_left);
    xi_left_unfolded_dw = (TH1F*)unfold_left_xi_dw.Hreco();

    RooUnfoldBayes unfold_right_log_x_up (&logx_minus_response, log_x_right_cut_up, n_iter_logx_right);
    log_x_right_unfolded_up = (TH1F*)unfold_right_log_x_up.Hreco();
    RooUnfoldBayes unfold_right_log_x_dw (&logx_minus_response, log_x_right_cut_dw, n_iter_logx_right);
    log_x_right_unfolded_dw = (TH1F*)unfold_right_log_x_dw.Hreco();

    RooUnfoldBayes unfold_left_log_x_up (&logx_plus_response, log_x_left_cut_up, n_iter_logx_left);
    log_x_left_unfolded_up = (TH1F*)unfold_left_log_x_up.Hreco();
    RooUnfoldBayes unfold_left_log_x_dw (&logx_plus_response, log_x_left_cut_dw, n_iter_logx_left);
    log_x_left_unfolded_dw = (TH1F*)unfold_left_log_x_dw.Hreco();

    //pf
    RooUnfoldBayes unfold_right_t_up_pf (&t_minus_response, t_right_cut_up_pf, n_iter_t_right);
    t_right_unfolded_up_pf = (TH1F*)unfold_right_t_up_pf.Hreco();
    RooUnfoldBayes unfold_right_t_dw_pf (&t_minus_response, t_right_cut_dw_pf, n_iter_t_right);
    t_right_unfolded_dw_pf = (TH1F*)unfold_right_t_dw_pf.Hreco();

    RooUnfoldBayes unfold_left_t_up_pf (&t_plus_response, t_left_cut_up_pf, n_iter_t_left);
    t_left_unfolded_up_pf = (TH1F*)unfold_left_t_up_pf.Hreco();
    RooUnfoldBayes unfold_left_t_dw_pf (&t_plus_response, t_left_cut_dw_pf, n_iter_t_left);
    t_left_unfolded_dw_pf = (TH1F*)unfold_left_t_dw_pf.Hreco();

    RooUnfoldBayes unfold_right_xi_up_pf (&xi_minus_response, xi_right_cut_up_pf, n_iter_xi_right);
    xi_right_unfolded_up_pf = (TH1F*)unfold_right_xi_up_pf.Hreco();
    RooUnfoldBayes unfold_right_xi_dw_pf (&xi_minus_response, xi_right_cut_dw_pf, n_iter_xi_right);
    xi_right_unfolded_dw_pf = (TH1F*)unfold_right_xi_dw_pf.Hreco();

    RooUnfoldBayes unfold_left_xi_up_pf (&xi_plus_response, xi_left_cut_up_pf, n_iter_xi_left);
    xi_left_unfolded_up_pf = (TH1F*)unfold_left_xi_up_pf.Hreco();
    RooUnfoldBayes unfold_left_xi_dw_pf (&xi_plus_response, xi_left_cut_dw_pf, n_iter_xi_left);
    xi_left_unfolded_dw_pf = (TH1F*)unfold_left_xi_dw_pf.Hreco();

    RooUnfoldBayes unfold_right_log_x_up_pf (&logx_minus_response, log_x_right_cut_up_pf, n_iter_logx_right);
    log_x_right_unfolded_up_pf = (TH1F*)unfold_right_log_x_up_pf.Hreco();
    RooUnfoldBayes unfold_right_log_x_dw_pf (&logx_minus_response, log_x_right_cut_dw_pf, n_iter_logx_right);
    log_x_right_unfolded_dw_pf = (TH1F*)unfold_right_log_x_dw_pf.Hreco();

    RooUnfoldBayes unfold_left_log_x_up_pf (&logx_plus_response, log_x_left_cut_up_pf, n_iter_logx_left);
    log_x_left_unfolded_up_pf = (TH1F*)unfold_left_log_x_up_pf.Hreco();
    RooUnfoldBayes unfold_left_log_x_dw_pf (&logx_plus_response, log_x_left_cut_dw_pf, n_iter_logx_left);
    log_x_left_unfolded_dw_pf = (TH1F*)unfold_left_log_x_dw_pf.Hreco();


cout<<"pythia8"+id<<"  right: "<<t_gen_right_cut->Integral()<<"   left: "<<t_gen_left_cut->Integral()<<endl;




//    TFile * outfile = new TFile("pythia8_4c_bottom_line_test.root", "RECREATE");

// TH2F* chi2_unfolded_vs_iter_t = new TH2F("chi2_unfolded_vs_iter_t","",100,0,40,100,0.0001,10);
// TH2F* chi2_unfolded_vs_iter_t_left = new TH2F("chi2_unfolded_vs_iter_t_left","",100,0,40,100,0.0001,10);
// TH2F* chi2_unfolded_vs_iter_xi = new TH2F("chi2_unfolded_vs_iter_xi","",100,0,40,100,0.0001,10);
// TH2F* chi2_unfolded_vs_iter_xi_left = new TH2F("chi2_unfolded_vs_iter_xi_left","",100,0,40,100,0.0001,10);
// TH2F* chi2_unfolded_vs_iter_logx = new TH2F("chi2_unfolded_vs_iter_logx","",100,0,40,100,0.0001,10);
// TH2F* chi2_unfolded_vs_iter_logx_left = new TH2F("chi2_unfolded_vs_iter_logx_left","",100,0,40,100,0.0001,10);
// TH2F* delta_chi2_unfolded_vs_iter_t = new TH2F("t","",2000,0,40,2000,0.0001,1);
// TH2F* delta_chi2_unfolded_vs_iter_t_left = new TH2F("t_left","",12000,0,40,2000,0.0001,1);
// TH2F* delta_chi2_unfolded_vs_iter_xi = new TH2F("xi","",2000,0,40,2000,0.0001,1);
// TH2F* delta_chi2_unfolded_vs_iter_xi_left = new TH2F("xi_left","",2000,0,40,2000,0.0001,1);
// TH2F* delta_chi2_unfolded_vs_iter_logx = new TH2F("logx","",2000,0,40,2000,0.0001,1);
// TH2F* delta_chi2_unfolded_vs_iter_logx_left = new TH2F("logx_left","",2000,0,40,2000,0.0001,1);
// double chi2_unfolded_t [40]; 
// double chi2_unfolded_t_left[40]; 
// double chi2_unfolded_xi[40]; 
// double chi2_unfolded_xi_left[40]; 
// double chi2_unfolded_logx[40]; 
// double chi2_unfolded_logx_left[40]; 

// for (int i = 1; i<40; ++i){
//     RooUnfoldBayes t_right_unfold_per_iter (&t_minus_response, t_right_cut, i);
//     TH1F* t_right_unfolded_per_iter = (TH1F*)t_right_unfold_per_iter.Hreco();
//     RooUnfoldBayes t_left_unfold_per_iter (&t_plus_response, t_left_cut, i);
//     TH1F* t_left_unfolded_per_iter = (TH1F*)t_left_unfold_per_iter.Hreco();

//     RooUnfoldBayes xi_right_unfold_per_iter (&xi_minus_response, xi_right_cut, i);
//     TH1F* xi_right_unfolded_per_iter = (TH1F*)xi_right_unfold_per_iter.Hreco();
//     RooUnfoldBayes xi_left_unfold_per_iter (&xi_plus_response, xi_left_cut, i);
//     TH1F* xi_left_unfolded_per_iter = (TH1F*)xi_left_unfold_per_iter.Hreco();

//     RooUnfoldBayes logx_right_unfold_per_iter (&logx_minus_response, log_x_right_cut, i);
//     TH1F* logx_right_unfolded_per_iter = (TH1F*)logx_right_unfold_per_iter.Hreco();
//     RooUnfoldBayes logx_left_unfold_per_iter (&logx_plus_response, log_x_left_cut, i);
//     TH1F* logx_left_unfolded_per_iter = (TH1F*)logx_left_unfold_per_iter.Hreco();

//     chi2_unfolded_t[i] = 0; 
//     chi2_unfolded_t_left[i] = 0; 
//     chi2_unfolded_xi[i] = 0; 
//     chi2_unfolded_xi_left[i] = 0; 
//     chi2_unfolded_logx[i] = 0; 
//     chi2_unfolded_logx_left[i] = 0; 
//     for (int j = 1; j<=t_left_unfolded->GetNbinsX()-1; ++j){
//          if(t_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_t[i] += pow(t_right_unfolded_per_iter->GetBinContent(j) - t_gen_right_cut->GetBinContent(j), 2)/(pow(t_right_unfolded_per_iter->GetBinError(j),2));
//          if(t_left_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_t_left[i] += pow(t_left_unfolded_per_iter->GetBinContent(j) - t_gen_left_cut->GetBinContent(j), 2)/(pow(t_left_unfolded_per_iter->GetBinError(j),2));
//     } 
//     for (int j = 1; j<=xi_left_unfolded->GetNbinsX(); ++j){
//          if(xi_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_xi[i] += pow(xi_right_unfolded_per_iter->GetBinContent(j) - xi_gen_right_cut->GetBinContent(j), 2)/(pow(xi_right_unfolded_per_iter->GetBinError(j),2));
//          if(xi_left_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_xi_left[i] += pow(xi_left_unfolded_per_iter->GetBinContent(j) - xi_gen_left_cut->GetBinContent(j), 2)/(pow(xi_left_unfolded_per_iter->GetBinError(j),2));
//     }
//     for (int j = 4; j<11; ++j){
//          if(logx_right_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_logx[i] += pow(logx_right_unfolded_per_iter->GetBinContent(j) - log_x_gen_right_cut->GetBinContent(j), 2)/(pow(logx_right_unfolded_per_iter->GetBinError(j),2));
//          if(logx_left_unfolded_per_iter->GetBinError(j)>0) chi2_unfolded_logx_left[i] += pow(logx_left_unfolded_per_iter->GetBinContent(j) - log_x_gen_left_cut->GetBinContent(j), 2)/(pow(logx_left_unfolded_per_iter->GetBinError(j),2));
//     }
//     chi2_unfolded_vs_iter_t->Fill(i, log10(chi2_unfolded_t[i])); 
//     chi2_unfolded_vs_iter_t_left->Fill(i, log10(chi2_unfolded_t_left[i])); 
//     chi2_unfolded_vs_iter_xi->Fill(i, log10(chi2_unfolded_xi[i])); 
//     chi2_unfolded_vs_iter_xi_left->Fill(i, log10(chi2_unfolded_xi_left[i])); 
//     chi2_unfolded_vs_iter_logx->Fill(i, log10(chi2_unfolded_logx[i])); 
//     chi2_unfolded_vs_iter_logx_left->Fill(i, log10(chi2_unfolded_logx_left[i])); 
// }

// double chi2_smeared_t = 0;
// double chi2_smeared_t_left = 0;
// double chi2_smeared_xi = 0;
// double chi2_smeared_xi_left = 0;
// double chi2_smeared_logx = 0;
// double chi2_smeared_logx_left = 0;
//     for (int j = 1; j<=t_left_cut->GetNbinsX()-1; ++j){
//          if(t_right_cut->GetBinError(j)>0) chi2_smeared_t += pow(t_right_cut->GetBinContent(j) - t_rec_right_cut->GetBinContent(j), 2)/(pow(t_right_cut->GetBinError(j),2));
//          if(t_left_cut->GetBinError(j)>0) chi2_smeared_t_left += pow(t_left_cut->GetBinContent(j) - t_rec_left_cut->GetBinContent(j), 2)/(pow(t_left_cut->GetBinError(j),2));
//     } 
//     for (int j = 1; j<=xi_left_cut->GetNbinsX(); ++j){
//          if(xi_right_cut->GetBinError(j)>0) chi2_smeared_xi += pow(xi_right_cut->GetBinContent(j) - xi_rec_right_cut->GetBinContent(j), 2)/(pow(xi_right_cut->GetBinError(j),2));
//          if(xi_left_cut->GetBinError(j)>0) chi2_smeared_xi_left += pow(xi_left_cut->GetBinContent(j) - xi_rec_left_cut->GetBinContent(j), 2)/(pow(xi_left_cut->GetBinError(j),2));
//     }
//     for (int j = 4; j<11; ++j){
//          if(log_x_right_cut->GetBinError(j)>0) chi2_smeared_logx += pow(log_x_right_cut->GetBinContent(j) - log_x_rec_right_cut->GetBinContent(j), 2)/(pow(log_x_right_cut->GetBinError(j),2));
//          if(log_x_left_cut->GetBinError(j)>0) chi2_smeared_logx_left += pow(log_x_left_cut->GetBinContent(j) - log_x_rec_left_cut->GetBinContent(j), 2)/(pow(log_x_left_cut->GetBinError(j),2));
//     }

// cout<<"t_right: "<<chi2_smeared_t<<endl;
// cout<<"t_left: "<<chi2_smeared_t_left<<endl;
// cout<<"xi_right: "<<chi2_smeared_xi<<endl;
// cout<<"xi_left: "<<chi2_smeared_xi_left<<endl;
// cout<<"logx_right: "<<chi2_smeared_logx<<endl;
// cout<<"logx_left: "<<chi2_smeared_logx_left<<endl;



// TCanvas *c1 = new TCanvas("c1","t_right");
// chi2_unfolded_vs_iter_t->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_t->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_t->SetMarkerColor(1);
// chi2_unfolded_vs_iter_t->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_t->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_t->Draw();
// TLine *line_t = new TLine(0,log10(chi2_smeared_t),40,log10(chi2_smeared_t));
// line_t->SetLineColor(kRed);
// line_t->Draw("same");
// TPad *text_line_t = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_t->SetFillColor(0);
//  // text_lixne_t->Draw();

// TCanvas *c2 = new TCanvas("c2","t_left");
// chi2_unfolded_vs_iter_t_left->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_t_left->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_t_left->SetMarkerColor(1);
// chi2_unfolded_vs_iter_t_left->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_t_left->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_t_left->Draw();
// TLine *line_t_left = new TLine(0,log10(chi2_smeared_t_left),40,log10(chi2_smeared_t_left));
// line_t_left->SetLineColor(kRed);
// line_t_left->Draw("same");
// TPad *text_line_t_left = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_t_left->SetFillColor(0);
//  // text_line_t_left->Draw();


// TCanvas *c3 = new TCanvas("c3","xi_right");
// chi2_unfolded_vs_iter_xi->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_xi->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_xi->SetMarkerColor(1);
// chi2_unfolded_vs_iter_xi->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_xi->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_xi->Draw();
// TLine *line_xi = new TLine(0,log10(chi2_smeared_xi),40,log10(chi2_smeared_xi));
// line_xi->SetLineColor(kRed);
// line_xi->Draw("same");
// TPad *text_line_xi = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_xi->SetFillColor(0);
//  // text_line_t->Draw();

// TCanvas *c4 = new TCanvas("c4","xi_left");
// chi2_unfolded_vs_iter_xi_left->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_xi_left->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_xi_left->SetMarkerColor(1);
// chi2_unfolded_vs_iter_xi_left->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_xi_left->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_xi_left->Draw();
// TLine *line_xi_left = new TLine(0,log10(chi2_smeared_xi_left),40,log10(chi2_smeared_xi_left));
// line_xi_left->SetLineColor(kRed);
// line_xi_left->Draw("same");
// TPad *text_line_xi_left = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_xi_left->SetFillColor(0);
//  // text_line_t_left->Draw();


// TCanvas *c5 = new TCanvas("c5","logx_right");
// chi2_unfolded_vs_iter_logx->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_logx->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_logx->SetMarkerColor(1);
// chi2_unfolded_vs_iter_logx->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_logx->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_logx->Draw();
// TLine *line_logx = new TLine(0,log10(chi2_smeared_logx),40,log10(chi2_smeared_logx));
// line_logx->SetLineColor(kRed);
// line_logx->Draw("same");
// TPad *text_line_logx = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_logx->SetFillColor(0);
//  // text_line_t->Draw();

// TCanvas *c6 = new TCanvas("c6","logx_left");
// chi2_unfolded_vs_iter_logx_left->SetXTitle("N_{iter}");
// chi2_unfolded_vs_iter_logx_left->SetYTitle("log_{10} #chi^{2}_{unfolded}");
// chi2_unfolded_vs_iter_logx_left->SetMarkerColor(1);
// chi2_unfolded_vs_iter_logx_left->SetMarkerSize(0.5);
// chi2_unfolded_vs_iter_logx_left->SetMarkerStyle(20);
// chi2_unfolded_vs_iter_logx_left->Draw();
// TLine *line_logx_left = new TLine(0,log10(chi2_smeared_logx_left),40,log10(chi2_smeared_logx_left));
// line_logx_left->SetLineColor(kRed);
// line_logx_left->Draw("same");
// TPad *text_line_logx_left = new TPad("","#chi^{2}_{smeared}",0.05,0.02,0.95,0.47);
//  text_line_logx_left->SetFillColor(0);
//  // text_line_t_left->Draw();



// for (int i = 2; i<40; ++i){
//     delta_chi2_unfolded_vs_iter_t->Fill(i, fabs(chi2_unfolded_t[i]-chi2_unfolded_t[i-1])/chi2_unfolded_t[i]); 
//     delta_chi2_unfolded_vs_iter_t_left->Fill(i, fabs(chi2_unfolded_t_left[i]-chi2_unfolded_t_left[i-1])/chi2_unfolded_t_left[i]); 
//     delta_chi2_unfolded_vs_iter_xi->Fill(i, fabs(chi2_unfolded_xi[i]-chi2_unfolded_xi[i-1])/chi2_unfolded_xi[i]);
//     delta_chi2_unfolded_vs_iter_xi_left->Fill(i, fabs(chi2_unfolded_xi_left[i]-chi2_unfolded_xi_left[i-1])/chi2_unfolded_xi_left[i]);
//     delta_chi2_unfolded_vs_iter_logx->Fill(i, fabs(chi2_unfolded_logx[i]-chi2_unfolded_logx[i-1])/chi2_unfolded_logx[i]);
//     delta_chi2_unfolded_vs_iter_logx_left->Fill(i, fabs(chi2_unfolded_logx_left[i]-chi2_unfolded_logx_left[i-1])/chi2_unfolded_logx_left[i]);
// }

// TLine *line_delta_chi2_2 = new TLine(0,0.02,40,0.02);
// line_delta_chi2_2->SetLineColor(kRed);
// TLine *line_delta_chi2_5 = new TLine(0,0.05,40,0.05);
// line_delta_chi2_5->SetLineColor(kRed);

// TCanvas *c7 = new TCanvas("c7","t_right");
// delta_chi2_unfolded_vs_iter_t->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_t->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_t->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_t->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_t->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_t->Draw();
// line_delta_chi2_2->Draw("same");
// line_delta_chi2_5->Draw("same");

// TCanvas *c8 = new TCanvas("c8","t_left");
// delta_chi2_unfolded_vs_iter_t_left->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_t_left->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_t_left->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_t_left->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_t_left->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_t_left->Draw();
// line_delta_chi2_2->Draw("same");
// line_delta_chi2_5->Draw("same");


// TCanvas *c9 = new TCanvas("c9","xi_right");
// delta_chi2_unfolded_vs_iter_xi->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_xi->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_xi->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_xi->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_xi->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_xi->Draw();
// line_delta_chi2_2->Draw("same");
// line_delta_chi2_5->Draw("same");

// TCanvas *c10 = new TCanvas("c10","xi_left");
// delta_chi2_unfolded_vs_iter_xi_left->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_xi_left->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_xi_left->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_xi_left->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_xi_left->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_xi_left->Draw();
// line_delta_chi2_2->Draw("same");
// line_delta_chi2_5->Draw("same");


// TCanvas *c11 = new TCanvas("c11","logx_right");
// delta_chi2_unfolded_vs_iter_logx->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_logx->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_logx->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_logx->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_logx->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_logx->Draw();
// line_delta_chi2_2->Draw("same");
// line_delta_chi2_5->Draw("same");

// TCanvas *c12 = new TCanvas("c12","logx_left");
// delta_chi2_unfolded_vs_iter_logx_left->SetXTitle("N_{iter}");
// delta_chi2_unfolded_vs_iter_logx_left->SetYTitle("#Delta#chi^{2}/#chi^{2}");
// delta_chi2_unfolded_vs_iter_logx_left->SetMarkerColor(1);
// delta_chi2_unfolded_vs_iter_logx_left->SetMarkerSize(0.5);
// delta_chi2_unfolded_vs_iter_logx_left->SetMarkerStyle(20);
// delta_chi2_unfolded_vs_iter_logx_left->Draw("hist");

// line_delta_chi2_2->Draw("same");
// line_delta_chi2_5->Draw("same");

//    outfile->Write();


// TH2F* response_t_minus_th2 = new TH2F("","",8,tbins,8,tbins);
// TH2F* response_t_plus_th2 = new TH2F("","",8,tbins,8,tbins);

// response_t_minus_th2->GetXaxis()->SetTitle("|t|_{rec} (GeV^{2})");
// response_t_minus_th2->GetYaxis()->SetTitle("|t|_{gen} (GeV^{2})");
// response_t_plus_th2->GetXaxis()->SetTitle("|t|_{rec} (GeV^{2})");
// response_t_plus_th2->GetYaxis()->SetTitle("|t|_{gen} (GeV^{2})");
// std::vector<double> sum_weights_t; sum_weights_t.clear(); sum_weights_t.resize( t_minus_th2->GetYaxis()->GetNbins() );
// std::vector<double> sum_weights_t_plus; sum_weights_t_plus.clear(); sum_weights_t_plus.resize( t_minus_th2->GetYaxis()->GetNbins() );
// for(int i=1; i <= t_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= t_minus_th2->GetYaxis()->GetNbins(); ++j){
//       sum_weights_t[j] += t_minus_th2->GetBinContent(i,j); 
//       sum_weights_t_plus[j] += t_plus_th2->GetBinContent(i,j); 
//    } 
// }

// for(int i=1; i <= t_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= t_minus_th2->GetYaxis()->GetNbins(); ++j){
//       double val = t_minus_th2->GetBinContent(i,j); 
//       double val_plus = t_plus_th2->GetBinContent(i,j); 
//       cout << "val: " << val << " - sum weights t: " << sum_weights_t[j] << endl;   
//       if( sum_weights_t[j] > 0. ) response_t_minus_th2->SetBinContent( i, j, val/sum_weights_t[j] ); 
//       if( sum_weights_t_plus[j] > 0. ) response_t_plus_th2->SetBinContent( i, j, val_plus/sum_weights_t_plus[j] ); 
//    } 
// }

//     TH2F* response_xi_minus_th2 = new TH2F("","",11,xi_bins,11,xi_bins);
//     TH2F* response_xi_plus_th2 = new TH2F("","",11,xi_bins,11,xi_bins);
// response_xi_minus_th2->GetXaxis()->SetTitle("#xi_{rec}");
// response_xi_minus_th2->GetYaxis()->SetTitle("#xi_{gen}");
// response_xi_plus_th2->GetXaxis()->SetTitle("#xi_{rec}");
// response_xi_plus_th2->GetYaxis()->SetTitle("#xi_{gen}");
//    //  TMatrix xi_right_covariance = (TMatrix)unfold_right_xi.Ereco();
//    //  TH2F* xi_right_covariance_matrix_unfold = new TH2F(xi_right_covariance);
//    //  TH2F* xi_right_covariance_matrix = new TH2F("","",11,xi_bins,11,xi_bins);
//    //  for (int i = 1; i<=xi_right_cut->GetNbinsX(); ++i){
//    //       for (int j = 1; j<=xi_right_cut->GetNbinsX(); ++j){
//    //            xi_right_covariance_matrix->SetBinContent(i,j,xi_right_covariance_matrix_unfold->GetBinContent(i,j)); 
//    //      }
//    //  }    
//    // xi_right_covariance_matrix->GetXaxis()->SetTitle("#xi");
//    // xi_right_covariance_matrix->GetYaxis()->SetTitle("#xi");
//    // xi_right_covariance_matrix->Draw();

// std::vector<double> sum_weights_xi; sum_weights_xi.clear(); sum_weights_xi.resize( xi_minus_th2->GetYaxis()->GetNbins() );
// std::vector<double> sum_weights_xi_plus; sum_weights_xi_plus.clear(); sum_weights_xi_plus.resize( xi_plus_th2->GetYaxis()->GetNbins() );
// for(int i=1; i <= xi_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= xi_minus_th2->GetYaxis()->GetNbins(); ++j){
//       sum_weights_xi[j] += xi_minus_th2->GetBinContent(i,j); 
//       sum_weights_xi_plus[j] += xi_plus_th2->GetBinContent(i,j); 
//    } 
// }

// for(int i=1; i <= xi_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= xi_minus_th2->GetYaxis()->GetNbins(); ++j){
//       double val = xi_minus_th2->GetBinContent(i,j); 
//       double val_plus = xi_plus_th2->GetBinContent(i,j); 
//       if( sum_weights_xi[j] > 0. ) response_xi_minus_th2->SetBinContent( i, j, val/sum_weights_xi[j] ); 
//       if( sum_weights_xi_plus[j] > 0. ) response_xi_plus_th2->SetBinContent( i, j, val_plus/sum_weights_xi_plus[j] ); 
//    } 
// }


// // // logx_minus_th2->Scale(1/logx_minus_th2->GetSumOfWeights());
//     TH2F* response_logx_minus_th2 = new TH2F("","",15,-4,0,15,-4,0);
// response_logx_minus_th2->GetXaxis()->SetTitle("log_{10} x_{rec}");
// response_logx_minus_th2->GetYaxis()->SetTitle("log_{10} x_{gen}");
//     TH2F* response_logx_plus_th2 = new TH2F("","",15,-4,0,15,-4,0);
// response_logx_plus_th2->GetXaxis()->SetTitle("log_{10} x_{rec}");
// response_logx_plus_th2->GetYaxis()->SetTitle("log_{10} x_{gen}");

// std::vector<double> sum_weights_x; sum_weights_x.clear(); sum_weights_x.resize( logx_minus_th2->GetYaxis()->GetNbins() );
// std::vector<double> sum_weights_x_plus; sum_weights_x_plus.clear(); sum_weights_x_plus.resize( logx_plus_th2->GetYaxis()->GetNbins() );
// for(int i=1; i <= logx_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= logx_minus_th2->GetYaxis()->GetNbins(); ++j){
//       sum_weights_x[j] += logx_minus_th2->GetBinContent(i,j); 
//       sum_weights_x_plus[j] += logx_plus_th2->GetBinContent(i,j); 
//    } 
// }

// for(int i=1; i <= logx_minus_th2->GetXaxis()->GetNbins(); ++i){
//    for(int j=1; j <= logx_minus_th2->GetYaxis()->GetNbins(); ++j){
//       double val = logx_minus_th2->GetBinContent(i,j); 
//       double val_plus = logx_plus_th2->GetBinContent(i,j); 
//       if( sum_weights_x[j] > 0. ) response_logx_minus_th2->SetBinContent( i, j, val/sum_weights_x[j] ); 
//       if( sum_weights_x_plus[j] > 0. ) response_logx_plus_th2->SetBinContent( i, j, val_plus/sum_weights_x_plus[j] ); 
//    } 
// }

//   TCanvas *c1_1 = new TCanvas("c1_1", "newpad");
//  response_t_minus_th2->Draw("colz");

//   TCanvas *c1_1_2 = new TCanvas("c1_1_2", "newpad");
//  response_t_plus_th2->Draw("colz");

//   TCanvas *c1_2 = new TCanvas("c1_2", "newpad");
// response_xi_minus_th2->Draw("colz");
//   TCanvas *c1_2_2 = new TCanvas("c1_2_2", "newpad");
// response_xi_plus_th2->Draw("colz");


//   TCanvas *c1_3 = new TCanvas("c1_3", "newpad");
// response_logx_minus_th2->Draw("colz");
//   TCanvas *c1_3_2 = new TCanvas("c1_3_2", "newpad");
// response_logx_plus_th2->Draw("colz");


//   TCanvas *c1 = new TCanvas("c1", "newpad");
// log_t_vs_log_xi_right->GetXaxis()->SetTitle("log_{10}(|t|)");
// log_t_vs_log_xi_right->GetYaxis()->SetTitle("log_{10}(#xi)");
// log_t_vs_log_xi_right->Draw("col");

//   TCanvas *c2 = new TCanvas("c2", "newpad");
// log_t_vs_log_xi_left->GetXaxis()->SetTitle("log_{10}(|t|)");
// log_t_vs_log_xi_left->GetYaxis()->SetTitle("log_{10}(#xi)");
// log_t_vs_log_xi_left->Draw("col");


//     TFile * outfile = new TFile("pythia8_xi_corr_rec.root", "RECREATE");
//     xi_rec_right_cut_nojet->Write();
//     t_rec_right_cut_nojet->Write();
//     t_rec_vs_xi_rec_right_cut->Write();
//     xi_gen_right_cut_nojet->Write();
//     t_gen_right_cut_nojet->Write();
//     t_gen_vs_xi_gen_right_cut->Write();
//     xi_cms_gen_right_cut_bin_sasha->Write();
//     xi_cms_gen_right_cut_bin_sasha_pt40->Write();
//     theta_y_vs_x_rec_nojet->Write();
//     theta_y_vs_x_gen_nojet->Write();
//     theta_x_rec_nojet->Write();
//     theta_y_rec_nojet->Write();
//     theta_x_gen_nojet->Write();
//     theta_y_gen_nojet->Write();
//     x_pos_top_right->Write();
//     y_pos_top_right->Write();
//     x_pos_bottom_right->Write();
//     y_pos_bottom_right->Write();
//     x_pos_top_left->Write();
//     y_pos_top_left->Write();
//     x_pos_bottom_left->Write();
//     y_pos_bottom_left->Write();
//     outfile->Close();
// // pt_jet1_gen_right_cut->Scale(1/luminosity, "width");
// pt_jet1_gen_right_cut->GetXaxis()->SetTitle("p_{T} (GeV)");
// pt_jet1_gen_right_cut->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} (nb/GeV)");
// pt_jet1_gen_right_cut->Draw();
    // outfile->Close();
// t_right_unfolded->Scale(1/luminosity,"width");
// t_left_unfolded->Scale(1/luminosity,"width");
// t_right_unfolded->Draw("e1");
// t_left_unfolded->Draw("e1sames");
t_right_unfolded->Scale(1/luminosity,"width");
t_left_unfolded->Scale(1/luminosity,"width");
t_right_unfolded->Draw("e1");
t_left_unfolded->Draw("e1sames");
}



