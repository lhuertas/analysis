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
// #include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldResponse.h"
// #include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldBayes.h"
// #include "/Users/lina/Dropbox/doctorado/RooUnfold/src/RooUnfoldBinByBin.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
// #include "/home/lina/Downloads/RooUnfold-1.1.1/src/RooUnfoldBayes.h"



 

// TFile* pomwig_pom = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pomwig_minus_ntuple_pt30_jec.root","READ");
// TFile* pomwig_regg = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pomwig_reggeon_minus_ntuple_newtheta_jec.root","READ");
// TFile* pomwig_pom_plus = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pomwig_plus_ntuple_pt30_jec.root","READ");
// TFile* pomwig_regg_plus = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pomwig_reggeon_plus_ntuple_newtheta_jec.root","READ");
// TFile* pomwig_pom = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/pomwig_minus_ntuple_pt20.root","READ");
TFile* pomwig_pom = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_minus_ntuple_beam_smearing.root","READ");
TFile* pomwig_regg = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_reggeon_minus_ntuple_beam_smearing.root","READ");
TFile* pomwig_pom_plus = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_plus_ntuple_beam_smearing.root","READ");
TFile* pomwig_regg_plus = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/beam_smearing/pomwig_reggeon_plus_ntuple_beam_smearing.root","READ");
// string treeName = "small_tree";  
map<string,TH1F*> reggeon_histos;
 TF1* func = new TF1("func", beta_fit, 0., 1., 3);

/////////// reweight slope
double t_slope_right = 5.35443;
double t_slope_left = 5.35443;

////// normalization factors
   double norm_pom = (0.049154*4.04e+06)/1.90724e+06;//1.11737;
    double norm_regg = 0.458476;
    double norm_pom_plus = (0.049154*4.04e+06)/1.83134e+06;//1.04824;
    double norm_regg_plus = 0.303177;




void pomwig(TH1F* &xi_cms_minus_totem_rec_right, TH1F* &xi_cms_minus_totem_rec_left, TH1F* &xi_cms_minus_totem_rec_right_bin, TH1F* &xi_cms_minus_totem_rec_left_bin, TH1F* &xi_rec_right, TH1F* &xi_rec_left, 
            TH1F* &xi_rec_right_cut, TH1F* &xi_rec_right_cut_sasha,TH1F* &xi_rec_left_cut_sasha,TH1F* &xi_rec_left_cut, TH1F* &xi_gen_right, TH1F* &xi_gen_left,TH1F* &t_rec_right_cut, TH1F* &t_gen_right_cut, 
            TH1F* &t_rec_left_cut, TH1F* &t_gen_left_cut, TH1F* &beta_right_cut, TH1F* &beta_left_cut, TH1F* &xi_cms_rec_right, TH1F* &xi_cms_gen_right,TH1F* &xi_cms_rec_left, TH1F* &xi_cms_gen_left, 
            TH1F* &xi_cms_rec_right_sasha, TH1F* &xi_cms_gen_right_sasha,TH1F* &xi_cms_rec_left_sasha, TH1F* &xi_cms_gen_left_sasha, 
            TH1F* &log_x_rec_right_cut, 
            TH1F* &log_x_gen_right_cut, TH1F* &log_x_rec_left_cut,  TH1F* &log_x_gen_left_cut, TH1F* &pt_jet1_right_cut, TH1F* &pt_jet1_left_cut, TH1F* &pt_jet2_right_cut, TH1F* &pt_jet2_left_cut, 
            TH1F* &eta_jet1_right_cut, TH1F* &eta_jet1_left_cut, TH1F* &eta_jet2_right_cut, TH1F* &eta_jet2_left_cut,  TH1F* &delta_eta_jets_right_cut, TH1F* &delta_eta_jets_left_cut, TH1F* &delta_phi_jets_right_cut, TH1F* &delta_phi_jets_left_cut,
            TH1F* &th_x_rec_right, TH1F* &th_y_rec_right, TH1F* &th_x_rec_left, TH1F* &th_y_rec_left,
            TH1F* &th_x_gen_right, TH1F* &th_y_gen_right, TH1F* &th_x_gen_left, TH1F* &th_y_gen_left,
            TH1F* &mass_jj_rec_right, TH1F* &mass_jj_rec_left, TH1F* &mass_x_rec_right, TH1F* &mass_x_rec_left, TH1F* &r_jj_rec_right, TH1F* &r_jj_rec_left,
            TH1F* &t_right_unfolded,  TH1F* &t_left_unfolded, TH1F* &xi_right_unfolded,  TH1F* &xi_left_unfolded,
            TH1F* &logx_right_unfolded,  TH1F* &logx_left_unfolded, TH1F* &t_right_unfolded_noeff,  TH1F* &t_left_unfolded_noeff, TH1F* &xi_right_unfolded_noeff,  TH1F* &xi_left_unfolded_noeff,
            TH1F* &logx_right_unfolded_noeff,  TH1F* &logx_left_unfolded_noeff, TH1F* &logx_right_noeff_unfolded_noeff,  TH1F* &logx_left_noeff_unfolded_noeff,  TH1F* &t_right_unfolded_binbybin,  
            TH1F* &t_left_unfolded_binbybin, TH1F* &xi_right_unfolded_binbybin,  TH1F* &xi_left_unfolded_binbybin,
            TH1F* &logx_right_unfolded_binbybin,  TH1F* &logx_left_unfolded_binbybin,  TH1F* &t_right_unfolded_up,  TH1F* &t_right_unfolded_dw, TH1F* &t_left_unfolded_up,  TH1F* &t_left_unfolded_dw,  
            TH1F* &xi_right_unfolded_up,  TH1F* &xi_right_unfolded_dw, TH1F* &xi_left_unfolded_up,  TH1F* &xi_left_unfolded_dw,
            TH1F* &log_x_right_unfolded_up,  TH1F* &log_x_right_unfolded_dw, TH1F* &log_x_left_unfolded_up,  TH1F* &log_x_left_unfolded_dw, 
            TH1F* &t_right_unfolded_up_pf,  TH1F* &t_right_unfolded_dw_pf, TH1F* &t_left_unfolded_up_pf,  TH1F* &t_left_unfolded_dw_pf, TH1F* &xi_right_unfolded_up_pf,  TH1F* &xi_right_unfolded_dw_pf,
            TH1F* &xi_left_unfolded_up_pf,  TH1F* &xi_left_unfolded_dw_pf, TH1F* &log_x_right_unfolded_up_pf,  TH1F* &log_x_right_unfolded_dw_pf, TH1F* &log_x_left_unfolded_up_pf,  TH1F* &log_x_left_unfolded_dw_pf, 
            TH1F* &t_right_unfolded_nobackgsub, TH1F* &t_left_unfolded_nobackgsub, TH1F* &xi_right_unfolded_nobackgsub, TH1F* &xi_left_unfolded_nobackgsub, TH1F* &logx_right_unfolded_nobackgsub, TH1F* &logx_left_unfolded_nobackgsub,
            TH1F* &t_right_pythia_unfolded, TH1F* &t_left_pythia_unfolded, TH1F* &xi_right_pythia_unfolded, TH1F* &xi_left_pythia_unfolded, TH1F* &logx_right_pythia_unfolded, TH1F* &logx_left_pythia_unfolded, 
            bool reweight = false, bool reweight_slope = false, bool unc_rp = false, bool unc_gauss = false,  bool sys_xi_up = false, bool sys_xi_dw = false, bool rp_x_min = false, bool rp_x_max = false,
            bool single_vertex = false, bool iter_variation = false);


