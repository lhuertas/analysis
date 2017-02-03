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

TFile* jet_up_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/data_SDdijet_CMSTOTEM_jet_up.root","READ");
TFile* jet_dw_file = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/data_SDdijet_CMSTOTEM_jet_dw.root","READ");
TFile* jet_up_file_rereco = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_jet_up_jecwinter.root","READ");
TFile* jet_dw_file_rereco = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_jet_dw_jecwinter.root","READ");
TFile* pf_up_file_rereco = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_pf_up_jecwinter.root","READ");
TFile* pf_dw_file_rereco = TFile::Open("~/cernbox/doctorado/note/scripts/root_files/rereco/data_ntuple_pf_dw_jecwinter.root","READ");


void unc(TH1F* &t_right_cut, TH1F* &t_left_cut, TH1F* &xi_right_cut, TH1F* &xi_left_cut, TH1F* &log_x_minus, TH1F* &log_x_plus, TH1F* &log_x_right_cut, TH1F* &log_x_left_cut, bool jet = false, bool up = false, bool rereco = false);

Double_t fFermiLike(Double_t *x, Double_t *par) {
  Double_t result = 0;
  result = 1.0/(TMath::Exp((par[0]-TMath::Sqrt(x[0]))/par[1]) + 1);
  return result;
} 

//     float tbins[9] = {0.03, 0.08, 0.13, 0.21, 0.31, 0.41,  0.55, 0.75, 1.};


//     t_right_cut = new TH1F("pt_jet1_right_cut_pomwig","", 8, tbins);
// //     pt_jet2_right_cut = new TH1F("pt_jet2_right_cut_pomwig","", 15, 0, 200);
// //     eta_jet1_right_cut = new TH1F("eta_jet1_right_cut_pomwig","", 20, -5.2, 5.2);
// //     eta_jet2_right_cut = new TH1F("eta_jet2_right_cut_pomwig","", 20, -5.2, 5.2);
// //     pt_jet1_left_cut = new TH1F("pt_jet1_left_cut_pomwig","", 15, 0, 200);
// //     pt_jet2_left_cut = new TH1F("pt_jet2_left_cut_pomwig","", 15, 0, 200);
// //     eta_jet1_left_cut = new TH1F("eta_jet1_left_cut_pomwig","", 20, -5.2, 5.2);
// //     eta_jet2_left_cut = new TH1F("eta_jet2_left_cut_pomwig","", 20, -5.2, 5.2);
// //     TH2F* rp_pos_left = new TH2F("rp_pos_left","",100,-0.1,0.1,100,-0.4,0.4);
// //     TH2F* rp_pos_rigth = new TH2F("rp_pos_rigth","",100,-0.1,0.1,100,-0.4,0.4);


//     TTree* tree_jetup = (TTree*) jet_up->Get( treeName.c_str() );
//     int nev_jetup = int(tree_jetup->GetEntriesFast());
 
//     double jet1_pt, jet1_eta, jet1_phi, jet2_pt, jet2_eta, jet2_phi;
//     double xi_cms_minus_data, xi_proton_right_data, xi_cms_plus_data, xi_proton_left_data, x_right_data, x_left_data;
//     double t_proton_right_data, t_proton_left_data, beta_proton_right_data, beta_proton_left_data, eff_trigger;
//     bool valid_proton_right_data, valid_proton_left_data;
//     bool rp_right_top_data, rp_right_bottom_data, rp_left_top_data, rp_left_bottom_data;
//     double x_pos_024, x_pos_025, x_pos_124, x_pos_125, y_pos_024, y_pos_025, y_pos_124, y_pos_125, vtx_x, vtx_y, vtx_z;
//     tree_jetup->SetBranchAddress("xi_cms_minus",&xi_cms_minus_data);
//     tree_jetup->SetBranchAddress("xi_cms_plus",&xi_cms_plus_data);
//     tree_jetup->SetBranchAddress("xi_totem_right",&xi_proton_right_data);
//     tree_jetup->SetBranchAddress("xi_totem_left",&xi_proton_left_data);
//     tree_jetup->SetBranchAddress("t_totem_right",&t_proton_right_data);
//     tree_jetup->SetBranchAddress("t_totem_left",&t_proton_left_data);
//     tree_jetup->SetBranchAddress("beta_proton_right",&beta_proton_right_data);
//     tree_jetup->SetBranchAddress("beta_proton_left",&beta_proton_left_data);
//     tree_jetup->SetBranchAddress("x_right",&x_right_data);
//     tree_jetup->SetBranchAddress("x_left",&x_left_data);
//     tree_jetup->SetBranchAddress("valid_proton_right",&valid_proton_right_data);
//     tree_jetup->SetBranchAddress("valid_proton_left",&valid_proton_left_data);
//     tree_jetup->SetBranchAddress("rp_right_top",&rp_right_top_data);
//     tree_jetup->SetBranchAddress("rp_right_bottom",&rp_right_bottom_data);
//     tree_jetup->SetBranchAddress("rp_left_top",&rp_left_top_data);
//     tree_jetup->SetBranchAddress("rp_left_bottom",&rp_left_bottom_data);
//     tree_jetup->SetBranchAddress("jet1_pt",&jet1_pt);
//     tree_jetup->SetBranchAddress("jet1_eta",&jet1_eta);
//     tree_jetup->SetBranchAddress("jet1_phi",&jet1_phi);
//     tree_jetup->SetBranchAddress("jet2_pt",&jet2_pt);
//     tree_jetup->SetBranchAddress("jet2_eta",&jet2_eta);
//     tree_jetup->SetBranchAddress("jet2_phi",&jet2_phi);
//     tree_jetup->SetBranchAddress("x_pos_024",&x_pos_024);
//     tree_jetup->SetBranchAddress("y_pos_024",&y_pos_024);
//     tree_jetup->SetBranchAddress("x_pos_025",&x_pos_025);
//     tree_jetup->SetBranchAddress("y_pos_025",&y_pos_025);
//     tree_jetup->SetBranchAddress("x_pos_124",&x_pos_124);
//     tree_jetup->SetBranchAddress("y_pos_124",&y_pos_124);
//     tree_jetup->SetBranchAddress("x_pos_125",&x_pos_125);
//     tree_jetup->SetBranchAddress("y_pos_125",&y_pos_125);
//     tree_jetup->SetBranchAddress("vtx_x",&vtx_x);
//     tree_jetup->SetBranchAddress("vtx_y",&vtx_y);
//     tree_jetup->SetBranchAddress("vtx_z",&vtx_z);


//     for(int i_evt = 0; i_evt < nev_jetup; ++i_evt){
//         tree_jetup->GetEntry(i_evt);
 
//         bool jet_sel = jet1_pt>pt_threshold && jet2_pt>pt_threshold && fabs(jet1_eta)<4.4 && fabs(jet2_eta)<4.4;
//         bool proton_sel = xi_proton_right_data<=0.1 && fabs(t_proton_right_data)>=0.03 && fabs(t_proton_right_data)<=1.;
//         bool rp_right = rp_right_bottom_data || rp_right_top_data;
//   	    eff_trigger = func40->Eval(jet2_pt);
           
//         if (jet_sel && proton_sel && rp_right){
//             if (xi_cms_minus_data - xi_proton_right_data<0){
//                 // xi_rec_right_cut->Fill(xi_rec_proton_right_pom, event_weight_pom_minus);
//                 // xi_rec_right_cut_sasha->Fill(xi_rec_proton_right_pom, event_weight_pom_minus);
//                 t_right_cut->Fill(fabs(t_proton_right_data), 1/(eff_trigger*0.94));
//                 // beta_right_cut->Fill(beta_rec_proton_right_pom, event_weight_pom_minus);
//                 // t_rec_right_response->Fill(fabs(t_rec_proton_right_pom), event_weight_pom_minus);
//                 // xi_rec_right_response->Fill(xi_rec_proton_right_pom, event_weight_pom_minus);
//                 // log_x_rec_right_cut->Fill(log10(x_rec_right_pom), event_weight_pom_minus);
// // 		pt_jet1_right_cut->Fill(jet1_pt, 1.);
// 		// pt_jet2_right_cut->Fill(jet2_rec_pt_pom, event_weight_pom_minus);
// 		// eta_jet1_right_cut->Fill(jet1_rec_eta_pom, event_weight_pom_minus);
// 		// eta_jet2_right_cut->Fill(jet2_rec_eta_pom, event_weight_pom_minus);
//   //               rp_pos_rigth->Fill(rp_xpos_124_pom, rp_ypos_124_pom, 1.);
//   //               if (fabs(t_rec_proton_right_pom) && !fabs(t_gen_proton_right_pom)) t_minus_response.Fake(fabs(t_rec_proton_right_pom), event_weight_pom_minus);
//   //               if (jet_gen_sel_pom && proton_gen_sel_pom){
//   //                  t_minus_response.Fill(fabs(t_rec_proton_right_pom), fabs(t_gen_proton_right_pom), event_weight_pom_minus);
//   //                  xi_minus_response.Fill(xi_rec_proton_right_pom, xi_gen_proton_right_pom, event_weight_pom_minus);
//   //                  response_xi_minus->Fill(xi_rec_proton_right_pom, xi_gen_proton_right_pom, event_weight_pom_minus);
//   //                  response_t_minus->Fill(fabs(t_rec_proton_right_pom),fabs(t_gen_proton_right_pom), event_weight_pom_minus);
// 		// }  
//             }

//         }

//     }
    
    

// }    