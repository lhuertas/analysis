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


// Double_t beta_fit(Double_t *x, Double_t *par ){
//   Double_t result = 0;
//   // result = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0];
//   result = par[0]+par[1]*x[0]+ par[2]*pow(x[0],2) + par[3]*pow(x[0],3) +par[4]*pow(x[0],4);
//   return result;
// }
 

TFile* pythia8_nondiffmc = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/pythia8_nondiff_QCD_pt_15_3000_ntuple.root","READ");
// string treeName = "small_tree";  
// double survival_prob = 0.1;
double scale_pythia8 = 0.889191;//0.815821;


void pythia8_nondiff(TH1F* &pt_jet1, TH1F* &pt_jet2, TH1F* &eta_jet1, TH1F* &eta_jet2, TH1F* &log_x_minus_rec,  TH1F* &log_x_minus_gen, TH1F* &log_x_plus_rec, TH1F* &log_x_plus_gen, 
                     TH1F* &logx_minus_unfolded, TH1F* &logx_plus_unfolded){



////////////// difff mc
    double norm = 57.5102;
    string treeName = "small_tree"; 


    pt_jet1 = new TH1F("pt_jet1_pythia8","", 15, 0, 200);
    pt_jet2 = new TH1F("pt_jet2_pythia8","", 15, 0, 200);
    eta_jet1 = new TH1F("eta_jet1_pythia8","", 20, -5.2, 5.2);
    eta_jet2 = new TH1F("eta_jet2_pythia8","", 20, -5.2, 5.2);
    log_x_minus_rec = new TH1F("log_x_minus_rec_pythia8","",15, -4, 0);
    log_x_minus_gen = new TH1F("log_x_minus_gen_pythia8","",15, -4, 0);
    log_x_plus_rec = new TH1F("log_x_plus_rec_pythia8","",15, -4, 0);
    log_x_plus_gen = new TH1F("log_x_plus_gen_pythia8","",15, -4, 0);

    RooUnfoldResponse logx_minus_response (log_x_minus_rec, log_x_minus_gen, "logx_minus_unfolded", "logx_minus_unfolded");
    RooUnfoldResponse logx_plus_response (log_x_plus_rec, log_x_plus_gen, "logx_plus_unfolded", "logx_plus_unfolded");

    TTree* tree= (TTree*) pythia8_nondiffmc->Get( treeName.c_str() );
    int nev = int(tree->GetEntriesFast());
    cout <<"The pythia8 non_diffractive file has " << nev << " entries  " << endl;
 
    double jet1_rec_pt, jet1_rec_eta, jet1_rec_phi, jet2_rec_pt, jet2_rec_eta, jet2_rec_phi;
    double jet1_gen_pt, jet1_gen_eta, jet1_gen_phi, jet2_gen_pt, jet2_gen_eta, jet2_gen_phi;
    double x_rec_right, x_rec_left, x_gen_right, x_gen_left;
    double weight;
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
    tree->SetBranchAddress("x_rec_left",&x_rec_left);
    tree->SetBranchAddress("x_gen_right",&x_gen_right);
    tree->SetBranchAddress("x_gen_left",&x_gen_left);
    tree->SetBranchAddress("weight",&weight);


    for(int i_evt = 0; i_evt < nev; ++i_evt){
        tree->GetEntry(i_evt);

        bool jet_rec_sel = jet1_rec_pt>pt_theshold && jet2_rec_pt>pt_theshold && fabs(jet1_rec_eta)<4.4 && fabs(jet2_rec_eta)<4.4;
        bool jet_gen_sel = jet1_gen_pt>pt_theshold && jet2_gen_pt>pt_theshold && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4;
           
        if (jet_rec_sel){
            pt_jet1->Fill(jet1_rec_pt, weight);
            pt_jet2->Fill(jet2_rec_pt, weight);
            eta_jet1->Fill(jet1_rec_eta, weight);
            eta_jet2->Fill(jet2_rec_eta, weight);
            log_x_minus_rec->Fill(log10(x_rec_right),weight);
            log_x_plus_rec->Fill(log10(x_rec_left),weight);
            if (!jet_gen_sel){
                logx_minus_response.Fake(log10(x_rec_right), weight*norm*scale_pythia8);
                logx_plus_response.Fake(log10(x_rec_left), weight*norm*scale_pythia8);
            }
            if (jet_gen_sel){
                logx_minus_response.Fill(log10(x_rec_right), log10(x_gen_right), weight*norm*scale_pythia8);
                logx_plus_response.Fill(log10(x_rec_left), log10(x_gen_left), weight*norm*scale_pythia8);
            }
        } 

        if (jet_gen_sel){
            log_x_minus_gen->Fill(log10(x_gen_right),weight);
            log_x_plus_gen->Fill(log10(x_gen_left),weight);           
            if (!jet_rec_sel){
                logx_minus_response.Miss(log10(x_gen_right), weight*norm*scale_pythia8);
                logx_plus_response.Miss(log10(x_gen_left), weight*norm*scale_pythia8);
            }
        }

    }


    pt_jet1->Scale(norm*scale_pythia8);
    pt_jet2->Scale(norm*scale_pythia8);
    eta_jet1->Scale(norm*scale_pythia8);
    eta_jet2->Scale(norm*scale_pythia8);
    log_x_minus_rec->Scale(norm*scale_pythia8);
    log_x_minus_gen->Scale(norm*scale_pythia8);
    log_x_plus_rec->Scale(norm*scale_pythia8);
    log_x_plus_gen->Scale(norm*scale_pythia8);


    TH1F* log_x_minus = new TH1F("log_x_minus","",15, -4, 0);
    TH1F* log_x_plus = new TH1F("log_x_plus","",15, -4, 0);
    TH1F* h1;
    TH1F* h2; 
    TH1F* h3; TH1F* h4; TH1F* h5; TH1F* h6; TH1F* h7; TH1F* h8; TH1F* h9; TH1F* h10; TH1F* h11; TH1F* h12; TH1F* h13; TH1F* h14; TH1F* h15; TH1F* h16; TH1F* h17; TH1F* h18; TH1F* h19; TH1F* h20; 
    TH1F* h21; TH1F* h22; TH1F* h23; TH1F* h24; TH1F* h25; TH1F* h26; TH1F* h27; TH1F* h28; TH1F* h29; TH1F* h30; TH1F* h31; TH1F* h32; TH1F* h33; TH1F* h34; TH1F* h35; TH1F* h36; TH1F* h37; TH1F* h38;
     TH1F* h39; TH1F* h40; TH1F* h41; 
    TH1F* h42; TH1F* h43; TH1F* h44; TH1F* h45; TH1F* h46; TH1F* h47; TH1F* h48; TH1F* h49; TH1F* h50; TH1F* h51; TH1F* h52; TH1F* h53; TH1F* h54; TH1F* h55; TH1F* h56; TH1F* h57; TH1F* h58;
    TH1F* th_x_right; TH1F* th_y_right; TH1F* th_x_left; TH1F* th_y_left;

    data(h1, h2, h3, h4, h5, 
          h6, h7, h8, h9, h10, h11, h12, h13,
          h14, h15, h16, h17, h17, h18, h19, h20, 
          h21, h22, h23, h23, h24, h25, h26, h27,
          h28, h29, h30, h31, h32, h33, h34, h35, 
          h36, h37, h38, h39, h40, h41, h42,
          log_x_minus, h58, h43, h44, h45, log_x_plus, h57, h46, h47, h48, 
          h51, h52, h53, h54, h55, h56, th_x_right, th_y_right, th_x_left, th_y_left, false, true, false, false, false, false, false);


    RooUnfoldBayes unfold_minus_logx (&logx_minus_response, log_x_minus, 4);
    logx_minus_unfolded = (TH1F*)unfold_minus_logx.Hreco();

    RooUnfoldBayes unfold_plus_logx (&logx_plus_response, log_x_plus, 4);
    logx_plus_unfolded = (TH1F*)unfold_plus_logx.Hreco();


}



