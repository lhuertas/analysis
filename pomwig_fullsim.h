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

// ///ROOUNFOLD
// #include "/Users/uerj/Downloads/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
// #include "/Users/uerj/Downloads/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
// #include "/Users/uerj/Downloads/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

// Double_t beta_fit(Double_t *x, Double_t *par ){
//   Double_t result = 0;
//   // result = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0];
//   result = par[0]+par[1]*x[0]+ par[2]*pow(x[0],2) + par[3]*pow(x[0],3) +par[4]*pow(x[0],4);
//   return result;
// }
 

TFile* fullsim_file = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/pomwig_fullsim_ntuple.root","READ");
TFile* pomwig_gen_file = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/pomwig_gen_ntuple.root","READ");
TFile* pomwig_gen_histos = TFile::Open("~/Dropbox/doctorado/note/scripts/root_files/pomwig_gen.root","READ");
// string treeName = "small_tree";  
// double survival_prob = 0.1;


void pomwig_fullsim(TH2F* &pt2_xi_gen_minus_fullsim, TH2F* &pt2_xi_gen_minus_rp_fullsim, TH2F* &pt2_xi_gen_minus_pomwig, TH2F* &pt2_xi_gen_minus_rp_pomwig, 
                    TH1F* &t_rec_right_fullsim, TH1F* &t_gen_right_fullsim, TH1F* &t_gen_right_fullsim_accep, TH1F* &t_rec_right_pomwig, TH1F* &t_gen_right_pomwig, TH1F* &t_gen_right_pomwig_accep, TH1F* &xi_rec_right_fullsim,  
                    TH1F* &xi_gen_right_fullsim, TH1F* &xi_gen_right_fullsim_accep, TH1F* &xi_rec_right_pomwig, TH1F* &xi_gen_right_pomwig, TH1F* &xi_gen_right_pomwig_accep){


    string treeName = "small_tree"; 
    double norm = 0.00319507;

    float tbins[9] = {0.03, 0.07, 0.12, 0.21, 0.31, 0.42,  0.55, 0.75, 1.};
//     float tbins[9] = {0.03, 0.08, 0.13, 0.21, 0.31, 0.41,  0.55, 0.75, 1.};
    float xi_bins[12] = {-0.05, 0, 0.028, 0.048, 0.065, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2};
     

    // xi_rec_right = new TH1F("xi_rec_right_pythia8diff","",100,-0.04,0.2);
    // xi_rec_left = new TH1F("xi_rec_left_pythia8diff","",100,-0.04,0.2);
//     TH1F* xi_gen_right_fullsim = new TH1F("xi_gen_right_fullsim","",11,xi_bins);
    xi_gen_right_fullsim_accep = new TH1F("xi_gen_right_fullsim_accep","",11,xi_bins);
    TH1F* pt2_gen_right_fullsim = new TH1F("pt2_gen_right_fullsim","",25,0,1);
    TH1F* pt2_gen_right_fullsim_accep = new TH1F("pt2_gen_right_fullsim_accep","",25,0,1);
    xi_rec_right_fullsim = new TH1F("xi_rec_right_fullsim","",11,xi_bins);
    xi_gen_right_fullsim = new TH1F("xi_gen_right_fullsim","",11,xi_bins);
    TH1F*xi_gen_right_fullsim_rp = new TH1F("xi_gen_right_fullsim_rp","",11,xi_bins);
    xi_rec_right_pomwig = new TH1F("xi_rec_proton_right_pomwig","",11,xi_bins);
    xi_gen_right_pomwig = new TH1F("xi_gen_proton_right_pomwig","",11,xi_bins);
   TH1F* xi_gen_right_pomwig_rp = new TH1F("xi_gen_proton_right_pomwig_rp","",11,xi_bins);
//     TH1F* xi_gen_right_pomwig = new TH1F("xi_gen_proton_right_pomwig","",11,xi_bins);
    xi_gen_right_pomwig_accep = new TH1F("xi_gen_proton_right_pomwig_accep","",11,xi_bins);
    TH1F* pt2_gen_right_pomwig = new TH1F("pt2_gen_proton_right_pomwig","",25,0,1);
    TH1F* pt2_gen_right_pomwig_accep = new TH1F("pt2_gen_proton_right_pomwig_accep","",25,0,1);
    t_rec_right_pomwig = new TH1F("t_rec_right_pomwig","", 8, tbins);
    t_gen_right_pomwig = new TH1F("t_gen_right_pomwig","", 8, tbins);
    t_gen_right_pomwig_accep = new TH1F("t_gen_right_pomwig_accep","", 8, tbins);
    t_rec_right_fullsim = new TH1F("t_rec_right_fullsim","", 8, tbins);
    t_gen_right_fullsim = new TH1F("t_gen_right_fullsim","", 8, tbins);
    t_gen_right_fullsim_accep = new TH1F("t_gen_right_fullsim_accep","", 8, tbins);
    // t_gen_right_cut = new TH1F("t_gen_right_cut_pythia8diff","", 8, tbins);
    pt2_xi_gen_minus_fullsim = new TH2F("pt2_xi_gen_minus_fullsim","",30,0,1,30,0,0.1);
    pt2_xi_gen_minus_rp_fullsim= new TH2F("pt2_xi_gen_minus_rp_fullsim","",30,0,1,30,0,0.1);
    pt2_xi_gen_minus_pomwig = new TH2F("pt2_xi_gen_minus_pomwig","",30,0,1,30,0,0.1);
    pt2_xi_gen_minus_rp_pomwig= new TH2F("pt2_xi_gen_minus_rp_pomwig","",30,0,1,30,0,0.1);
    TH2F* rp_pos_x_vs_y_120_121 = new TH2F("rp_pos_x_vs_y_120_121","",100,-1,10,100,-40,40);
    TH2F* rp_pos_x_vs_y_120_121_bf = new TH2F("rp_pos_x_vs_y_120_121","",100,-1,10,100,-40,40);
    TH2F* rp_pos_x_vs_y_124_125 = new TH2F("rp_pos_x_vs_y_124_125","",100,-1,7,100,-40,40);
    TH2F* rp_pos_x_vs_y_124_125_bf = new TH2F("rp_pos_x_vs_y_124_125","",100,-1,7,100,-40,40);
    TH2F* rp_pos_x_vs_y_120_121_pomwig = new TH2F("rp_pos_x_vs_y_120_121_pomwig","",100,-1,10,100,-40,40);
    TH2F* rp_pos_x_vs_y_120_121_pomwig_bf = new TH2F("rp_pos_x_vs_y_120_121_pomwig","",100,-1,10,100,-40,40);
    TH2F* rp_pos_x_vs_y_124_125_pomwig = new TH2F("rp_pos_x_vs_y_124_125_pomwig","",100,-1,7,100,-40,40);
    TH2F* rp_pos_x_vs_y_124_125_pomwig_bf = new TH2F("rp_pos_x_vs_y_124_125_pomwig","",100,-1,7,100,-40,40);



    TTree* tree= (TTree*) fullsim_file->Get( treeName.c_str() );
    int nev = int(tree->GetEntriesFast());
    cout <<"The fullsim file has " << nev << " entries  " << endl;
 
    double jet1_rec_pt, jet1_rec_eta, jet1_rec_phi, jet2_rec_pt, jet2_rec_eta, jet2_rec_phi;
    double xi_rec_cms_minus, xi_rec_proton_right, xi_rec_cms_plus, xi_rec_proton_left;
    double xi_gen_cms_minus, xi_gen_proton_right, xi_gen_cms_plus, xi_gen_proton_left;
    double t_rec_proton_right, beta_rec_proton_right, t_rec_proton_left, beta_rec_proton_left;
    double t_gen_proton_right, beta_gen_proton_right, t_gen_proton_left, beta_gen_proton_left;
    double weight, rp_x_120, rp_x_121, rp_x_124, rp_x_125, rp_y_120, rp_y_121, rp_y_124, rp_y_125;
    bool rp_right_top, rp_right_bottom, rp_left_top, rp_left_bottom, good_rec_proton_right, good_gen_proton_right;
    tree->SetBranchAddress("xi_rec_proton_right",&xi_rec_proton_right);
    tree->SetBranchAddress("xi_gen_proton_right",&xi_gen_proton_right);
    tree->SetBranchAddress("t_rec_proton_right",&t_rec_proton_right);
    tree->SetBranchAddress("t_gen_proton_right",&t_gen_proton_right);
    tree->SetBranchAddress("rp_right_top",&rp_right_top);
    tree->SetBranchAddress("rp_x_120",&rp_x_120);
    tree->SetBranchAddress("rp_y_120",&rp_y_120);
    tree->SetBranchAddress("rp_x_121",&rp_x_121);
    tree->SetBranchAddress("rp_y_121",&rp_y_121);
    tree->SetBranchAddress("rp_x_124",&rp_x_124);
    tree->SetBranchAddress("rp_y_124",&rp_y_124);
    tree->SetBranchAddress("rp_x_125",&rp_x_125);
    tree->SetBranchAddress("rp_y_125",&rp_y_125);
    tree->SetBranchAddress("rp_right_bottom",&rp_right_bottom);
    tree->SetBranchAddress("good_rec_proton_right",&good_rec_proton_right);
    tree->SetBranchAddress("good_gen_proton_right",&good_gen_proton_right);


    for(int i_evt = 0; i_evt < nev; ++i_evt){
        tree->GetEntry(i_evt);

        bool proton_right_rec_sel = -xi_rec_proton_right<0.1 && -xi_rec_proton_right<0.1 && fabs(t_rec_proton_right)>0.03 && fabs(t_rec_proton_right)<1;
        bool proton_right_gen_sel = -xi_gen_proton_right<0.1 && fabs(t_gen_proton_right)>0.03 && fabs(t_gen_proton_right)<1;
        bool rp_right = rp_right_top || rp_right_bottom;
        bool fid_cut_right_top = rp_x_124>0 && rp_x_124<6 && rp_y_124 >8.4 && rp_y_124<27 ;//&& rp_x_120>0 && rp_x_120<6 && rp_y_120>8.4 && rp_y_120<27;
        bool fid_cut_right_bottom = rp_x_125>0 && rp_x_125<6 && rp_y_125<-8.4 && rp_y_125>-27;// && rp_x_121>0 && rp_x_121<6 && rp_y_121<-8.4 && rp_y_121>-27;


        // pt2_xi_gen_minus_fullsim->Fill( -t_gen_proton_right*(1 + xi_gen_proton_right), -xi_gen_proton_right, 1); 
         pt2_xi_gen_minus_fullsim->Fill( fabs(t_gen_proton_right), -xi_gen_proton_right, 1); 

        rp_pos_x_vs_y_120_121_bf->Fill(rp_x_120,rp_y_120,1.);
        rp_pos_x_vs_y_124_125_bf->Fill(rp_x_124,rp_y_124,1.);
        rp_pos_x_vs_y_120_121_bf->Fill(rp_x_121,rp_y_121,1.);
        rp_pos_x_vs_y_124_125_bf->Fill(rp_x_125,rp_y_125,1.);

        xi_gen_right_fullsim->Fill(-xi_gen_proton_right, 1);
        if (!rp_right_bottom && rp_right_top /* && proton_right_gen_sel*/){
            rp_pos_x_vs_y_120_121->Fill(rp_x_120,rp_y_120,1.);
            rp_pos_x_vs_y_124_125->Fill(rp_x_124,rp_y_124,1.);
            pt2_xi_gen_minus_rp_fullsim->Fill( fabs(t_gen_proton_right), -xi_gen_proton_right, 1);
            xi_gen_right_fullsim_rp->Fill(-xi_gen_proton_right, 1);
            if (proton_right_rec_sel){
                xi_rec_right_fullsim->Fill(-xi_rec_proton_right, 1);  
                t_rec_right_fullsim->Fill(fabs(t_rec_proton_right), 1.);
            }         
        }

        if (rp_right_bottom && !rp_right_top  /* && proton_right_gen_sel*/){
            rp_pos_x_vs_y_120_121->Fill(rp_x_121,rp_y_121,1.);
            rp_pos_x_vs_y_124_125->Fill(rp_x_125,rp_y_125,1.);
            pt2_xi_gen_minus_rp_fullsim->Fill( fabs(t_gen_proton_right), -xi_gen_proton_right, 1);       
            xi_gen_right_fullsim_rp->Fill(-xi_gen_proton_right, 1);
            if (proton_right_rec_sel){
                xi_rec_right_fullsim->Fill(-xi_rec_proton_right, 1);  
                t_rec_right_fullsim->Fill(fabs(t_rec_proton_right), 1.);
            }

        }


     }   
    

// xi_rec_right_fullsim->Divide(xi_gen_right_fullsim);
// xi_rec_right_fullsim->Draw("E1");

    TTree* tree_pomwig= (TTree*) pomwig_gen_file->Get( treeName.c_str() );
    int nev_gen = int(tree_pomwig->GetEntriesFast());
    cout <<"The pomwig gen file has " << nev_gen << " entries  " << endl;
 
    double xi_rec_proton_right_pomwig_gen;
    double xi_gen_proton_right_pomwig_gen;
    double t_rec_proton_right_pomwig_gen;
    double t_gen_proton_right_pomwig_gen;
    double jet1_gen_pt, jet1_gen_eta, jet2_gen_pt, jet2_gen_eta;
    bool rp_right_top_pomwig_gen, rp_right_bottom_pomwig_gen;
    double rp_x_120_pomwig, rp_x_121_pomwig, rp_x_124_pomwig, rp_x_125_pomwig, rp_y_120_pomwig, rp_y_121_pomwig, rp_y_124_pomwig, rp_y_125_pomwig;
    tree_pomwig->SetBranchAddress("xi_rec_proton_right",&xi_rec_proton_right_pomwig_gen);
    tree_pomwig->SetBranchAddress("xi_gen_proton_right",&xi_gen_proton_right_pomwig_gen);
    tree_pomwig->SetBranchAddress("t_rec_proton_right",&t_rec_proton_right_pomwig_gen);
    tree_pomwig->SetBranchAddress("t_gen_proton_right",&t_gen_proton_right_pomwig_gen);
    tree_pomwig->SetBranchAddress("jet1_pt",&jet1_gen_pt);
    tree_pomwig->SetBranchAddress("jet1_eta",&jet1_gen_eta);
    tree_pomwig->SetBranchAddress("jet2_pt",&jet2_gen_pt);
    tree_pomwig->SetBranchAddress("jet2_eta",&jet2_gen_eta);
    tree_pomwig->SetBranchAddress("rp_right_top",&rp_right_top_pomwig_gen);
    tree_pomwig->SetBranchAddress("rp_right_bottom",&rp_right_bottom_pomwig_gen);
    tree_pomwig->SetBranchAddress("x_pos_120",&rp_x_120_pomwig);
    tree_pomwig->SetBranchAddress("y_pos_120",&rp_y_120_pomwig);
    tree_pomwig->SetBranchAddress("x_pos_121",&rp_x_121_pomwig);
    tree_pomwig->SetBranchAddress("y_pos_121",&rp_y_121_pomwig);
    tree_pomwig->SetBranchAddress("x_pos_124",&rp_x_124_pomwig);
    tree_pomwig->SetBranchAddress("y_pos_124",&rp_y_124_pomwig);
    tree_pomwig->SetBranchAddress("x_pos_125",&rp_x_125_pomwig);
    tree_pomwig->SetBranchAddress("y_pos_125",&rp_y_125_pomwig);

    for(int i_evt = 0; i_evt < nev_gen; ++i_evt){
        tree_pomwig->GetEntry(i_evt);
        bool jet_gen_sel = jet1_gen_pt>30 && jet2_gen_pt>30 && fabs(jet1_gen_eta)<4.4 && fabs(jet2_gen_eta)<4.4;
        bool proton_right_rec_sel_pomwig_gen = fabs(t_rec_proton_right_pomwig_gen)>0.03 && fabs(t_rec_proton_right_pomwig_gen)<1. && xi_rec_proton_right_pomwig_gen<0.1;
        bool proton_right_gen_sel_pomwig_gen = fabs(t_gen_proton_right_pomwig_gen)>0.03 && fabs(t_gen_proton_right_pomwig_gen)<1. && xi_gen_proton_right_pomwig_gen<0.1;
        bool fid_cut_right_top_pomwig = rp_x_124_pomwig>0 && rp_x_124_pomwig<0.006 && rp_y_124_pomwig >0.0084 && rp_y_124_pomwig<0.027;
        bool fid_cut_right_bottom_pomwig = rp_x_125_pomwig>0 && rp_x_125_pomwig<0.006 && rp_y_125_pomwig<-0.0084 && rp_y_125_pomwig>-0.027;


        pt2_xi_gen_minus_pomwig->Fill(fabs(t_gen_proton_right_pomwig_gen), xi_gen_proton_right_pomwig_gen, 1);
        rp_pos_x_vs_y_120_121_pomwig_bf->Fill(rp_x_120_pomwig*1000,rp_y_120_pomwig*1000,1.);
        rp_pos_x_vs_y_124_125_pomwig_bf->Fill(rp_x_124_pomwig*1000,rp_y_124_pomwig*1000,1.);
            xi_gen_right_pomwig->Fill(xi_gen_proton_right_pomwig_gen, 1);

        if (!rp_right_bottom_pomwig_gen && rp_right_top_pomwig_gen  /*&& proton_right_rec_sel_pomwig_gen*/){
            rp_pos_x_vs_y_120_121_pomwig->Fill(rp_x_120_pomwig*1000,rp_y_120_pomwig*1000,1.);
            rp_pos_x_vs_y_124_125_pomwig->Fill(rp_x_124_pomwig*1000,rp_y_124_pomwig*1000,1.);
            xi_gen_right_pomwig_rp->Fill(xi_gen_proton_right_pomwig_gen, 1);
            pt2_xi_gen_minus_rp_pomwig->Fill( fabs(t_gen_proton_right_pomwig_gen), xi_gen_proton_right_pomwig_gen, 1);
            if (proton_right_rec_sel_pomwig_gen) {
                xi_rec_right_pomwig->Fill(xi_rec_proton_right_pomwig_gen,1);
                t_rec_right_pomwig->Fill(fabs(t_rec_proton_right_pomwig_gen),1);
            }
        }

        if (rp_right_bottom_pomwig_gen && !rp_right_top_pomwig_gen /*&& proton_right_rec_sel_pomwig_gen*/){
            rp_pos_x_vs_y_120_121_pomwig->Fill(rp_x_121_pomwig*1000,rp_y_121_pomwig*1000,1.);
            rp_pos_x_vs_y_124_125_pomwig->Fill(rp_x_125_pomwig*1000,rp_y_125_pomwig*1000,1.);
            pt2_xi_gen_minus_rp_pomwig->Fill( fabs(t_gen_proton_right_pomwig_gen), xi_gen_proton_right_pomwig_gen, 1);
            xi_gen_right_pomwig_rp->Fill(xi_gen_proton_right_pomwig_gen, 1);
            if (proton_right_rec_sel_pomwig_gen) {
                xi_rec_right_pomwig->Fill(xi_rec_proton_right_pomwig_gen,1);
                t_rec_right_pomwig->Fill(fabs(t_rec_proton_right_pomwig_gen),1);
            }
        }

   }

pt2_xi_gen_minus_rp_fullsim->SetXTitle("|t| GeV^{2}");
pt2_xi_gen_minus_rp_fullsim->SetYTitle("#xi");
pt2_xi_gen_minus_rp_fullsim->Divide(pt2_xi_gen_minus_fullsim);
pt2_xi_gen_minus_rp_fullsim->Draw();

pt2_xi_gen_minus_rp_pomwig->SetXTitle("|t| GeV^{2}");
pt2_xi_gen_minus_rp_pomwig->SetYTitle("#xi");
pt2_xi_gen_minus_rp_pomwig->Divide(pt2_xi_gen_minus_pomwig);
// rp_pos_x_vs_y_124_125_pomwig->Draw();
 pt2_xi_gen_minus_rp_pomwig->Divide(pt2_xi_gen_minus_rp_fullsim);
pt2_xi_gen_minus_rp_pomwig->Draw();


xi_gen_right_pomwig_rp->Divide(xi_gen_right_pomwig);
xi_gen_right_pomwig_rp->Draw();
xi_gen_right_fullsim_rp->Divide(xi_gen_right_fullsim);
xi_gen_right_fullsim_rp->Draw("histsame");


rp_pos_x_vs_y_124_125->SetXTitle("x(mm)");
rp_pos_x_vs_y_124_125->SetYTitle("y(mm)");
// rp_pos_x_vs_y_124_125->Divide(rp_pos_x_vs_y_124_125_bf);
rp_pos_x_vs_y_124_125->Draw();

// rp_pos_x_vs_y_124_125_pomwig->SetXTitle("x(mm)");
// rp_pos_x_vs_y_124_125_pomwig->SetYTitle("y(mm)");
// // rp_pos_x_vs_y_124_125_pomwig->Divide(rp_pos_x_vs_y_124_125_pomwig_bf);
// // rp_pos_x_vs_y_124_125_pomwig->Divide(rp_pos_x_vs_y_124_125);
// rp_pos_x_vs_y_124_125_pomwig->Draw();

}



