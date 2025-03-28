
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <filesystem>  // C++17

#include <TTree.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include "RtypesCore.h"
#include "TEfficiency.h"
#include "TTreeReaderValue.h"
#include "file_list.hpp"
#include "Utils.hpp"



// --- HARD CODE HERE ----------------
int rebinning = 2;

const size_t n_opdet = 480; // 480

TString buildmu_output = "./Mu_Expected_seed2000_ana_hit_thr.root";
// -----------------------------------

void compute_likelihood(){
  // --- INPUT ----------------------------------------------------------------
  TFile* exp_reco_file = TFile::Open(buildmu_output, "READ");
  TH2D* h2_exp_reco = nullptr;
  exp_reco_file->GetObject("h2_exp_reco_large", h2_exp_reco);
  TTree* MuRecoTree  = (TTree*)exp_reco_file->Get("MuRecoTree");
  TTreeReader MuRecoReader(MuRecoTree);
  TTreeReaderValue<Int_t> ifile(MuRecoReader, "ifile");
  TTreeReaderValue<Int_t> iev(MuRecoReader, "iev");
  TTreeReaderArray<Double_t> mu(MuRecoReader, "mu");
  TTreeReaderArray<Double_t> reco_pe(MuRecoReader, "reco_pe");
  TTreeReaderArray<Double_t> hit_t(MuRecoReader, "hit_t");
  TTreeReaderArray<Bool_t> hit(MuRecoReader, "hit");

  h2_exp_reco->Rebin2D(rebinning, rebinning);

  // for each bin in the y-axis of h2_exp_reco project on the x-axis and normalize
  // the projection to have area=1
  for(int idx_y=0; idx_y<=h2_exp_reco->GetNbinsY(); idx_y++){
    TH1D* h1_proj = h2_exp_reco->ProjectionX("h1_proj", idx_y, idx_y);
    double h1_proj_integral = h1_proj->Integral();
    if(h1_proj_integral > 0.){
      h1_proj->Scale(1./h1_proj->Integral());
      for(int idx_x=0; idx_x<=h2_exp_reco->GetNbinsX(); idx_x++){
        h2_exp_reco->SetBinContent(idx_x, idx_y, h1_proj->GetBinContent(idx_x));
      }
    }
  }

  // TFile* debug_file = TFile::Open("./debug_xy_buildmu.root", "RECREATE");
  // debug_file->cd();
  // h2_exp_reco->Write();
  // debug_file->Close();
  // return;

  // TFile* exp_reco_file2 = TFile::Open("./out_buildmu.root", "READ");
  TEfficiency* he_ProbHit_ExpPe = nullptr;
  exp_reco_file->GetObject("HitProb_ExpPe", he_ProbHit_ExpPe);

 
  // --- LOOP OVER EVENTS ------------------------------------------------------
  size_t fake_counter = 0;
  size_t evt_counter = 0;
  size_t zeros_term = 0;
  size_t zeros_term_fake = 0;
  int idx_event = 0;
  Long64_t entry = 0;
  Long64_t nentries = MuRecoTree->GetEntries();
  Double_t old_mu[n_opdet];
  while(MuRecoReader.Next()){
    entry++;
    if (idx_event==0) {
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        old_mu[idx_opdet] = mu[idx_opdet];
      }
      idx_event++;
      continue;
    }
    double Fq = 0.;
    double Fq_fake = 0.;
    
    for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
      double exp_ph = mu[idx_opdet];
      double exp_ph_fake = old_mu[idx_opdet];

      // Compute a single term of the Fq summatory
      double term = 0.;
      double term_fake = 0.;
      // double P_hit_mu = 0.5; 
      // double P_hit_mu_fake = 0.5;
      double P_hit_mu = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(exp_ph));
      double P_hit_mu_fake = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(exp_ph_fake));

      if(hit[idx_opdet]){
        term = P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(reco_pe[idx_opdet], exp_ph));
        term_fake = P_hit_mu_fake*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(reco_pe[idx_opdet], exp_ph_fake));
      } else {
        term = 1.;
        term_fake = 1.;
      }

      if(term == 0.) zeros_term++;
      if(term_fake == 0.) zeros_term_fake++;

      Fq += -log(term);
      Fq_fake += -log(term_fake);
    } // end loop over opdets
    
    for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
      old_mu[idx_opdet] = mu[idx_opdet];
    }

    if(Fq_fake < Fq) fake_counter++;
    evt_counter++;

    if(entry % 1000 == 0)
      std::cout << "Analyzed " << double(entry)/double(nentries)*100. << "% of the events. Mis-matches: " << double(fake_counter)/double(evt_counter)*100. << "%\r" << std::flush;
  } // end loop over events
  
  std::cout << "\nMismatch/Evts:   " <<fake_counter << "/" <<evt_counter
            << "\nMismatch/Evts %: " <<fake_counter/double(evt_counter)*100. 
            << "\nZeros terms: " << zeros_term << "/" << evt_counter*n_opdet
            << "\nZeros terms %: " << zeros_term/double(evt_counter*n_opdet)*100.
            << "\nFake Zeros terms: " << zeros_term_fake << "/" << evt_counter*n_opdet
            << "\nFake Zeros terms %: " << zeros_term_fake/double(evt_counter*n_opdet)*100.
            << std::endl;  

}
