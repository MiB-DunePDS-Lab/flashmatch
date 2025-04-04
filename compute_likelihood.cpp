// GENERAL
#include <cmath>
#include <cstddef>
#include <iostream>
// ROOT
#include <TTree.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <RtypesCore.h>
#include <TEfficiency.h>
#include <TTreeReaderValue.h>
// UTILS
#include "file_list.hpp"
#include "Utils.hpp"



// --- HARD CODE HERE ---------------------------------------------------------
int rebinning = 4;              
bool rebin_h2 = false;

const size_t n_opdet = 480;

bool ignore_noreco_terms = false;

TString buildmu_output  = "./Mu_Expected_FullTPC_norefl.root";
TString output_filename = "./Compute_Likelihood_FullTPC_norefl.root";
// ----------------------------------------------------------------------------

void compute_likelihood(){
  // --- INPUT ----------------------------------------------------------------
  TFile* exp_reco_file = TFile::Open(buildmu_output, "READ");
  TH2D* h2_exp_reco = nullptr;
  exp_reco_file->GetObject("h2_exp_reco_large", h2_exp_reco);
  TTree* MuRecoTree  = (TTree*)exp_reco_file->Get("MuRecoTree");
  TTreeReader MuRecoReader(MuRecoTree);
  TTreeReaderValue<Int_t> ifile(MuRecoReader, "ifile");
  TTreeReaderValue<Int_t> iev(MuRecoReader, "iev");
  TTreeReaderValue<Double_t> true_energy(MuRecoReader, "true_energy");
  TTreeReaderArray<Double_t> vertex_coor(MuRecoReader, "vertex_coor");
  TTreeReaderArray<Double_t> mu(MuRecoReader, "mu");
  TTreeReaderArray<Double_t> reco_pe(MuRecoReader, "reco_pe");
  TTreeReaderArray<Double_t> hit_t(MuRecoReader, "hit_t");
  TTreeReaderArray<Bool_t> hit(MuRecoReader, "hit");

  if (rebin_h2) h2_exp_reco->Rebin2D(rebinning, rebinning);

  // For each bin in the y-axis of h2_exp_reco project on the x-axis and normalize
  // the projection to have area=1
  for(int idx_y=0; idx_y<=h2_exp_reco->GetNbinsY()+1; idx_y++){
    TH1D* h1_proj = h2_exp_reco->ProjectionX("h1_proj", idx_y, idx_y);
    double h1_proj_integral = h1_proj->Integral();
    if(h1_proj_integral > 0.){
      h1_proj->Scale(1./h1_proj_integral);
      for(int idx_x=0; idx_x<=h2_exp_reco->GetNbinsX()+1; idx_x++){
        h2_exp_reco->SetBinContent(idx_x, idx_y, h1_proj->GetBinContent(idx_x));
      }
    }
  }
  // Get X and Y axis maximum values
  double reco_max = h2_exp_reco->GetXaxis()->GetXmax();
  double exp_max = h2_exp_reco->GetYaxis()->GetXmax();

  TEfficiency* he_ProbHit_ExpPe = nullptr;
  exp_reco_file->GetObject("HitProb_ExpPe", he_ProbHit_ExpPe);

  // --- DEBUG ----------------------------------------------------------------
  TH1D* h_term_hit  = new TH1D("h_term_hit ",Form("%s;%s;%s","h_term_hit ","-log(term)","count"),
                               100, 0., 50.);
  TH1D* h_term_miss = new TH1D("h_term_miss",Form("%s;%s;%s","h_term_miss","-log(term)","count"),
                               100, 0., 50.);
  TH1D* h_term_hit_fake  = new TH1D("h_term_hit_fake ",Form("%s;%s;%s","h_term_hit_fake ","-log(term)","count"),
                                    100, 0., 50.);
  TH1D* h_term_miss_fake = new TH1D("h_term_miss_fake",Form("%s;%s;%s","h_term_miss_fake","-log(term)","count"),
                                    100, 0., 50.);
  TH1D* h_Fq = new TH1D("h_Fq",Form("%s;%s;%s","h_Fq","Fq","counts"),
                        200, 0., 600.);
  TH1D* h_Fq_fake = new TH1D("h_Fq_fake",Form("%s;%s;%s","h_Fq_fake","Fq_fake","counts"),
                             200, 0., 600.);

  TFile* debug_file = TFile::Open(output_filename, "RECREATE");
  debug_file->cd();
  h2_exp_reco->Write();
 
  // --- LOOP OVER EVENTS -----------------------------------------------------
  double min_term = +1.e+60;  double min_term_fake = +1.e+60;
  double max_term = -1.e+60;  double max_term_fake = -1.e+60;
  size_t evt_counter = 0;     size_t fake_counter = 0;
  size_t zeros_term = 0;      size_t zeros_term_fake = 0;
  size_t evts_with_zeros = 0;
  int idx_event = 0;
  Long64_t entry = 0;
  Long64_t nentries = MuRecoTree->GetEntries();
  Double_t old_mu[n_opdet]; Double_t old_true_energy = 0.; Double_t old_vertex_coor[3] = {0.};
  while(MuRecoReader.Next()){
    entry++;
    if (idx_event==0) {
      old_true_energy = *true_energy;
      old_vertex_coor[0] = vertex_coor[0];
      old_vertex_coor[1] = vertex_coor[1];
      old_vertex_coor[2] = vertex_coor[2];
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        old_mu[idx_opdet] = mu[idx_opdet];
      }
      idx_event++;
      continue;
    }
    double Fq = 0.;
    double Fq_fake = 0.;
    bool evt_with_zeros = false;
    
    // --- LOOP OVER OPDETS ---------------------------------------------------
    for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
      double exp_ph = mu[idx_opdet];
      double exp_ph_fake = old_mu[idx_opdet];

      // Compute a single term of the Fq summatory ----------------------------
      double term = 0.;
      double term_fake = 0.;
      double P_hit_mu = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(exp_ph));
      double P_hit_mu_fake = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(exp_ph_fake));

      if(hit[idx_opdet]){
        term = P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(reco_pe[idx_opdet], exp_ph));
        if(term == 0.) zeros_term++;
        term_fake = P_hit_mu_fake*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(reco_pe[idx_opdet], exp_ph_fake));
        if(term_fake == 0.){
          zeros_term_fake++;
          evt_with_zeros = true;
        }
        
        term = -log(term);
        term_fake = -log(term_fake);

        h_term_hit->Fill(term);
        h_term_hit_fake->Fill(term_fake);
      } else {
        term = 1. - P_hit_mu;
        if (ignore_noreco_terms) term = 1.;
        if (term == 0.) zeros_term++;
        term_fake = 1. - P_hit_mu_fake;
        if (ignore_noreco_terms) term_fake = 1.;
        if(term_fake == 0.){
          zeros_term_fake++;
          evt_with_zeros = true;
        }

        term = -log(term);
        term_fake = -log(term_fake);

        h_term_miss->Fill(term);
        h_term_miss_fake->Fill(term_fake);
      } // single term of the Fq summatory


      if (term < min_term) min_term = term;
      if (term > max_term) max_term = term;
      if (term_fake < min_term_fake) min_term_fake = term_fake;
      if (term_fake > max_term_fake) max_term_fake = term_fake;

      Fq += term;
      Fq_fake += term_fake;
    } // end loop over opdets
    

    if(evt_with_zeros) evts_with_zeros++;
    if(Fq_fake < Fq){
      fake_counter++;
      std::cout << "\n\t" << old_true_energy << "\t" << old_vertex_coor[0] << "\t" << old_vertex_coor[1] << "\t" << old_vertex_coor[2]
                << "\n\t" << *true_energy << "\t" << vertex_coor[0] << "\t" << vertex_coor[1] << "\t" << vertex_coor[2]
                <<  std::endl;
    }
    h_Fq->Fill(Fq);
    h_Fq_fake->Fill(Fq_fake);
    evt_counter++;

    if(entry % 1000 == 0)
      std::cout << "Analyzed " << double(entry)/double(nentries)*100. << "% of the events. Mis-matches: " << double(fake_counter)/double(evt_counter)*100. << "%\r" << std::flush;
    
    old_true_energy = *true_energy; old_vertex_coor[0] = vertex_coor[0]; old_vertex_coor[1] = vertex_coor[1]; old_vertex_coor[2] = vertex_coor[2];
    for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
      old_mu[idx_opdet] = mu[idx_opdet];
    }
  } // end loop over events
  
  std::cout << "\nMismatch/Evts:   " <<fake_counter << "/" <<evt_counter
            << "\nMismatch/Evts %: " <<fake_counter/double(evt_counter)*100. 
            << "\nZeros terms: " << zeros_term << "/" << evt_counter*n_opdet
            << "\nZeros terms %: " << zeros_term/double(evt_counter*n_opdet)*100.
            << "\nFake Zeros terms: " << zeros_term_fake << "/" << evt_counter*n_opdet
            << "\nFake Zeros terms %: " << zeros_term_fake/double(evt_counter*n_opdet)*100.
            << "\nEvents with zeros: " << evts_with_zeros << "/" << evt_counter
            << "\nEvents with zeros %: " << evts_with_zeros/double(evt_counter)*100.
            << "\nMin term: " << min_term
            << "\nMax term: " << max_term
            << "\nMin term fake: " << min_term_fake
            << "\nMax term fake: " << max_term_fake
            << std::endl;  

  h_term_hit->Write();
  h_term_miss->Write();
  h_term_hit_fake->Write();
  h_term_miss_fake->Write();
  h_Fq->Write();
  h_Fq_fake->SetLineColor(kRed);
  h_Fq_fake->Write();
  debug_file->Close();
  return;
}
