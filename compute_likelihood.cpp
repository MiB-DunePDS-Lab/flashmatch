// GENERAL
#include <algorithm>
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
#include <numeric>
#include <vector>
// UTILS





// --- HARD CODE HERE ---------------------------------------------------------
int rebinning = 4;              
bool rebin_h2 = false;

#define NOPDET 480;

const size_t n_opdet = NOPDET;

bool ignore_noreco_terms = false;

TString buildmu_output  = "./Mu_Expected_x30cm_yz50cm_fiducial_norefl.root";
TString output_filename = "./Compute_Likelihood_fiducial_norefl.root";
// ----------------------------------------------------------------------------

class Event {
  public:
    static constexpr size_t n_opdet = NOPDET;
    double true_energy;
    double vertex_coor[3];
    double mu[n_opdet];

    Event& operator=(const Event& other) {
      if (this != &other) {
        true_energy = other.true_energy;
        std::copy(std::begin(other.vertex_coor), std::end(other.vertex_coor), vertex_coor);
        std::copy(std::begin(other.mu), std::end(other.mu), mu);
      }
      return *this;
    }

  Event() : true_energy(0), vertex_coor{0, 0, 0}, mu{0} {}
  
  Event(const Double_t true_energy,
        const TTreeReaderArray<Double_t>& vertex_coor,
        const TTreeReaderArray<Double_t>& mu) {
    this->true_energy = true_energy;
    for (size_t i = 0; i < 3; ++i) {
      this->vertex_coor[i] = vertex_coor[i];
    }
    for (size_t i = 0; i < n_opdet; ++i) {
      this->mu[i] = mu[i];
    }
  }
};

class Flash {
  public:
    static constexpr size_t n_opdet = NOPDET;
    double reco_pe[n_opdet];
    double hit_t[n_opdet];
    bool hit[n_opdet];

    Flash& operator=(const Flash& other) {
      if (this != &other) {
        std::copy(std::begin(other.reco_pe), std::end(other.reco_pe), reco_pe);
        std::copy(std::begin(other.hit_t), std::end(other.hit_t), hit_t);
        std::copy(std::begin(other.hit), std::end(other.hit), hit);
      }
      return *this;
    }

    Flash() : reco_pe{0}, hit_t{0}, hit{false} {}
    Flash(const TTreeReaderArray<Double_t>& reco_pe,
          const TTreeReaderArray<Double_t>& hit_t,
          const TTreeReaderArray<Bool_t>& hit) {
      for (size_t i = 0; i < n_opdet; ++i) {
        this->reco_pe[i] = reco_pe[i];
        this->hit_t[i] = hit_t[i];
        this->hit[i] = hit[i];
      }
    }
};

std::vector<double> compute_loglikelihood_terms(Event& event,
                                                Flash& flash,
                                                TEfficiency* he_ProbHit_ExpPe,
                                                TH2D* h2_exp_reco) {
  std::vector<double> terms(event.n_opdet);
  double term = 0.;
  for(size_t idx_opdet=0; idx_opdet<event.n_opdet; idx_opdet++){
    double P_hit_mu = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(event.mu[idx_opdet]));
    if(flash.hit[idx_opdet]){
      term = P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(flash.reco_pe[idx_opdet], event.mu[idx_opdet]));
    } else {
      term = 1. - P_hit_mu;
      if (ignore_noreco_terms) term = 1.;
    }
    terms[idx_opdet] = -log(term);
  }
  return terms;
}

bool flash_matcher(Event& event_a, Flash& flash_a, Event& event_b, Flash& flash_b, TEfficiency* he_ProbHit_ExpPe, TH2D* h2_exp_reco) {
  std::vector<double> ea_fa_vector = compute_loglikelihood_terms(event_a, flash_a, he_ProbHit_ExpPe, h2_exp_reco); 
  double ea_fa = std::accumulate(ea_fa_vector.begin(), ea_fa_vector.end(), 0.0);
  std::vector<double> eb_fb_vector = compute_loglikelihood_terms(event_b, flash_b, he_ProbHit_ExpPe, h2_exp_reco);
  double eb_fb = std::accumulate(eb_fb_vector.begin(), eb_fb_vector.end(), 0.0);
  std::vector<double> ea_fb_vector = compute_loglikelihood_terms(event_a, flash_b, he_ProbHit_ExpPe, h2_exp_reco);
  double ea_fb = std::accumulate(ea_fb_vector.begin(), ea_fb_vector.end(), 0.0);
  std::vector<double> eb_fa_vector = compute_loglikelihood_terms(event_b, flash_a, he_ProbHit_ExpPe, h2_exp_reco);
  double eb_fa = std::accumulate(eb_fa_vector.begin(), eb_fa_vector.end(), 0.0);

  double min = std::min({ea_fa, eb_fb, ea_fb, eb_fa});

  if (ea_fa == min || eb_fb == min) {
    return true;
  } else {
    return false;
  }
}

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
  size_t evt_counter = 0;     size_t fake_counter = 0;        size_t mismatch_counter = 0;
  size_t zeros_term = 0;      size_t zeros_term_fake = 0;
  size_t evts_with_zeros = 0;
  int idx_event = 0;
  Long64_t entry = 0;
  Long64_t nentries = MuRecoTree->GetEntries();
  Double_t old_mu[n_opdet]; Double_t old_true_energy = 0.; Double_t old_vertex_coor[3] = {0.};

  Event fake_event, true_event;
  Flash true_flash, fake_flash;

  while(MuRecoReader.Next()){
    // if (entry == 30000) break;
    entry++;
    if (idx_event==0) {
      fake_event = Event(*true_energy, vertex_coor, mu);
      fake_flash = Flash(reco_pe, hit_t, hit);
      idx_event++;
      continue;
    }

    true_event = Event(*true_energy, vertex_coor, mu);
    true_flash = Flash(reco_pe, hit_t, hit);

    bool evt_with_zeros = false;

    std::vector<double> true_terms = compute_loglikelihood_terms(true_event, true_flash, he_ProbHit_ExpPe, h2_exp_reco);
    std::vector<double> fake_terms = compute_loglikelihood_terms(fake_event, true_flash, he_ProbHit_ExpPe, h2_exp_reco);

    double Fq = std::accumulate(true_terms.begin(), true_terms.end(), 0.0);
    double Fq_fake = std::accumulate(fake_terms.begin(), fake_terms.end(), 0.0);

    for (size_t i = 0; i < n_opdet; ++i) {
      if (true_flash.hit[i]) {
        h_term_hit->Fill(true_terms[i]);
        h_term_hit_fake->Fill(fake_terms[i]);
      } else {
        h_term_miss->Fill(true_terms[i]);
        h_term_miss_fake->Fill(fake_terms[i]);
      }
    }

    if (max_term < *std::max_element(true_terms.begin(), true_terms.end())) {
      max_term = *std::max_element(true_terms.begin(), true_terms.end());
    }
    if (min_term > *std::min_element(true_terms.begin(), true_terms.end())) {
      min_term = *std::min_element(true_terms.begin(), true_terms.end());
    }
    if (max_term_fake < *std::max_element(fake_terms.begin(), fake_terms.end())) {
      max_term_fake = *std::max_element(fake_terms.begin(), fake_terms.end());
    }
    if (min_term_fake > *std::min_element(fake_terms.begin(), fake_terms.end())) {
      min_term_fake = *std::min_element(fake_terms.begin(), fake_terms.end());
    }

    // if(evt_with_zeros) evts_with_zeros++;
    
    h_Fq->Fill(Fq);
    h_Fq_fake->Fill(Fq_fake);

    if(entry % 1000 == 0)
      std::cout << "Analyzed " << double(entry)/double(nentries)*100. << "% of the events. Fake-matches: " << double(fake_counter)/double(evt_counter)*100. << "%\r" << std::flush;
   
    if(Fq_fake < Fq) fake_counter++;
    if(!flash_matcher(true_event, true_flash, fake_event, fake_flash, he_ProbHit_ExpPe, h2_exp_reco)) mismatch_counter++;
    evt_counter++;
    fake_event = true_event;
    fake_flash = true_flash;
  } // end loop over events
  
  std::cout << "\nFake/Evts:   " <<fake_counter << "/" <<evt_counter
            << "\nFake/Evts %: " <<fake_counter/double(evt_counter)*100. 
            << "\nMismatch/Evts:   " <<mismatch_counter << "/" <<evt_counter
            << "\nMismatch/Evts %: " <<mismatch_counter/double(evt_counter)*100. 
            // << "\nZeros terms: " << zeros_term << "/" << evt_counter*n_opdet
            // << "\nZeros terms %: " << zeros_term/double(evt_counter*n_opdet)*100.
            // << "\nFake Zeros terms: " << zeros_term_fake << "/" << evt_counter*n_opdet
            // << "\nFake Zeros terms %: " << zeros_term_fake/double(evt_counter*n_opdet)*100.
            // << "\nEvents with zeros: " << evts_with_zeros << "/" << evt_counter
            // << "\nEvents with zeros %: " << evts_with_zeros/double(evt_counter)*100.
            // << "\nMin term: " << min_term
            // << "\nMax term: " << max_term
            // << "\nMin term fake: " << min_term_fake
            // << "\nMax term fake: " << max_term_fake
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
