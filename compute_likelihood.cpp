
#define NOPDET 480;

// GENERAL
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>
// ROOT
#include <TTree.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <RtypesCore.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TTreeReaderValue.h>
#include <Math/MinimizerOptions.h>
#include "TGraphErrors.h"
// UTILS
#include "Utils.hpp"


std::vector<double> compute_loglikelihood_terms(Event& event,
                                                Flash& flash,
                                                TF1*  f_ProbHit_ExpPe,
                                                TF1* f_lognormal,
                                                TF1*  f_logms_trend,
                                                TF1*  f_sigmas_trend,
                                                TGraphErrors* g_logms,
                                                TGraphErrors* g_sigmas,
                                                const double& trend_thr){
  const size_t n_opdet = event.n_opdet;
  std::vector<double> terms(n_opdet);
  double term = 0., logm = 0., sigma = 0., P_hit_mu = 0., reco_pe = 0.;
  for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
    P_hit_mu = f_ProbHit_ExpPe->Eval(event.mu[idx_opdet]);
    if(flash.hit[idx_opdet]){
      reco_pe = flash.reco_pe[idx_opdet];
      if (reco_pe > trend_thr){
        f_lognormal->SetParameters(f_logms_trend->Eval(reco_pe), f_sigmas_trend->Eval(reco_pe));
      } else {
        logm  = g_logms->Eval(reco_pe);
        sigma = g_sigmas->Eval(reco_pe);
        f_lognormal->SetParameters(logm, sigma);
      }
      term = P_hit_mu*f_lognormal->Eval(reco_pe);
    } else {
      term = 1. - P_hit_mu;
    }
    terms[idx_opdet] = -log(term);
  }
  return terms;
}



// --- HARD CODE HERE ---------------------------------------------------------
int rebinning = 2;              
bool rebin_h2 = true;


const size_t n_opdet = NOPDET;

bool ignore_noreco_terms = false;

bool use_fit_efficiency = true; // Use fit of TEfficiency "he_hit_prob" or the TEfficiency itself
bool use_fit_lognormal = true; // Use fit of lognormal estrapolation

double trend_thr = 25.;

TString buildmu_output  = "./Mu_Expected_x30cm_yz50cm_fiducial_norefl.root";
TString output_filename = "./Compute_Likelihood_fiducial_norefl.root";
// ----------------------------------------------------------------------------



void compute_likelihood(){
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
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
  
  
  TF1* f_lognormal = new TF1("f_lognormal", "ROOT::Math::lognormal_pdf(x,[0],[1])", 0, 2000);
  f_lognormal->SetParNames("log(m)", "#sigma");
  f_lognormal->SetNpx(3000);
  std::vector<double> logms, sigmas, reco_pes;
  std::vector<double> err_logms, err_sigmas, err_reco_pes;

  
  // For each bin in the y-axis of h2_exp_reco project on the x-axis and normalize
  // the projection to have area=1
  if (rebin_h2) h2_exp_reco->Rebin2D(rebinning, rebinning);
  for(int idx_y=0; idx_y<=h2_exp_reco->GetNbinsY()+1; idx_y++){
    std::cout << idx_y << "/" << h2_exp_reco->GetNbinsY() << "\r" << std::flush;
    TH1D* h1_proj = h2_exp_reco->ProjectionX("h1_proj", idx_y, idx_y);
    double h1_proj_integral = h1_proj->Integral();
    
    if(h1_proj_integral > 0.){
      h1_proj->Scale(1./(h1_proj_integral*h1_proj->GetBinWidth(1)));
      
      double mean = h1_proj->GetBinCenter(h1_proj->GetMaximumBin());
      double stddev = h1_proj->GetStdDev();

      f_lognormal->SetParameters(log(mean), 0.25);
      f_lognormal->SetParLimits(0, log(mean-5*stddev), log(mean+5*stddev));

      TFitResultPtr fit_res = nullptr;
      if (h1_proj->GetEntries() > 1500) fit_res = h1_proj->Fit(f_lognormal, "Q", "", std::max(0.,mean-5*stddev), mean+5*stddev);
      h1_proj->SetName(Form("h1_proj_%i_pe", (int)h2_exp_reco->GetYaxis()->GetBinCenter(idx_y)));
      // h1_proj->Write();

      for(int idx_x=0; idx_x<=h2_exp_reco->GetNbinsX()+1; idx_x++){
        h2_exp_reco->SetBinContent(idx_x, idx_y, h1_proj->GetBinContent(idx_x));
      } 

      
      double reco = h2_exp_reco->GetXaxis()->GetBinCenter(idx_y);
      // if (reco > 400) break; // Stop if reco is greater than 40
      if (reco > 3 && h1_proj->GetEntries() > 1500 && fit_res==0){
        reco_pes.push_back(reco);                        err_reco_pes.push_back(0);
        logms.push_back(f_lognormal->GetParameter(0));   err_logms.push_back(f_lognormal->GetParError(0));
        sigmas.push_back(f_lognormal->GetParameter(1));  err_sigmas.push_back(f_lognormal->GetParError(1));
      }
    }
  } // handle h2_exp_reco and fit lognormal_pdf
  // --- TGRAPHS ------------------------------------------------
  TGraphErrors* g_logms = new TGraphErrors(logms.size(), &reco_pes[0], &logms[0], &err_reco_pes[0], &err_logms[0]);
  g_logms->SetTitle("log(MPV) vs reco_pe;#Reco Pe;log(MPV)");
  g_logms->SetName("g_logms");
  TGraphErrors* g_sigmas = new TGraphErrors(sigmas.size(), &reco_pes[0], &sigmas[0], &err_reco_pes[0], &err_sigmas[0]);
  g_sigmas->SetTitle("#sigma vs reco_pe;#Reco Pe;#sigma");
  g_sigmas->SetName("g_sigmas");

  // Fit the parameter trends
  TF1* f_logms_trend = new TF1("f_logms_trend", "[0]+log(x)", trend_thr, 2000.);
  f_logms_trend->SetParNames("const");
  f_logms_trend->SetParameters(.5);
  f_logms_trend->SetNpx(2000);

  g_logms->Fit(f_logms_trend, "", "", trend_thr, 80.);

  TF1* f_sigmas_trend = new TF1("f_sigmas_trend", "[3]+[2]*exp(-[1]*(x-[0]))", trend_thr, 2000.);
  f_sigmas_trend->SetParNames("x0", "lambda", "A", "const");
  f_sigmas_trend->SetParameters(5., 0.1, 1., 0.2);
  f_sigmas_trend->SetNpx(2000);

  g_sigmas->Fit(f_sigmas_trend, "", "", trend_thr, 80.);











  // Get X and Y axis maximum values
  double reco_max = h2_exp_reco->GetXaxis()->GetXmax();
  double exp_max = h2_exp_reco->GetYaxis()->GetXmax();

  TEfficiency* he_hit_prob = nullptr;
  exp_reco_file->GetObject("he_hit_prob", he_hit_prob);
  TF1* f_hit_prob_fit = new TF1("f_hit_prob_fit", sigmoid_sigmoid_erf, 0., 40., 6); // To fit only in a short range
  f_hit_prob_fit->SetParameters(7, 1.7, -5.8, 7., 2.8, 1.2);
  he_hit_prob->Fit(f_hit_prob_fit, "R");
  TF1* f_hit_prob = new TF1("f_hit_prob", sigmoid_sigmoid_erf, 0., 5000., 6); // Extend the function range
  f_hit_prob->SetParameters(f_hit_prob_fit->GetParameters());


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

  // --- EXTRA VARIABLES ------------------------------------------------------
  std::vector<double> x_mis, y_mis, z_mis;
  std::vector<double> dx_mis, dy_mis, dz_mis, de_mis;
 
  // --- LOOP OVER EVENTS -----------------------------------------------------
  double min_term = +1.e+60;  double min_term_fake = +1.e+60;
  double max_term = -1.e+60;  double max_term_fake = -1.e+60;
  size_t evt_counter = 0;     size_t fake_counter = 0;        size_t mismatch_counter = 0;
  size_t combo_counter = 0;
  size_t zeros_term = 0;      size_t zeros_term_fake = 0;
  size_t evts_with_zeros = 0;
  int idx_event = 0;
  Long64_t entry = 0;
  Long64_t nentries = MuRecoTree->GetEntries();
  Double_t old_mu[n_opdet]; Double_t old_true_energy = 0.; Double_t old_vertex_coor[3] = {0.};

  Event true_event, t_event;
  Flash true_flash, t_flash, null_flash;
  std::vector<Event> fake_events;
  std::vector<Flash> fake_flashes;

  while(MuRecoReader.Next()){
    // if (entry == 30000) break;
    entry++;
    if (fake_events.size()<30) {
      t_event = Event(*true_energy, vertex_coor, mu);
      t_flash = Flash(reco_pe, hit_t, hit);
      if (t_flash == null_flash) continue;
      fake_events.push_back(t_event);
      fake_flashes.push_back(t_flash);
      continue;
    }

    t_event = Event(*true_energy, vertex_coor, mu);
    t_flash = Flash(reco_pe, hit_t, hit);
    if (t_flash == null_flash) continue;

    true_flash = fake_flashes[0];
    true_event = fake_events[0];
    std::rotate(fake_flashes.begin(), fake_flashes.begin()+1, fake_flashes.end());
    std::rotate(fake_events.begin(),  fake_events.begin()+1,  fake_events.end());
    fake_flashes[fake_flashes.size()-1] = t_flash;
    fake_events[fake_events.size()-1]   = t_event;

    for (size_t idx_fake=0; idx_fake<fake_events.size(); idx_fake++){
      Event fake_event = fake_events[idx_fake];
      Flash fake_flash = fake_flashes[idx_fake];
      std::vector<double> true_terms, fake_terms;

      if (use_fit_lognormal) {
        true_terms = compute_loglikelihood_terms(true_event, true_flash, f_hit_prob, f_lognormal, f_logms_trend, f_sigmas_trend, g_logms, g_sigmas, trend_thr);
        fake_terms = compute_loglikelihood_terms(fake_event, true_flash, f_hit_prob, f_lognormal, f_logms_trend, f_sigmas_trend, g_logms, g_sigmas, trend_thr);
      }
      else if (use_fit_efficiency) {
        true_terms = compute_loglikelihood_terms(true_event, true_flash, f_hit_prob, h2_exp_reco);
        fake_terms = compute_loglikelihood_terms(fake_event, true_flash, f_hit_prob, h2_exp_reco);
      } else {
        true_terms = compute_loglikelihood_terms(true_event, true_flash, he_hit_prob, h2_exp_reco);
      fake_terms = compute_loglikelihood_terms(fake_event, true_flash, he_hit_prob, h2_exp_reco);
      }

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


      h_Fq->Fill(Fq);
      h_Fq_fake->Fill(Fq_fake);

      if(entry % 1000 == 0)
        std::cout << "Analyzed " << double(entry)/double(nentries)*100. << "% of the events. Fake-matches: " << double(fake_counter)/double(evt_counter)*100. 
          << " % Mismatches: " << double(mismatch_counter)/double(evt_counter)*100.
          << " % Combos: " << double(combo_counter)/double(evt_counter)*100.
          << "%\r" << std::flush;

      bool match;
      if (use_fit_efficiency) match = flash_matcher(true_event, true_flash, fake_event, fake_flash, f_hit_prob, h2_exp_reco);
      else match = flash_matcher(true_event, true_flash, fake_event, fake_flash, he_hit_prob, h2_exp_reco);

      if(Fq_fake < Fq) {
        x_mis.push_back(true_event.vertex_coor[0]);
        y_mis.push_back(true_event.vertex_coor[1]);
        z_mis.push_back(true_event.vertex_coor[2]);
        dx_mis.push_back(true_event.vertex_coor[0] - fake_event.vertex_coor[0]);
        dy_mis.push_back(true_event.vertex_coor[1] - fake_event.vertex_coor[1]);
        dz_mis.push_back(true_event.vertex_coor[2] - fake_event.vertex_coor[2]);
        de_mis.push_back(true_event.true_energy - fake_event.true_energy);
        fake_counter++;
      }
      if(!match) mismatch_counter++;
      if(Fq_fake < Fq && !match) combo_counter++;
      evt_counter++;
    }
  } // end loop over events
  
  // --- OUTPUT ---------------------------------------------------------------
  TH1D* h_x_mis = new TH1D("h_x_mis",Form("%s;%s;%s","h_x_mis","x [cm]","counts"),
                          100, *std::min_element(x_mis.begin(), x_mis.end()), *std::max_element(x_mis.begin(), x_mis.end()));

  TH1D* h_y_mis = new TH1D("h_y_mis",Form("%s;%s;%s","h_y_mis","y [cm]","counts"),
                          100, *std::min_element(y_mis.begin(), y_mis.end()), *std::max_element(y_mis.begin(), y_mis.end()));

  TH1D* h_z_mis = new TH1D("h_z_mis",Form("%s;%s;%s","h_z_mis","z [cm]","counts"),
                          100, *std::min_element(z_mis.begin(), z_mis.end()), *std::max_element(z_mis.begin(), z_mis.end()));

  TH1D* h_dx_mis = new TH1D("h_dx_mis",Form("%s;%s;%s","h_dx_mis","#Deltax [cm]","counts"),
                          100, *std::min_element(dx_mis.begin(), dx_mis.end()), *std::max_element(dx_mis.begin(), dx_mis.end()));

  TH1D* h_dy_mis = new TH1D("h_dy_mis",Form("%s;%s;%s","h_dy_mis","#Deltay [cm]","counts"),
                          100, *std::min_element(dy_mis.begin(), dy_mis.end()), *std::max_element(dy_mis.begin(), dy_mis.end()));

  TH1D* h_dz_mis = new TH1D("h_dz_mis",Form("%s;%s;%s","h_dz_mis","#Deltaz [cm]","counts"),
                          100, *std::min_element(dz_mis.begin(), dz_mis.end()), *std::max_element(dz_mis.begin(), dz_mis.end()));

  TH1D* h_de_mis = new TH1D("h_de_mis",Form("%s;%s;%s","h_de_mis","#DeltaE [MeV]","counts"),
                          100, *std::min_element(de_mis.begin(), de_mis.end()), *std::max_element(de_mis.begin(), de_mis.end()));

  TH2D* h2_de_dx_mis = new TH2D("h2_de_dx_mis",Form("%s;%s;%s","h2_de_dx_mis","#Deltax [cm]","#DeltaE [MeV]"),
                          100, *std::min_element(dx_mis.begin(), dx_mis.end()), *std::max_element(dx_mis.begin(), dx_mis.end()),
                          100, *std::min_element(de_mis.begin(), de_mis.end()), *std::max_element(de_mis.begin(), de_mis.end()));

  TH2D* h2_de_dy_mis = new TH2D("h2_de_dy_mis",Form("%s;%s;%s","h2_de_dy_mis","#Deltay [cm]","#DeltaE [MeV]"),
                          100, *std::min_element(dy_mis.begin(), dy_mis.end()), *std::max_element(dy_mis.begin(), dy_mis.end()),
                          100, *std::min_element(de_mis.begin(), de_mis.end()), *std::max_element(de_mis.begin(), de_mis.end()));

  TH2D* h2_de_dz_mis = new TH2D("h2_de_dz_mis",Form("%s;%s;%s","h2_de_dz_mis","#Deltaz [cm]","#DeltaE [MeV]"),
                          100, *std::min_element(dz_mis.begin(), dz_mis.end()), *std::max_element(dz_mis.begin(), dz_mis.end()),
                          100, *std::min_element(de_mis.begin(), de_mis.end()), *std::max_element(de_mis.begin(), de_mis.end()));

  for (size_t i = 0; i < x_mis.size(); ++i) {
    h_x_mis->Fill(x_mis[i]);
    h_y_mis->Fill(y_mis[i]);
    h_z_mis->Fill(z_mis[i]);
    h_dx_mis->Fill(dx_mis[i]);
    h_dy_mis->Fill(dy_mis[i]);
    h_dz_mis->Fill(dz_mis[i]);
    h_de_mis->Fill(de_mis[i]);
    h2_de_dx_mis->Fill(dx_mis[i], de_mis[i]);
    h2_de_dy_mis->Fill(dy_mis[i], de_mis[i]);
    h2_de_dz_mis->Fill(dz_mis[i], de_mis[i]);
  }

  std::cout << "\nFake/Evts:   " <<fake_counter << "/" <<evt_counter
            << "\nFake/Evts %: " <<fake_counter/double(evt_counter)*100. 
            << "\nMismatch/Evts:   " <<mismatch_counter << "/" <<evt_counter
            << "\nMismatch/Evts %: " <<mismatch_counter/double(evt_counter)*100.
            << "\nCombo/Evts:   " <<combo_counter << "/" <<evt_counter
            << "\nCombo/Evts %: " <<combo_counter/double(evt_counter)*100.
            << "\ndE dx correlation:\t" << h2_de_dx_mis->GetCorrelationFactor()
            << "\ndE dy correlation:\t" << h2_de_dy_mis->GetCorrelationFactor()
            << "\ndE dz correlation:\t" << h2_de_dz_mis->GetCorrelationFactor()
            << std::endl;  

  g_logms->Write();
  g_sigmas->Write();
  h_x_mis->Write();
  h_y_mis->Write();
  h_z_mis->Write();
  h_dx_mis->Write();
  h_dy_mis->Write();
  h_dz_mis->Write();
  h_de_mis->Write();
  h2_de_dx_mis->Write();
  h2_de_dy_mis->Write();
  h2_de_dz_mis->Write();
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
