#include "RtypesCore.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"


#include "Utils.hpp"
#include "flash_matcher.hpp"
#include <cstdio>
#include <vector>

void fm_Offline(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLcconfigs f = load_ana_config("./config.json");
  int max_nfiles               = f.max_nfiles;
  std::string input_dir        = f.input_dir;
  TString visibility_file_name = f.visibility_file_name;
  double trend_thr             = f.trend_thr;
  float light_yield            = f.light_yield;
  float arapuca_pde            = f.arapuca_pde;
  const size_t n_opdet         = f.n_opdet;
  float LY_times_PDE           = light_yield * arapuca_pde;
  
  // TFile* ana_file = TFile::Open("../MLL_FM/debugs/files/solar_ana_dune10kt_1x2x6_hist.root", "READ");
  TFile* ana_file = TFile::Open((input_dir+"files/SolarNuAnaMerged.root").c_str(), "READ");
  TTree* tree = static_cast<TTree*>(ana_file->Get("solarnuana/SolarNuAnaTree"));
  std::vector<size_t> MaxChargeIndxs = take_max_charge_indices(tree, "Event", "Charge");
  TTreeReader treeReader(tree);
  // True info
  TTreeReaderValue<float> SignalParticleX(treeReader, "SignalParticleX");
  // Cluster stuff
  TTreeReaderValue<float> RecoY(treeReader, "RecoY");
  TTreeReaderValue<float> RecoZ(treeReader, "RecoZ");
  TTreeReaderValue<float> Charge(treeReader,"Charge");
  TTreeReaderValue<float> Time(treeReader,  "Time");
  // Adj. flashes stuff
  TTreeReaderArray<float> MAdjFlashResidual(treeReader, "AdjOpFlashResidual");
  TTreeReaderArray<float> MAdjFlashTime(treeReader, "AdjOpFlashTime");
  TTreeReaderArray<float> MAdjFlashPur(treeReader, "AdjOpFlashPur");
  TTreeReaderArray<int>   MAdjFlashNHits(treeReader, "AdjOpFlashNHits");
  TTreeReaderArray<float> MAdjFlashPE(treeReader, "AdjOpFlashPE");
  TTreeReaderArray<float> MAdjFlashPEperOpDet(treeReader, "AdjOpFlashPEperOpDet");
  // Matched flash stuff
  TTreeReaderValue<float> MatchedOpFlashTime(treeReader, "MatchedOpFlashTime");
  TTreeReaderValue<float> MatchedOpFlashPE(treeReader, "MatchedOpFlashPE");
  TTreeReaderArray<float> MatchedOpFlashPEperOpDet(treeReader, "MatchedOpFlashPEperOpDet");
  TTreeReaderValue<bool> MatchedOpFlashCorrectly(treeReader, "MatchedOpFlashCorrectly");
  
  // --- LikelihoodComputer -----------------------------------------------------
  TFile* calib_file = TFile::Open((input_dir+"MLL_Calibrator.root").c_str(), "READ");
  TTree* calib_tree = static_cast<TTree*>(calib_file->Get("calib_tree"));
  Float_t calib_c = 0.;         Float_t calib_slope = 0.;
  Float_t drift_velocity = 0.0; Float_t corr_lambda = 0.0;
  calib_tree->SetBranchAddress("calib_c", &calib_c); 
  calib_tree->SetBranchAddress("calib_slope", &calib_slope);
  calib_tree->SetBranchAddress("drift_velocity", &drift_velocity);
  calib_tree->SetBranchAddress("corr_lambda", &corr_lambda);
  calib_tree->GetEntry(0);
  
  TFile* parametrizer_file  = TFile::Open((input_dir+"MLL_Parametrizer.root").c_str(), "READ");
  TF1* f_reco_prob          = static_cast<TF1*>(parametrizer_file->Get("f_reco_prob"));
  TF1* f_RecoExpDistr       = static_cast<TF1*>(parametrizer_file->Get("f_RecoExpDistr"));
  TF1* f_par1_trend         = static_cast<TF1*>(parametrizer_file->Get("f_par1_trend"));
  TF1* f_par2_trend         = static_cast<TF1*>(parametrizer_file->Get("f_par2_trend"));
  TH2D* h2_exp_reco         = static_cast<TH2D*>(parametrizer_file->Get("h2_exp_reco"));
  TGraphErrors* g_par1      = static_cast<TGraphErrors*>(parametrizer_file->Get("g_par1")); 
  TGraphErrors* g_par2      = static_cast<TGraphErrors*>(parametrizer_file->Get("g_par2"));
 
  LikelihoodComputer likelihood_computer(
    visibility_file_name, // Visibility file name
    n_opdet,              // Number of optical detectors
    drift_velocity,       // Drift velocity
    LY_times_PDE,         // Light yield times photo detector efficiency
    f_reco_prob,          // Reconstruction probability function
    f_RecoExpDistr,          // Lognormal function for extrapolation
    f_par1_trend,        // Trend function for logm
    f_par2_trend,       // Trend function for sigma
    g_par1,              // Graph for logm values
    g_par2,             // Graph for sigma values
    trend_thr,            // Threshold for trend
    calib_c,              // Calibration constant
    calib_slope,          // Calibration slope
    corr_lambda,          // Correction lambda value
    h2_exp_reco
  );

  // --- OUTPUT PLOTS ---------------------------------------------------------
  TFile* output_file = TFile::Open((input_dir+"/MLL_Offline.root").c_str(), "RECREATE");
  Long64_t nn = tree->Draw("SignalParticleX", "", "goff");
  float max_drift = TMath::MaxElement(tree->GetSelectedRows(), tree->GetV1());
  TEfficiency* he_EffvsDrift       = new TEfficiency("he_EffvsDrift",       "Efficiency vs Drift; Drift [cm]; Efficiency", 20, 0., double(max_drift));
  TEfficiency* he_EffvsDrift_cheat = new TEfficiency("he_EffvsDrift_cheat", "Efficiency vs Drift; Drift [cm]; Efficiency", 20, 0., double(max_drift));
  TEfficiency* he_EffvsDrift_maxpe = new TEfficiency("he_EffvsDrift_maxpe", "Efficiency vs Drift; Drift [cm]; Efficiency", 20, 0., double(max_drift));
  
  TH1D* h_LL = new TH1D("h_LL",Form("%s;%s;%s","h_LL","LL","counts"),
                        200, 0., 20.);
  TH1D* h_LL_target = new TH1D("h_LL_target",Form("%s;%s;%s","h_LL_target","LL","counts"),
                        200, 0., 20.);
  TH1D* h_LL_scale = new TH1D("h_LL_scale",Form("%s;%s;%s","h_LL_scale","LL","counts"),
                        200, 0., 20.);
  TH1D* h_LL_bkg = new TH1D("h_LL_bkg",Form("%s;%s;%s","h_LL_bkg","LL_bkg","counts"),
                             200, 0., 20.);
  TH1D* h_TargetRecoTerms = new TH1D("h_TargetRecoTerms",Form("%s;%s;%s","h_TargetRecoTerms","Reco Terms for Target Match","counts"),
                                    200, 0, 100);

  TH1D* h_TargetNoRecoTerms = new TH1D("h_TargetNoRecoTerms",Form("%s;%s;%s","h_TargetNoRecoTerms","NoReco Terms for Target Match","counts"),
                                    200, 0, 100);

  TH1D* h_SignalRecoTerms = new TH1D("h_SignalRecoTerms",Form("%s;%s;%s","h_SignalRecoTerms","Reco Terms for Signal Match","counts"),
                                    200, 0, 100);

  TH1D* h_SignalNoRecoTerms = new TH1D("h_SignalNoRecoTerms",Form("%s;%s;%s","h_SignalNoRecoTerms","NoReco Terms for Signal Match","counts"),
                                    200, 0, 100);

  TH1D* h_bkgRecoTerms = new TH1D("h_bkgRecoTerms",Form("%s;%s;%s","h_bkgRecoTerms","Reco Terms for bkg Match","counts"),
                                    200, 0, 100);

  TH1D* h_bkgNoRecoTerms = new TH1D("h_bkgNoRecoTerms",Form("%s;%s;%s","h_bkgNoRecoTerms","NoReco Terms for bkg Match","counts"),
                                    200, 0, 100);

 
  // --- LOOP OVER ENTRIES ----------------------------------------------------
  std::vector<float> reco_pes(n_opdet, 0.);
  std::vector<float> reco_terms;
  std::vector<float> noreco_terms;
  size_t n_try = 0;
  size_t n_match = 0;
  size_t n_match_cheating = 0;
  size_t n_noflash = 0;
  for (auto& idx_entry : MaxChargeIndxs){
    treeReader.SetEntry(idx_entry);
  // while (treeReader.Next()) {
    if (MAdjFlashPE.GetSize()==0) {
      // std::cout << "No flash PE information available!" << std::endl;
      he_EffvsDrift->Fill(false, abs(*SignalParticleX));
      he_EffvsDrift_cheat->Fill(false, abs(*SignalParticleX));
      he_EffvsDrift_maxpe->Fill(false, abs(*SignalParticleX));
      n_try++; n_noflash++;
      continue;
    }

    int idx_selected_max = -1;
    int idx_selected_cheat = -1;
    int idx_selected = -1;
    float max_pe = -10000.;
    float value_selected = std::numeric_limits<float>::max();
    float purity = -1.;
    ClusterTPC tpc_cluster = ClusterTPC(*Charge, *Time, *RecoY, *RecoZ);

    for (size_t idx_flash=0; idx_flash<MAdjFlashResidual.GetSize(); idx_flash++){
      for(size_t j=0; j<n_opdet; j++) reco_pes[j] = MAdjFlashPEperOpDet.At(idx_flash*n_opdet + j);

      float flash_time = MAdjFlashTime.At(idx_flash);
      float flash_residual = MAdjFlashResidual.At(idx_flash);
      float flash_pe = MAdjFlashPE.At(idx_flash);
      float flash_nhits = MAdjFlashNHits.At(idx_flash);
      ClusterPDS pds_cluster = ClusterPDS(flash_time, &reco_pes);
    
      float offline_likelihood = -likelihood_computer.GetLikelihoodMatch(tpc_cluster, pds_cluster, reco_terms, noreco_terms);

      if (MAdjFlashPE.At(idx_flash) == *MatchedOpFlashPE && MAdjFlashTime.At(idx_flash) == *MatchedOpFlashTime){
        h_LL_target->Fill(offline_likelihood);
        for (auto& term : reco_terms) h_TargetRecoTerms->Fill(-term);
        for (auto& term : noreco_terms) h_TargetNoRecoTerms->Fill(-term);
      }

      if (MAdjFlashPur.At(idx_flash)){
        h_LL->Fill(offline_likelihood);
        for (auto& term : reco_terms) h_SignalRecoTerms->Fill(-term);
        for (auto& term : noreco_terms) h_SignalNoRecoTerms->Fill(-term);
      }
      else{
        h_LL_bkg->Fill(offline_likelihood);
        for (auto& term : reco_terms) h_bkgRecoTerms->Fill(-term);
        for (auto& term : noreco_terms) h_bkgNoRecoTerms->Fill(-term);
      }

      // Selection ------------------------------------------------------------
      // if (
      //   flash_pe < 5. ||
      //   flash_nhits < 4
      // ) continue;

      // Sergio's
      if (flash_pe > max_pe){
        max_pe = flash_pe;
        idx_selected_max = idx_flash;
      }
      // Cheating: select the flash that matches the true flash
      if (*MatchedOpFlashCorrectly) idx_selected_cheat = idx_flash;

      // Selection criteria for the matched flash
      float value = offline_likelihood;
      // float value = offline_likelihood/(flash_nhits*flash_nhits); // 60%
      if (value_selected == std::numeric_limits<float>::max() || value < value_selected){
        value_selected = value;
        idx_selected = idx_flash;
        purity = MAdjFlashPur.At(idx_flash);
      }
    } // End loop over adjacent flashes
    
    n_try++;
    
    he_EffvsDrift->Fill((purity>0.), abs(*SignalParticleX));
    he_EffvsDrift_cheat->Fill(*MatchedOpFlashCorrectly, abs(*SignalParticleX));
    if (idx_selected_max>=0) he_EffvsDrift_maxpe->Fill((MAdjFlashPur.At(idx_selected_max))>0., abs(*SignalParticleX));
    else he_EffvsDrift_maxpe->Fill(false, abs(*SignalParticleX));

    if (purity > 0.) n_match++;
    if(*MatchedOpFlashCorrectly) n_match_cheating++;
  } // End loop over entries
  std::cout << "Method  : " << float(n_match)/float(n_try) << std::endl;
  std::cout << "Cheating: " << float(n_match_cheating)/float(n_try) << std::endl;
  std::cout << "N tries : " << n_try << std::endl;
  std::printf("No flash: %zu (%.3f%%)\n", n_noflash, float(n_noflash)/float(n_try)*100);
  std::printf("Target matches: %.0f (%.3f%%) of reco_terms and %.0f (%.3f%%) of non-reco\n",
              h_TargetRecoTerms->GetEntries(), h_TargetRecoTerms->GetEntries()/(n_opdet*h_LL_target->GetEntries())*100,
              h_TargetNoRecoTerms->GetEntries(), h_TargetNoRecoTerms->GetEntries()/(n_opdet*h_LL_target->GetEntries())*100);
  std::printf("True matches : %.0f (%.3f%%) of reco_terms and %.0f (%.3f%%) of non-reco\n",
              h_SignalRecoTerms->GetEntries(), h_SignalRecoTerms->GetEntries()/(n_opdet*h_LL->GetEntries())*100,
              h_SignalNoRecoTerms->GetEntries(), h_SignalNoRecoTerms->GetEntries()/(n_opdet*h_LL->GetEntries())*100);
  std::printf("False matches %.0f (%.3f%%) of reco_terms and %.0f (%.3f%%) of non-reco\n",
              h_bkgRecoTerms->GetEntries(), h_bkgRecoTerms->GetEntries()/(n_opdet*h_LL_bkg->GetEntries())*100,
              h_bkgNoRecoTerms->GetEntries(), h_bkgNoRecoTerms->GetEntries()/(n_opdet*h_LL_bkg->GetEntries())*100);

  output_file->WriteTObject(he_EffvsDrift);
  output_file->WriteTObject(he_EffvsDrift_cheat);
  output_file->WriteTObject(he_EffvsDrift_maxpe);
  h_LL->Write();
  h_LL_target->Write();
  h_LL_scale->Write();
  h_LL_bkg->Write();
  h_TargetRecoTerms->Write();
  h_TargetNoRecoTerms->Write();
  h_SignalRecoTerms->Write();
  h_SignalNoRecoTerms->Write();
  h_bkgRecoTerms->Write();
  h_bkgNoRecoTerms->Write();
  output_file->Close();

  return;
}
