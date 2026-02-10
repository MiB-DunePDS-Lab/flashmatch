
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

#include "TTreeReader.h"
#include "Utils.hpp"
void fm_debug(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLcconfigs f = load_ana_config("./config.json");
  int max_nfiles                 = f.max_nfiles;
  std::string input_dir          = f.input_dir;
  double fit_Qcorr_Etrue_low     = f.fit_Qcorr_Etrue_low;
  double fit_Qcorr_Etrue_up      = f.fit_Qcorr_Etrue_up;

  // --- EXTRA VARIABLES -------------------------------------------------------
  std::vector<std::tuple<float, float, float, float>> true_calib_info; // 
  float drift_velocity = 360./2244.44; // HARD CODED: Drift velocity in cm/tick

  // --- LOOP OVER ANA FILES ---------------------------------------------------
  std::vector<std::string> ana_files = get_list_of_files_in_folder(input_dir, ".root");
  int nfile_to_analyze = std::min(int(ana_files.size()), max_nfiles);
  int nfile_analyzed = 0; int idx_file = 0;
  while(nfile_analyzed < nfile_to_analyze){
    // --- ANA STUFF -----------------------------------------------------------
    std::string ana_file_name = input_dir+ana_files[idx_file];
    idx_file++;
    if(!std::filesystem::exists(ana_file_name)) continue;
    nfile_analyzed++;
    /* if (idx_file % 10 == 0) */ std::cout <<nfile_analyzed<<"--"<< ana_file_name << "\r" << std::flush;

    TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(ana_file->Get("solarnuana/SolarNuAnaTree"));
    TTreeReader treeReader(tree);

    TTreeReaderValue<float> E_true(treeReader, "SignalParticleE");
    TTreeReaderValue<float> x_true(treeReader, "SignalParticleX");
    TTreeReaderValue<float> Charge(treeReader, "Charge");
    TTreeReaderValue<bool> MatchedOpFlashCorrectly(treeReader, "MatchedOpFlashCorrectly");
   
    while (treeReader.Next()) {
      if (!(*MatchedOpFlashCorrectly)) continue;
      float driftTime = abs((*x_true)/drift_velocity);
      true_calib_info.push_back(std::make_tuple(*Charge, *E_true, driftTime, *Charge/(*E_true)));
    }
    ana_file->Close();
  }

  float min_charge = get_minimum_from_tuples(true_calib_info, 0); float max_charge = get_maximum_from_tuples(true_calib_info, 0);
  float min_etrue  = get_minimum_from_tuples(true_calib_info, 1); float max_etrue  = get_maximum_from_tuples(true_calib_info, 1);
  float min_drift  = get_minimum_from_tuples(true_calib_info, 2); float max_drift  = get_maximum_from_tuples(true_calib_info, 2);
  float min_QperE  = get_minimum_from_tuples(true_calib_info, 3); float max_QperE  = get_maximum_from_tuples(true_calib_info, 3);

  // Create charge-per-energy vs drift-time plots
  TH2D* h2_QperE_driftTime = new TH2D("h2_QperE_driftTime",
                                      Form("%s;%s;%s", "h2_QperE_driftTime", "Drift Time", "QperE"),
                                      100, min_drift, max_drift, 50, min_QperE, max_QperE);
  for (const auto& [charge, etrue, drifttime, QperE] : true_calib_info){
    h2_QperE_driftTime->Fill(drifttime, QperE);
  }

  TGraphErrors* g_QperE_driftTime = th2d_to_tgraph_mpv(h2_QperE_driftTime, "g_QperE_driftTime");
  g_QperE_driftTime->SetTitle("QperE vs Drift Time;Drift Time [ticks];QperE [Charge/MeV]");
  
  TF1* f_Qcorr = new TF1("f_Qcorr", "[0]*exp(-x/[1])");
  f_Qcorr->SetParNames("A", "#tau");
  f_Qcorr->SetParameters(TMath::MaxElement(g_QperE_driftTime->GetN(), g_QperE_driftTime->GetY()), 20000);
  
  float dDrift = max_drift - min_drift;
  g_QperE_driftTime->Fit("f_Qcorr", "", "", min_drift+0.05*dDrift, max_drift-0.05*dDrift);


  TH2D* h2_Qcorr_Etrue = new TH2D("h2_Qcorr_Etrue",
                                  Form("%s;%s;%s", "h2_Qcorr_Etrue", "E_{true} [MeV]", "Qcorr"),
                                  100, min_etrue, max_etrue, 50, min_charge, max_charge);
  TH2D* h2_Qcorr_Etrue_low = new TH2D("h2_Qcorr_Etrue_low",
                                  Form("%s;%s;%s", "h2_Qcorr_Etrue_low", "E_{true} [MeV]", "Qcorr"),
                                  100, min_etrue, max_etrue, 50, min_charge, max_charge);
  TH2D* h2_Qcorr_Etrue_high = new TH2D("h2_Qcorr_Etrue_high",
                                  Form("%s;%s;%s", "h2_Qcorr_Etrue_high", "E_{true} [MeV]", "Qcorr"),
                                  100, min_etrue, max_etrue, 50, min_charge, max_charge);
  TH2D* h2_Qcorr_Etrue_sum = new TH2D("h2_Qcorr_Etrue_sum",
                                  Form("%s;%s;%s", "h2_Qcorr_Etrue_sum", "E_{true} [MeV]", "Qcorr"),
                                  100, min_etrue, max_etrue, 50, min_charge, max_charge);
  TH2D* h2_duplicate_Etrue_xtrue = new TH2D("h2_duplicate_Etrue_xtrue",
                                  Form("%s;%s;%s", "h2_duplicate_Etrue_xtrue", "E_{true} [MeV]", "x_{true} [cm]"),
                                  100, min_etrue, max_etrue, 100, min_drift*drift_velocity, max_drift*drift_velocity);
  float corr_lambda = 1/f_Qcorr->GetParameter(1);
  float old_E = 0., old_charge = 0.;
  for (const auto& [charge, etrue, drifttime, QperE] : true_calib_info){
    h2_Qcorr_Etrue->Fill(etrue, charge*exp(drifttime*corr_lambda));
    if(old_E == etrue){
      h2_duplicate_Etrue_xtrue->Fill(etrue, drifttime*drift_velocity);
      if (old_E < etrue){
        h2_Qcorr_Etrue_low->Fill(old_E,  old_charge*exp(drifttime*corr_lambda));
        h2_Qcorr_Etrue_high->Fill(etrue, charge*exp(drifttime*corr_lambda));
        h2_Qcorr_Etrue_sum->Fill(etrue, (charge+old_charge)*exp(drifttime*corr_lambda));
      }
      else{
        h2_Qcorr_Etrue_low->Fill(etrue, charge*exp(drifttime*corr_lambda));
        h2_Qcorr_Etrue_high->Fill(old_E,  old_charge*exp(drifttime*corr_lambda));
        h2_Qcorr_Etrue_sum->Fill(etrue, (charge+old_charge)*exp(drifttime*corr_lambda));
      }
    }
    old_E = etrue; old_charge = charge;
  }

  TGraphErrors* g_Qcorr_Etrue = th2d_to_tgraph_mpv(h2_Qcorr_Etrue, "g_Qcorr_Etrue");
  g_Qcorr_Etrue->SetTitle("Calib Graph;E_{true} [MeV];Q_{corr}");
  TF1* f_Calib = new TF1("f_Calib", "pol1");
  f_Calib->SetParNames("c", "slope");
  g_Qcorr_Etrue->Fit(f_Calib, "", "", fit_Qcorr_Etrue_low, fit_Qcorr_Etrue_up);

  TTree* calib_tree = new TTree("calib_tree", "calib_tree");
  Float_t calib_c = 0.; Float_t calib_slope = 0.;
  calib_tree->Branch("calib_c",        &calib_c, "calib_c/F"); 
  calib_tree->Branch("calib_slope",    &calib_slope, "calib_slope/F");
  calib_tree->Branch("drift_velocity", &drift_velocity, "drift_velocity/F");
  calib_tree->Branch("corr_lambda",    &corr_lambda, "corr_lambda/F");
  calib_c = f_Calib->GetParameter(0); calib_slope = f_Calib->GetParameter(1);
  calib_tree->Fill();

  TFile* calibrator_file = TFile::Open("MLL_debug.root", "RECREATE");
  calibrator_file->cd();
  calib_tree->Write();
  h2_QperE_driftTime->Write();
  g_QperE_driftTime->Write();
  h2_duplicate_Etrue_xtrue->Write();
  h2_Qcorr_Etrue->Write();
  h2_Qcorr_Etrue_low->Write();
  h2_Qcorr_Etrue_high->Write();
  h2_Qcorr_Etrue_sum->Write();
  g_Qcorr_Etrue->Write();
  calibrator_file->Close();

  return;
}
