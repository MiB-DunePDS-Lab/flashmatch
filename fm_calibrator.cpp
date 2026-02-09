#include "RtypesCore.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeReader.h"
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <numeric>
#include <tuple>
#include <vector>

#include "Utils.hpp"


// --- HARD CODE HERE ----------------
int max_nfiles = 8000; // Maximum number of files to analyze
std::string ana_folder_name =  "./ophitfinder_edep_debug/ana/";
std::string base_ana_file_name = "solar_ana_dune10kt_1x2x6_hist_";

// -----------------------------------
void fm_calibrator(){
  // --- EXTRA VARIABLES -------------------------------------------------------
  std::vector<std::tuple<float, float, float, float>> true_calib_info; // 
  float drift_velocity = 0.0; // Drift velocity in cm/tick, will be set later from the config_tree

  // --- LOOP OVER ANA FILES ---------------------------------------------------
  std::filesystem::path ana_folder(ana_folder_name);
  int nfile_in_folder = std::distance(std::filesystem::directory_iterator(ana_folder), std::filesystem::directory_iterator{});
  int nfile_to_analyze = std::min(nfile_in_folder, max_nfiles);
  int nfile_analyzed = 0;
  int idx_file = 0;
  while(nfile_analyzed < nfile_to_analyze){
    // --- ANA STUFF -----------------------------------------------------------
    std::string ana_file_name = std::string(ana_folder_name+base_ana_file_name+idx_file+".root"); 
    // Check whether both files exist
    idx_file++;
    if(!std::filesystem::exists(ana_file_name)) continue;
    nfile_analyzed++;
    std::cout <<nfile_analyzed<<"--"<< ana_file_name << std::endl;

    TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(ana_file->Get("solarnuana/MCTruthTree"));
    TTreeReader treeReader(tree);

    TTreeReaderValue<float> E_true(treeReader, "SignalParticleE");
    TTreeReaderValue<float> x_true(treeReader, "SignalParticleX");
    TTreeReaderValue<std::vector<int>> EDepElectron(treeReader, "TSignalElectronDepList");
   
    if (nfile_analyzed==1){
      TTree* config_tree = static_cast<TTree*>(ana_file->Get("solarnuana/ConfigTree"));
      int DetectorSizeX = 6000; // Hardcoded for now, should be read from config_tree
      int DetectorDriftTime = 2; // Hardcoded for now, should be read from config_tree
      config_tree->SetBranchAddress("DetectorSizeX", &DetectorSizeX);
      config_tree->SetBranchAddress("DetectorDriftTime", &DetectorDriftTime);
      config_tree->GetEntry(0);
      drift_velocity = float(DetectorSizeX) / float(DetectorDriftTime);
    }
    
    while (treeReader.Next()) {
      float total_charge = float(std::accumulate((*EDepElectron).begin(), (*EDepElectron).end(), 0));
      float driftTime = (*x_true)/drift_velocity;
      true_calib_info.push_back(std::make_tuple(total_charge, *E_true, driftTime, total_charge/(*E_true)));
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
                                      200, min_drift, max_drift, 200, min_QperE, max_QperE);
  for (const auto& [charge, etrue, drifttime, QperE] : true_calib_info){
    h2_QperE_driftTime->Fill(drifttime, QperE);
  }
  TGraphErrors* g_QperE_driftTime = th2d_to_tgraph_mpv(h2_QperE_driftTime, "g_QperE_driftTime");
  TF1* f_Qcorr = new TF1("f_Qcorr", "[0]*exp(-x/[1])");
  f_Qcorr->SetParNames("A", "#tau");
  f_Qcorr->SetParameters(TMath::MaxElement(g_QperE_driftTime->GetN(), g_QperE_driftTime->GetY()), 20000);
  float dDrift = max_drift - min_drift;
  g_QperE_driftTime->Fit("f_Qcorr", "", "", min_drift+0.05*dDrift, max_drift-0.05*dDrift);


  TH2D* h2_Qcorr_Etrue = new TH2D("h2_Qcorr_Etrue",
                                      Form("%s;%s;%s", "h2_Qcorr_Etrue", "E_{true} [MeV]", "Qcorr"),
                                      200, min_etrue, max_etrue, 200, min_charge, max_charge);

  float corr_lambda = 1/f_Qcorr->GetParameter(1);
  for (const auto& [charge, etrue, drifttime, QperE] : true_calib_info){
    h2_Qcorr_Etrue->Fill(etrue, charge*exp(drifttime*corr_lambda));
  }
  TGraphErrors* g_Qcorr_Etrue = th2d_to_tgraph_mpv(h2_Qcorr_Etrue, "g_Qcorr_Etrue");
  TF1* f_Calib = new TF1("f_Calib", "pol1");
  f_Calib->SetParNames("c", "slope");
  g_Qcorr_Etrue->Fit(f_Calib, "", "", min_etrue, max_etrue);

  TTree* calib_tree = new TTree("calib_tree", "calib_tree");
  Float_t calib_c = 0.; Float_t calib_slope = 0.;
  calib_tree->Branch("calib_c", &calib_c, "calib_c/F"); 
  calib_tree->Branch("calib_slope", &calib_slope, "calib_slope/F");
  calib_tree->Branch("drift_velocity", &drift_velocity, "drift_velocity/F");
  calib_tree->Branch("corr_lambda", &corr_lambda, "corr_lambda/F");
  calib_c = f_Calib->GetParameter(0); calib_slope = f_Calib->GetParameter(1);
  calib_tree->Fill();

  TFile* calibrator_file = TFile::Open("calibrator.root", "RECREATE");
  calibrator_file->cd();
  calib_tree->Write();
  h2_QperE_driftTime->Write();
  g_QperE_driftTime->Write();
  h2_Qcorr_Etrue->Write();
  g_Qcorr_Etrue->Write();
  calibrator_file->Close();

  return;
}
