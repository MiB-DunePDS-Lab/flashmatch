#include "RtypesCore.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeReader.h"
#include <cstddef>
#include <filesystem>
#include <tuple>
#include <vector>

#include "TTreeReaderValue.h"
#include "Utils.hpp"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void fm_calibrator(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLConfigs mll_conf = load_ana_config("./configs/ana_config.json");
  std::string ana_file_name      = mll_conf.ana_file_name;
  std::string sample_config_file = mll_conf.sample_config_file;
  std::string visibility_dir     = mll_conf.visibility_dir;
  double fit_Qcorr_Etrue_low     = mll_conf.fit_Qcorr_Etrue_low;
  double fit_Qcorr_Etrue_up      = mll_conf.fit_Qcorr_Etrue_up;
  float q_cut_high               = mll_conf.q_cut_high;

  SampleConfigs sample_conf = load_sample_config("./configs/"+sample_config_file);
  std::string input_dir          = sample_conf.input_dir;
  std::string geom_identifier    = sample_conf.geom_identifier;
  float q_cut_low                = sample_conf.q_cut_low;

  TString visibility_file_name     = TString(visibility_dir+"dunevis_"+geom_identifier+".root");
  
  // --- EXTRA VARIABLES -------------------------------------------------------
  std::vector<std::tuple<float, float, float, float>> true_calib_info; // 
  float drift_velocity = 360./2244.44; // HARD CODED: Drift velocity in cm/tick

  TFile* visibility_file = TFile::Open(visibility_file_name, "READ");
  TTree* tDimensions = (TTree*)visibility_file->Get("photovisAr/tDimensions");
  Double_t coor_dim[3] = {0.};
  tDimensions->SetBranchAddress("dimension", coor_dim);
  tDimensions->GetEntry(0);
  float tpc_min[3] = {float(coor_dim[0]), float(coor_dim[1]), float(coor_dim[2])}; // x,y,z
  tDimensions->GetEntry(1);
  float tpc_max[3] = {float(coor_dim[0]), float(coor_dim[1]), float(coor_dim[2])}; // x,y,z
  float anode_x = (geom_identifier == "dune10kt") ? tpc_max[0]+tpc_min[0] : tpc_max[0];
  

  ana_file_name = input_dir+ana_file_name;
  if(!std::filesystem::exists(ana_file_name)) printf("File %s does not exist. Exiting.\n", ana_file_name.c_str());

  TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
  TTree* tree = static_cast<TTree*>(ana_file->Get("solarnuana/SolarNuAnaTree"));
  std::vector<size_t> MaxChargeIndxs = take_max_charge_indices(tree, "Event", "Charge");
  printf("Total entries in tree: %lld\n", tree->GetEntries());
  printf("Entries passing max charge selection: %lu\n", MaxChargeIndxs.size());

  TTreeReader treeReader(tree);
  TTreeReaderValue<float> E_true(treeReader, "SignalParticleE");
  TTreeReaderValue<float> x_true(treeReader, "SignalParticleX");
  TTreeReaderValue<float> Charge(treeReader, "Charge");
  TTreeReaderValue<bool> MatchedOpFlashCorrectly(treeReader, "MatchedOpFlashCorrectly");

  for (auto& idx_entry : MaxChargeIndxs){
    treeReader.SetEntry(idx_entry);
    if (!(*MatchedOpFlashCorrectly) || *Charge<q_cut_low || *Charge>q_cut_high) continue;
    float driftTime = abs((*x_true-anode_x)/drift_velocity);
    true_calib_info.push_back(std::make_tuple(*Charge, *E_true, driftTime, *Charge/(*E_true)));
  }
  printf("Entries passing cuts: %lu\n", true_calib_info.size());
  ana_file->Close();

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
  //std::cout << "min_QperE"<< min_QperE << std::endl;
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
  float corr_lambda = 1/f_Qcorr->GetParameter(1);
  for (const auto& [charge, etrue, drifttime, QperE] : true_calib_info){
    h2_Qcorr_Etrue->Fill(etrue, charge*exp(drifttime*corr_lambda));
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

  TFile* calibrator_file = TFile::Open((input_dir+"MLL_Calibrator_"+geom_identifier+".root").c_str(), "RECREATE");
  calibrator_file->cd();
  calib_tree->Write();
  h2_QperE_driftTime->Write();
  g_QperE_driftTime->Write();
  h2_Qcorr_Etrue->Write();
  g_Qcorr_Etrue->Write();
  calibrator_file->Close();

  return;
}
