#include <cmath>
#include <iostream>
#include <vector>
#include "RtypesCore.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TTree.h"

#include "flash_matcher.hpp"

const size_t n_combinations = 50;
const double trend_thr = 20.;
// const TString visibility_file_name = "dunevis_fdhd_1x2x6_test_float.root"; // File with the visibility maps


void fm_ana(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLcconfigs f = load_ana_config("./config.json");
  // double trend_thr = f.trend_thr;
  TString visibility_file_name = f.visibility_file_name;


  // --- INPUTS ---------------------------------------------------------------
  TFile* calib_file = TFile::Open("./MLL_Calibrator.root", "READ");
  TTree* calib_tree = static_cast<TTree*>(calib_file->Get("calib_tree"));
  Float_t calib_c = 0.;         Float_t calib_slope = 0.;
  Float_t drift_velocity = 0.0; Float_t corr_lambda = 0.0;
  calib_tree->SetBranchAddress("calib_c", &calib_c); 
  calib_tree->SetBranchAddress("calib_slope", &calib_slope);
  calib_tree->SetBranchAddress("drift_velocity", &drift_velocity);
  calib_tree->SetBranchAddress("corr_lambda", &corr_lambda);
  calib_tree->GetEntry(0);
  // calib_file->Close();

  TFile* distribution_file = TFile::Open("./MLL_Distributions.root", "READ");
  TTree* tpc_pds_tree = static_cast<TTree*>(distribution_file->Get("tpc_pds_tree"));
  int ifile, iev;
  float charge, time_tpc, time_pds, y_reco, z_reco, e_reco;
  float x_true, y_true, z_true, e_true;
  std::vector<float>* reco_pes = nullptr;
  std::vector<float>* exp_phs = nullptr;
  tpc_pds_tree->SetBranchAddress("reco_pes", &reco_pes);
  tpc_pds_tree->SetBranchAddress("exp_phs", &exp_phs);
  tpc_pds_tree->SetBranchAddress("ifile", &ifile);
  tpc_pds_tree->SetBranchAddress("iev", &iev);
  tpc_pds_tree->SetBranchAddress("charge", &charge);
  tpc_pds_tree->SetBranchAddress("time_tpc", &time_tpc);
  tpc_pds_tree->SetBranchAddress("y_reco", &y_reco);
  tpc_pds_tree->SetBranchAddress("z_reco", &z_reco);
  tpc_pds_tree->SetBranchAddress("e_reco", &e_reco);
  tpc_pds_tree->SetBranchAddress("time_pds", &time_pds);
  tpc_pds_tree->SetBranchAddress("x_true", &x_true);
  tpc_pds_tree->SetBranchAddress("y_true", &y_true);
  tpc_pds_tree->SetBranchAddress("z_true", &z_true);
  tpc_pds_tree->SetBranchAddress("e_true", &e_true);

  TTree* pds_tree = static_cast<TTree*>(distribution_file->Get("pds_tree"));
  size_t n_opdet = 0; float LY_times_PDE = 0.;
  pds_tree->SetBranchAddress("n_opdet", &n_opdet);
  pds_tree->SetBranchAddress("LY_times_PDE", &LY_times_PDE);
  pds_tree->GetEntry(0); // Read the first entry to get the parameters
  // TH2D* h2_exp_reco = static_cast<TH2D*>(distribution_file->Get("h2_exp_reco_norm"));

  TFile* parametrizer_file = TFile::Open("./MLL_Parametrizer.root", "READ");
  TF1* f_reco_prob = static_cast<TF1*>(parametrizer_file->Get("f_reco_prob"));
  TF1* f_lognormal = static_cast<TF1*>(parametrizer_file->Get("f_lognormal"));
  TF1* f_logms_trend = static_cast<TF1*>(parametrizer_file->Get("f_logms_trend"));
  TF1* f_sigmas_trend = static_cast<TF1*>(parametrizer_file->Get("f_sigmas_trend"));
  TH2D* h2_exp_reco = static_cast<TH2D*>(parametrizer_file->Get("h2_exp_reco"));

  TGraphErrors* g_logms = static_cast<TGraphErrors*>(parametrizer_file->Get("g_logms")); 
  TGraphErrors* g_sigmas = static_cast<TGraphErrors*>(parametrizer_file->Get("g_sigmas"));

  // --- HISTOS -----------------------------------------------------------------
  TH1D* h_LL = new TH1D("h_LL",Form("%s;%s;%s","h_LL","LL","counts"),
                        200, 0., 600.);
  TH1D* h_LL_fake = new TH1D("h_LL_fake",Form("%s;%s;%s","h_LL_fake","LL_fake","counts"),
                             200, 0., 600.);

  TH1D* h_x_mis = new TH1D("h_x_mis",Form("%s;%s;%s","h_x_mis","x [cm]","counts"),
                          100, -350, 350);

  TH1D* h_y_mis = new TH1D("h_y_mis",Form("%s;%s;%s","h_y_mis","y [cm]","counts"),
                          100, -600, 600);

  TH1D* h_e_mis = new TH1D("h_e_mis",Form("%s;%s;%s","h_e_mis","E [MeV]","counts"),
                          100, 0, 30);

  TH1D* h_z_mis = new TH1D("h_z_mis",Form("%s;%s;%s","h_z_mis","z [cm]","counts"),
                          100, 0, 1400);
  
  TH1D* h_dx_mis = new TH1D("h_dx_mis",Form("%s;%s;%s","h_dx_mis","#Deltax [cm]","counts"),
                          100, -750, 750); 

  TH1D* h_dy_mis = new TH1D("h_dy_mis",Form("%s;%s;%s","h_dy_mis","#Deltay [cm]","counts"),
                          100, -1200, 1200);

  TH1D* h_dz_mis = new TH1D("h_dz_mis",Form("%s;%s;%s","h_dz_mis","#Deltaz [cm]","counts"),
                          100, -1400, 1400);

  TH1D* h_de_mis = new TH1D("h_de_mis",Form("%s;%s;%s","h_de_mis","#DeltaE [MeV]","counts"),
                          100, -10, 10);

  TH2D* h2_de_dx_mis = new TH2D("h2_de_dx_mis",Form("%s;%s;%s","h2_de_dx_mis","#Deltax [cm]","#DeltaE [MeV]"),
                          100, -750, 750,
                          100, -30, 30);

  TH2D* h2_de_dy_mis = new TH2D("h2_de_dy_mis",Form("%s;%s;%s","h2_de_dy_mis","#Deltay [cm]","#DeltaE [MeV]"),
                          100, -1200, 1200,
                          100, -30, 30);

  TH2D* h2_de_dz_mis = new TH2D("h2_de_dz_mis",Form("%s;%s;%s","h2_de_dz_mis","#Deltaz [cm]","#DeltaE [MeV]"),
                          100, -1400, 1400,
                          100, -30, 30);
  
  TH1D* h_dx = new TH1D("h_dx",Form("%s;%s;%s","h_dx","dX_{reco-true} [cm]","Counts"), 200, -100, 100);

  TH2D* h2_xreco_xtrue = new TH2D("h2_xreco_xtrue",
                                  Form("%s;%s;%s", "h2_xreco_xtrue", "X_{true}", "X_{reco}"),
                                  200, 0, 360,
                                  200, 0, 360);

  TH2D* h2_Ereco_Etrue = new TH2D("h2_Ereco_Etrue",
                                  Form("%s;%s;%s", "h2_Ereco_Etrue", "E_{true} [MeV]", "E_{reco} [MeV]"),
                                  200, 0, 30,
                                  200, 0, 30);
  
  TH2D* h2_Ydiff_Zdiff = new TH2D("h2_Ydiff_Zdiff",
                                 Form("%s;%s;%s", "h2_Ydiff_Zdiff", "Y_{True}-Y_{Reco} [cm]", "Z_{True}-Z_{Reco} [cm]"),
                                 200, -650, 650, 200, -650, 650);

  // --- LikelihoodComputer -----------------------------------------------------
  LikelihoodComputer likelihood_computer(
    visibility_file_name, // Visibility file name
    n_opdet,              // Number of optical detectors
    drift_velocity,       // Drift velocity
    LY_times_PDE,         // Light yield times photo detector efficiency
    f_reco_prob,          // Reconstruction probability function
    f_lognormal,          // Lognormal function for extrapolation
    f_logms_trend,        // Trend function for logm
    f_sigmas_trend,       // Trend function for sigma
    g_logms,              // Graph for logm values
    g_sigmas,             // Graph for sigma values
    trend_thr,            // Threshold for trend
    calib_c,              // Calibration constant
    calib_slope,          // Calibration slope
    corr_lambda,          // Correction lambda value
    h2_exp_reco
  );


  // --- LOOP OVER TPC-PDS CLUSTERS ---------------------------------------------
  float ntry = 0; float nmismatch = 0; float ninfinity = 0;
  std::vector<ClusterTPC> fake_tpc_clusters;
  std::vector<VertexInfo> fake_vertex_info;
  for (Long64_t entry = 0; entry < tpc_pds_tree->GetEntries(); entry++) {
    tpc_pds_tree->GetEntry(entry);
    // time_pds = -18;
    // if (time_pds > -13.8) continue;

    if (fake_tpc_clusters.size() < n_combinations) {
      ClusterPDS t_pds_cluster = ClusterPDS(time_pds, reco_pes);
      if (t_pds_cluster.IsEmpty()) {
        std::cerr << "This Error: PDS cluster vectors are empty!" << std::endl;
        continue; // Skip empty clusters
      }
      VertexInfo t_vertex_info = VertexInfo(x_true, y_true, z_true, e_true);
      fake_vertex_info.push_back(t_vertex_info);
      ClusterTPC t_tpc_cluster = ClusterTPC(charge, time_tpc, y_reco, z_reco);
      fake_tpc_clusters.push_back(t_tpc_cluster);
      continue;
    }


    ClusterPDS true_pds_cluster = ClusterPDS(time_pds, reco_pes);
    if (true_pds_cluster.IsEmpty()) {
      // std::cerr << "Error: PDS cluster vectors are empty!" << std::endl;
      continue; // Skip empty clusters
    }
    ClusterTPC true_tpc_cluster = ClusterTPC(charge, time_tpc, y_reco, z_reco);
    VertexInfo true_vertex_info = VertexInfo(x_true, y_true, z_true, e_true);

    float delta_time = true_tpc_cluster.time_tpc - true_pds_cluster.time_pds;
    if (delta_time < 0) {
      std::cout << "-! " << true_tpc_cluster.time_tpc << "-" << std::endl;
    }
    float true_loglikelihood = likelihood_computer.GetLikelihoodMatch(true_tpc_cluster, true_pds_cluster);
    if (std::isinf(true_loglikelihood)) ninfinity++;

    h2_xreco_xtrue->Fill(x_true, likelihood_computer.x_reco);
    h2_Ereco_Etrue->Fill(e_true, likelihood_computer.E_reco);
    h2_Ydiff_Zdiff->Fill(y_reco - y_true, z_reco - z_true);
    h_dx->Fill(likelihood_computer.x_reco - x_true);

    // --- Compute likelihood for each combination of TPC and PDS clusters
    for (size_t i = 0; i < fake_tpc_clusters.size(); i++) {
      // check if E_reco, y_reco, and z_reco are the same as the true cluster, if so skip this combination
      if (std::abs(fake_vertex_info[i].energy - true_vertex_info.energy) < 1e-3 &&
          std::abs(fake_vertex_info[i].y - true_vertex_info.y) < 1e-3 &&
          std::abs(fake_vertex_info[i].z - true_vertex_info.z) < 1e-3) {
        continue; // Skip the combination if it's essentially the same as the true cluster
      }
      delta_time = fake_tpc_clusters[i].time_tpc - true_pds_cluster.time_pds;
      if (delta_time < 0) continue; // Skip combinations where TPC time is before PDS time
      
      float fake_loglikelihood = likelihood_computer.GetLikelihoodMatch(fake_tpc_clusters[i], true_pds_cluster);
      // float fake_loglikelihood = 55; 
      h_LL->Fill(-true_loglikelihood);
      h_LL_fake->Fill(-fake_loglikelihood);
      ntry++;
      if (fake_loglikelihood>true_loglikelihood) {
        nmismatch++;
        h_x_mis->Fill(true_vertex_info.x);
        h_y_mis->Fill(true_vertex_info.y);
        h_z_mis->Fill(true_vertex_info.z);
        h_dx_mis->Fill(true_vertex_info.x - fake_vertex_info[i].x);
        h_dy_mis->Fill(true_vertex_info.y - fake_vertex_info[i].y);
        h_dz_mis->Fill(true_vertex_info.z - fake_vertex_info[i].z);
        h_de_mis->Fill(true_vertex_info.energy - fake_vertex_info[i].energy);
        h2_de_dx_mis->Fill(true_vertex_info.x - fake_vertex_info[i].x, true_vertex_info.energy - fake_vertex_info[i].energy);
        h2_de_dy_mis->Fill(true_vertex_info.y - fake_vertex_info[i].y, true_vertex_info.energy - fake_vertex_info[i].energy);
        h2_de_dz_mis->Fill(true_vertex_info.z - fake_vertex_info[i].z, true_vertex_info.energy - fake_vertex_info[i].energy);
      }
      else if (fake_loglikelihood == true_loglikelihood && !std::isinf(fake_loglikelihood)) {
        std::cout << "Warneeeng: Fake log-likelihood is equal to true log-likelihood!" << std::endl;
        std::cout << fake_loglikelihood << std::endl;
      }
    }

    // --- Rotate the vectors to keep the last element as the current true cluster
    std::rotate(fake_tpc_clusters.begin(), fake_tpc_clusters.begin()+1, fake_tpc_clusters.end());
    std::rotate(fake_vertex_info.begin(), fake_vertex_info.begin()+1, fake_vertex_info.end());
    fake_tpc_clusters[fake_tpc_clusters.size()-1] = true_tpc_cluster;
    fake_vertex_info[fake_vertex_info.size()-1] = true_vertex_info;
  }
  std::cout << "Mismatch: " << nmismatch << "/" << ntry << "\t" << nmismatch/ntry*100 << std::endl;
  std::cout << "of which "  << ninfinity << "/" << nmismatch << " are infinite" << std::endl;

  // --- WRITE OUTPUT ----------------------------------------------------------
  TFile* out_file = TFile::Open("MLL_AnaOutput.root", "RECREATE");
  out_file->cd();
  h_LL->Write();
  h_LL_fake->Write();
  h_dx->Write();
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
  h2_xreco_xtrue->Write();
  h2_Ereco_Etrue->Write();
  h2_Ydiff_Zdiff->Write();
  h2_exp_reco->Write();
  likelihood_computer.g_logms->Write();
  out_file->Close();

  return;
}
