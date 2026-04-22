#include <cmath>
#include <iostream>
#include <vector>
#include "RtypesCore.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TString.h"
#include "TTree.h"

#include "flash_matcher.hpp"



void fm_ana(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLConfigs mll_conf = load_ana_config("./configs/ana_config.json");
  std::string sample_config_file = mll_conf.sample_config_file;
  float light_yield            = mll_conf.light_yield;
  float arapuca_pde            = mll_conf.arapuca_pde;
  std::string visibility_dir   = mll_conf.visibility_dir;
  size_t n_combinations        = mll_conf.n_combinations;
  bool loop_on_tpc_clusters    = mll_conf.loop_on_tpc_clusters;
  
  SampleConfigs sample_conf = load_sample_config("./configs/"+sample_config_file);
  std::string input_dir        = sample_conf.input_dir;
  std::string geom_identifier  = sample_conf.geom_identifier;
  double trend_thr             = sample_conf.trend_thr;
  double q_cut_low             = sample_conf.q_cut_low;
  float LY_times_PDE           = light_yield * arapuca_pde;

  DuneGeom geom = load_dune_geom("./configs/"+geom_identifier+".json");
  

  // --- INPUTS ---------------------------------------------------------------
  TFile* calib_file = TFile::Open((input_dir+"MLL_Calibrator_"+geom_identifier+".root").c_str(), "READ");
  TTree* calib_tree = static_cast<TTree*>(calib_file->Get("calib_tree"));
  Float_t calib_c = 0.;         Float_t calib_slope = 0.;
  Float_t drift_velocity = 0.0; Float_t corr_lambda = 0.0;
  calib_tree->SetBranchAddress("calib_c", &calib_c); 
  calib_tree->SetBranchAddress("calib_slope", &calib_slope);
  calib_tree->SetBranchAddress("drift_velocity", &drift_velocity);
  calib_tree->SetBranchAddress("corr_lambda", &corr_lambda);
  calib_tree->GetEntry(0);

  TFile* distribution_file = TFile::Open((input_dir+"MLL_Distributions_"+geom_identifier+".root").c_str(), "READ");
  TTree* tpc_pds_tree = static_cast<TTree*>(distribution_file->Get("tpc_pds_tree"));
  std::vector<size_t> MaxChargeIndxs = take_max_charge_indices(tpc_pds_tree, "iev", "charge");
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

  TFile* parametrizer_file  = TFile::Open((input_dir+"MLL_Parametrizer_"+geom_identifier+".root").c_str(), "READ");
  TF1* f_RecoExpDistr       = static_cast<TF1*>(parametrizer_file->Get("f_RecoExpDistr"));
  TF1* f_par1_trend         = static_cast<TF1*>(parametrizer_file->Get("f_par1_trend"));
  TF1* f_par2_trend         = static_cast<TF1*>(parametrizer_file->Get("f_par2_trend"));
  TH2D* h2_exp_reco         = static_cast<TH2D*>(parametrizer_file->Get("h2_exp_reco"));
  TGraphErrors* g_par1      = static_cast<TGraphErrors*>(parametrizer_file->Get("g_par1")); 
  TGraphErrors* g_par2      = static_cast<TGraphErrors*>(parametrizer_file->Get("g_par2"));

  TEfficiency* he_hit_prob = static_cast<TEfficiency*>(parametrizer_file->Get("he_hit_prob"));

  // --- HISTOS -----------------------------------------------------------------
  std::string loop_string = (loop_on_tpc_clusters) ? "LoopOnTPCClusters" : "LoopOnPDSClusters";
  TFile* out_file = TFile::Open((input_dir+"MLL_AnaOutput_"+geom_identifier+"_"+loop_string+".root").c_str(), "RECREATE");
  float t_x_true, t_y_true, t_z_true, t_e_true, t_x_reco, t_y_reco, t_z_reco, t_e_reco, t_NLL, t_nhits;
  float f_x_true, f_y_true, f_z_true, f_e_true, f_x_reco, f_y_reco, f_z_reco, f_e_reco, f_NLL, f_nhits;
  std::vector<float>* t_reco_pes = nullptr; std::vector<float>* t_exp_phs = nullptr; std::vector<float>* t_reco_terms = nullptr; std::vector<float>* t_noreco_terms = nullptr;
  std::vector<float>* f_reco_pes = nullptr; std::vector<float>* f_exp_phs = nullptr; std::vector<float>* f_reco_terms = nullptr; std::vector<float>* f_noreco_terms = nullptr;
  // bool match = false;

  TTree* mismatch_tree = new TTree("mismatch_tree", "mismatch_tree");
  mismatch_tree->Branch("t_x_true", &t_x_true);
  mismatch_tree->Branch("f_x_true", &f_x_true);
  mismatch_tree->Branch("t_y_true", &t_y_true);
  mismatch_tree->Branch("f_y_true", &f_y_true);
  mismatch_tree->Branch("t_z_true", &t_z_true);
  mismatch_tree->Branch("f_z_true", &f_z_true);
  mismatch_tree->Branch("t_e_true", &t_e_true);
  mismatch_tree->Branch("f_e_true", &f_e_true);
  mismatch_tree->Branch("t_x_reco", &t_x_reco);
  mismatch_tree->Branch("f_x_reco", &f_x_reco);
  mismatch_tree->Branch("t_y_reco", &t_y_reco);
  mismatch_tree->Branch("f_y_reco", &f_y_reco);
  mismatch_tree->Branch("t_z_reco", &t_z_reco);
  mismatch_tree->Branch("f_z_reco", &f_z_reco);
  mismatch_tree->Branch("t_e_reco", &t_e_reco);
  mismatch_tree->Branch("f_e_reco", &f_e_reco);
  mismatch_tree->Branch("t_NLL", &t_NLL);
  mismatch_tree->Branch("f_NLL", &f_NLL);
  mismatch_tree->Branch("t_nhits", &t_nhits);
  mismatch_tree->Branch("f_nhits", &f_nhits);
  mismatch_tree->Branch("t_reco_pes", &t_reco_pes);
  mismatch_tree->Branch("f_reco_pes", &f_reco_pes);
  mismatch_tree->Branch("t_exp_phs", &t_exp_phs);
  mismatch_tree->Branch("f_exp_phs", &f_exp_phs);
  mismatch_tree->Branch("t_reco_terms", &t_reco_terms);
  mismatch_tree->Branch("f_reco_terms", &f_reco_terms);
  mismatch_tree->Branch("t_noreco_terms", &t_noreco_terms);
  mismatch_tree->Branch("f_noreco_terms", &f_noreco_terms);
  // mismatch_tree->Branch("match", &match);

  TH1D* h_nmismatch = new TH1D("h_nmismatch",Form("%s;%s;%s","h_nmismatch","#Mismatch","Counts"), n_combinations, -0.5, float(n_combinations)-0.5);

  TH1D* h_TrueRecoTerms = new TH1D("h_TrueRecoTerms",Form("%s;%s;%s","h_TrueRecoTerms","Reco Terms for True Match","counts"),
                                    200, 0, 100);

  TH1D* h_TrueNoRecoTerms = new TH1D("h_TrueNoRecoTerms",Form("%s;%s;%s","h_TrueNoRecoTerms","NoReco Terms for True Match","counts"),
                                    200, 0, 100);

  TH1D* h_FakeRecoTerms = new TH1D("h_FakeRecoTerms",Form("%s;%s;%s","h_FakeRecoTerms","Reco Terms for Fake Match","counts"),
                                    200, 0, 100);

  TH1D* h_FakeNoRecoTerms = new TH1D("h_FakeNoRecoTerms",Form("%s;%s;%s","h_FakeNoRecoTerms","NoReco Terms for Fake Match","counts"),
                                    200, 0, 100);

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

  TH1D* h_ErecoEtrue_mis = new TH1D("h_ErecoEtrue_mis",Form("%s;%s;%s","h_ErecoEtrue_mis","E_{reco}/E_{true}","counts"),
                                100, 0, 2);

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
                                  200, -360, 360,
                                  200, -360, 360);

  TH2D* h2_Ereco_Etrue = new TH2D("h2_Ereco_Etrue",
                                  Form("%s;%s;%s", "h2_Ereco_Etrue", "E_{true} [MeV]", "E_{reco} [MeV]"),
                                  200, 0, 30,
                                  200, 0, 30);
  
  TH2D* h2_Ydiff_Zdiff = new TH2D("h2_Ydiff_Zdiff",
                                 Form("%s;%s;%s", "h2_Ydiff_Zdiff", "Y_{True}-Y_{Reco} [cm]", "Z_{True}-Z_{Reco} [cm]"),
                                 200, -650, 650, 200, -650, 650);

  // --- LikelihoodComputer -----------------------------------------------------
  TString visibility_file_name = TString(visibility_dir+"dunevis_"+geom_identifier+".root");
  LikelihoodComputer likelihood_computer(
    visibility_file_name, // Visibility file name
    geom,                 // DUNE geometry
    drift_velocity,       // Drift velocity
    LY_times_PDE,         // Light yield times photo detector efficiency
    he_hit_prob,          // Hit probability function (TEfficiency)
    f_RecoExpDistr,       // PDF for extrapolation
    f_par1_trend,         // Trend function for par1
    f_par2_trend,         // Trend function for par2
    g_par1,               // Graph for par1 values
    g_par2,               // Graph for par2 values
    trend_thr,            // Threshold for trend
    calib_c,              // Calibration constant
    calib_slope,          // Calibration slope
    corr_lambda,          // Correction lambda value
    h2_exp_reco
  );


  // --- LOOP OVER TPC-PDS CLUSTERS ---------------------------------------------
  float ntry = 0; float nmismatch = 0; float ninfinity = 0; float n_equal = 0;
  float x_sign = 0.;
  std::vector<float> LLs_true, LLs_fake, LLs_true_scaled;
  std::vector<ClusterPDS> fake_pds_clusters;
  std::vector<ClusterTPC> fake_tpc_clusters;
  std::vector<VertexInfo> fake_vertex_infos;
  // for (Long64_t entry = 0; entry < tpc_pds_tree->GetEntries(); entry++) {
  // for (Long64_t entry = 0; entry < 10000; entry++) {
  for (size_t entry : MaxChargeIndxs) {
    if (entry % 100 == 0) std::cout <<entry<<"/"<< MaxChargeIndxs.size()<< "\r" << std::flush;
    // if (entry>10000) break;
    tpc_pds_tree->GetEntry(entry);
    x_sign = (geom_identifier == "dune10kt" && x_true <= 0) ? -1. : 1.;
    // time_pds = -18;
    // if (time_pds > -13.8) continue;

    if (fake_tpc_clusters.size() < n_combinations) {
      ClusterPDS tmp_pds_cluster = ClusterPDS(time_pds, *reco_pes);
      if (tmp_pds_cluster.IsEmpty()) {
        std::cerr << "This Error: PDS cluster vectors are empty!" << std::endl;
        continue; // Skip empty clusters
      }
      VertexInfo tmp_vertex_info = VertexInfo(x_true, y_true, z_true, e_true);
      ClusterTPC tmp_tpc_cluster = ClusterTPC(charge, time_tpc, y_reco, z_reco);
      fake_pds_clusters.push_back(tmp_pds_cluster);
      fake_tpc_clusters.push_back(tmp_tpc_cluster);
      fake_vertex_infos.push_back(tmp_vertex_info);
      continue;
    }


    ClusterPDS true_pds_cluster = ClusterPDS(time_pds, *reco_pes);
    if (true_pds_cluster.IsEmpty() || charge < q_cut_low) {
      // std::cerr << "Error: PDS cluster vectors are empty!" << std::endl;
      continue; // Skip empty clusters
    }
    ClusterTPC true_tpc_cluster = ClusterTPC(charge, time_tpc, y_reco, z_reco);
    VertexInfo true_vertex_info = VertexInfo(x_true, y_true, z_true, e_true);


    float delta_time = true_tpc_cluster.time_tpc - true_pds_cluster.time_pds;
    if (delta_time < 0) {
      std::cout << "-! " << true_tpc_cluster.time_tpc << "-" << std::endl;
    }

    std::vector<float> true_reco_terms, true_noreco_terms;
    float true_NLL = likelihood_computer.GetLikelihoodMatch(true_tpc_cluster, true_pds_cluster, true_reco_terms, true_noreco_terms, x_sign);
    // std::cout << "vvvv " << -true_NLL << std::endl;
    // push back the true log-likelihood to the vector if not inf or nan
    if (!std::isinf(true_NLL) && !std::isnan(true_NLL)) LLs_true.push_back(true_NLL);
    for (auto& term : true_reco_terms) h_TrueRecoTerms->Fill(term);
    for (auto& term : true_noreco_terms) h_TrueNoRecoTerms->Fill(term);
    float e_reco = likelihood_computer.E_reco;
    if (std::isinf(true_NLL)) ninfinity++;

    h2_xreco_xtrue->Fill(x_true, likelihood_computer.x_reco);
    h2_Ereco_Etrue->Fill(e_true, likelihood_computer.E_reco);
    h2_Ydiff_Zdiff->Fill(y_reco - y_true, z_reco - z_true);
    h_dx->Fill(likelihood_computer.x_reco - x_true);

    // Set mismatch_tree variables for the true match
    t_x_true = x_true; t_y_true = y_true; t_z_true = z_true; t_e_true = e_true;
    t_y_reco = y_reco; t_z_reco = z_reco; t_e_reco = e_reco;
    std::vector<float> tmp_exp_phs = likelihood_computer.exp_phs;
    t_reco_pes = reco_pes; t_exp_phs = &tmp_exp_phs;
    t_x_reco = likelihood_computer.x_reco; t_NLL = true_NLL; t_nhits = likelihood_computer.n_hit;
    t_reco_terms = &true_reco_terms; t_noreco_terms = &true_noreco_terms;

    // --- Compute likelihood for each combination of TPC and PDS clusters
    int n_mismatch = 0;
    for (size_t i = 0; i < fake_tpc_clusters.size(); i++) {
      // check if E_reco, y_reco, and z_reco are the same as the true cluster, if so skip this combination
      if (std::abs(fake_vertex_infos[i].energy - true_vertex_info.energy) < 1e-3 &&
          std::abs(fake_vertex_infos[i].y - true_vertex_info.y) < 1e-3 &&
          std::abs(fake_vertex_infos[i].z - true_vertex_info.z) < 1e-3) {
        continue; // Skip the combination if it's essentially the same as the true cluster
      }
      delta_time = (loop_on_tpc_clusters) ? fake_tpc_clusters[i].time_tpc - true_pds_cluster.time_pds :
                                            true_tpc_cluster.time_tpc - fake_pds_clusters[i].time_pds;
      if (delta_time < 0) continue; // Skip combinations where TPC time is before PDS time
      
      std::vector<float> fake_reco_terms, fake_noreco_terms;
      float fake_NLL = (loop_on_tpc_clusters) ? likelihood_computer.GetLikelihoodMatch(fake_tpc_clusters[i], true_pds_cluster, fake_reco_terms, fake_noreco_terms, x_sign) :
                                                likelihood_computer.GetLikelihoodMatch(true_tpc_cluster, fake_pds_clusters[i], fake_reco_terms, fake_noreco_terms, x_sign);

      
      LLs_true_scaled.push_back(true_NLL);
      LLs_fake.push_back(fake_NLL);
      
      for (auto& term : fake_reco_terms) h_FakeRecoTerms->Fill(term);
      for (auto& term : fake_noreco_terms) h_FakeNoRecoTerms->Fill(term);
      ntry++;
      if (fake_NLL<true_NLL) {
        nmismatch++;
        n_mismatch++;
        h_x_mis->Fill(true_vertex_info.x);
        h_y_mis->Fill(true_vertex_info.y);
        h_z_mis->Fill(true_vertex_info.z);
        h_dx_mis->Fill(true_vertex_info.x - fake_vertex_infos[i].x);
        h_dy_mis->Fill(true_vertex_info.y - fake_vertex_infos[i].y);
        h_dz_mis->Fill(true_vertex_info.z - fake_vertex_infos[i].z);
        h_de_mis->Fill(true_vertex_info.energy - fake_vertex_infos[i].energy);
        h_ErecoEtrue_mis->Fill(e_reco/true_vertex_info.energy);
        h2_de_dx_mis->Fill(true_vertex_info.x - fake_vertex_infos[i].x, true_vertex_info.energy - fake_vertex_infos[i].energy);
        h2_de_dy_mis->Fill(true_vertex_info.y - fake_vertex_infos[i].y, true_vertex_info.energy - fake_vertex_infos[i].energy);
        h2_de_dz_mis->Fill(true_vertex_info.z - fake_vertex_infos[i].z, true_vertex_info.energy - fake_vertex_infos[i].energy);
      
        // Set mismatch_tree variables for the fake match
        f_x_true = fake_vertex_infos[i].x; f_y_true = fake_vertex_infos[i].y; f_z_true = fake_vertex_infos[i].z; f_e_true = fake_vertex_infos[i].energy;
        f_x_reco = likelihood_computer.x_reco; f_y_reco = fake_tpc_clusters[i].y_reco; f_z_reco = fake_tpc_clusters[i].z_reco;
        f_e_reco = likelihood_computer.E_reco; f_NLL = fake_NLL; f_nhits = likelihood_computer.n_hit;
        f_reco_pes = reco_pes; f_exp_phs = &likelihood_computer.exp_phs;
        f_reco_terms = &fake_reco_terms; f_noreco_terms = &fake_noreco_terms;
        mismatch_tree->Fill();
      }
      else if (fake_NLL == true_NLL && !std::isinf(fake_NLL)) {
        n_equal++;
        std::cout << "Warneeeng: Fake log-likelihood is equal to true log-likelihood!" << std::endl;
        std::cout << fake_NLL << std::endl;
      }
      h_nmismatch->Fill(n_mismatch);
    }

    // --- Rotate the vectors to keep the last element as the current true cluster
    std::rotate(fake_pds_clusters.begin(), fake_pds_clusters.begin()+1, fake_pds_clusters.end());
    std::rotate(fake_tpc_clusters.begin(), fake_tpc_clusters.begin()+1, fake_tpc_clusters.end());
    std::rotate(fake_vertex_infos.begin(), fake_vertex_infos.begin()+1, fake_vertex_infos.end());
    fake_pds_clusters[fake_pds_clusters.size()-1] = true_pds_cluster;
    fake_tpc_clusters[fake_tpc_clusters.size()-1] = true_tpc_cluster;
    fake_vertex_infos[fake_vertex_infos.size()-1] = true_vertex_info;
  }

  // sort LLs_true
  std::sort(LLs_true.begin(), LLs_true.end());
  double hLL_min = LLs_true[0];
  double hLL_max = LLs_true[size_t(LLs_true.size()*0.97)];

  TH1D* h_LL       = new TH1D("h_LL",Form("%s;%s;%s","h_LL","LL","counts"), 200, hLL_min, hLL_max);
  TH1D* h_LL_scale = new TH1D("h_LL_scale",Form("%s;%s;%s","h_LL_scale","LL","counts"), 200, hLL_min, hLL_max);
  TH1D* h_LL_fake  = new TH1D("h_LL_fake",Form("%s;%s;%s","h_LL_fake","LL_fake","counts"), 200, hLL_min, hLL_max);

  for (auto& ll : LLs_true) h_LL->Fill(ll);
  for (auto& ll : LLs_true_scaled) h_LL_scale->Fill(ll);
  for (auto& ll : LLs_fake) h_LL_fake->Fill(ll);
  LLs_true.clear(); LLs_true_scaled.clear(); LLs_fake.clear();

  std::cout << "N equals: " << n_equal << "/" << ntry << std::endl;
  std::cout << "Mismatch: " << nmismatch << "/" << ntry << "\t" << nmismatch/ntry*100 << std::endl;
  std::cout << "of which "  << ninfinity << "/" << nmismatch << " are infinite" << std::endl;
  std::printf("True matches : %.0f (%.3f%%) of reco_terms and %.0f (%.3f%%) of non-reco\n",
              h_TrueRecoTerms->GetEntries(), h_TrueRecoTerms->GetEntries()/(geom.n_opdet*h_LL->GetEntries())*100,
              h_TrueNoRecoTerms->GetEntries(), h_TrueNoRecoTerms->GetEntries()/(geom.n_opdet*h_LL->GetEntries())*100);
  std::printf("False matches %.0f (%.3f%%) of reco_terms and %.0f (%.3f%%) of non-reco\n",
              h_FakeRecoTerms->GetEntries(), h_FakeRecoTerms->GetEntries()/(geom.n_opdet*h_LL_fake->GetEntries())*100,
              h_FakeNoRecoTerms->GetEntries(), h_FakeNoRecoTerms->GetEntries()/(geom.n_opdet*h_LL_fake->GetEntries())*100);
  

  // --- WRITE OUTPUT ----------------------------------------------------------
  out_file->cd();
  mismatch_tree->Write("", TObject::kOverwrite);
  h_LL->Write();
  h_LL_scale->Write();
  h_LL_fake->Write();
  h_TrueRecoTerms->Write();
  h_TrueNoRecoTerms->Write();
  h_FakeRecoTerms->Write();
  h_FakeNoRecoTerms->Write();
  h_dx->Write();
  h_x_mis->Write();
  h_y_mis->Write();
  h_z_mis->Write();
  h_dx_mis->Write();
  h_dy_mis->Write();
  g_par1->Write();
  g_par2->Write();
  h_dz_mis->Write();
  h_de_mis->Write();
  h_ErecoEtrue_mis->Write();
  h2_de_dx_mis->Write();
  h2_de_dy_mis->Write();
  h2_de_dz_mis->Write();
  h2_xreco_xtrue->Write();
  h2_Ereco_Etrue->Write();
  h2_Ydiff_Zdiff->Write();
  h2_exp_reco->Write();
  h_nmismatch->Write();
  out_file->Close();

  calib_file->Close();
  
  return;
}
