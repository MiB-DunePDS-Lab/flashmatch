#include "RtypesCore.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TTree.h"

#include "TTreeReader.h"
#include "Utils.hpp"
#include <iostream>
#include <vector>

const size_t n_combinations = 5;
const double trend_thr = 20.;
const TString visibility_file_name = "dunevis_fdhd_1x2x6_test_float.root"; // File with the visibility maps

class ClusterTPC{
  public:
    float charge;
    float time_tpc;
    float y_reco;
    float z_reco;

    ClusterTPC& operator=(const ClusterTPC& other) {
      if (this != &other) {
        charge = other.charge;
        time_tpc = other.time_tpc;
        y_reco = other.y_reco;
        z_reco = other.z_reco;
      }
      return *this;
    }

    bool operator==(const ClusterTPC& other) {
      return (charge == other.charge &&
          time_tpc == other.time_tpc &&
          y_reco == other.y_reco &&
          z_reco == other.z_reco);
    }

    ClusterTPC() : charge(0), time_tpc(0), y_reco(0), z_reco(0) {}

    ClusterTPC(const float charge,
        const float time_tpc,
        const float y_reco,
        const float z_reco) {
      this->charge = charge;
      this->time_tpc = time_tpc;
      this->y_reco = y_reco;
      this->z_reco = z_reco;    
    }

};

class ClusterPDS{
  public:
    float time_pds;
    std::vector<size_t>* flash_channels;
    std::vector<float>* exp_phs;
    std::vector<float>* reco_pes;

    ClusterPDS& operator=(const ClusterPDS& other) {
      if (this != &other) {
        time_pds = other.time_pds;
        flash_channels = other.flash_channels;
        exp_phs = other.exp_phs;
        reco_pes = other.reco_pes;
      }
      return *this;
    }

    bool operator==(const ClusterPDS& other) {
      return (time_pds == other.time_pds &&
          flash_channels == other.flash_channels &&
          exp_phs == other.exp_phs &&
          reco_pes == other.reco_pes &&
          *flash_channels == *other.flash_channels &&
          *exp_phs == *other.exp_phs &&
          *reco_pes == *other.reco_pes);
    }

    ClusterPDS() : time_pds(0), flash_channels(nullptr), exp_phs(nullptr), reco_pes(nullptr) {
      flash_channels = new std::vector<size_t>();
      exp_phs = new std::vector<float>();
      reco_pes = new std::vector<float>();
    }


    ClusterPDS(const float& time_pds,
        std::vector<size_t>* flash_channels,
        std::vector<float>* exp_phs,
        std::vector<float>* reco_pes) {
      this->time_pds = time_pds;
      this->flash_channels = flash_channels;
      this->exp_phs = exp_phs;
      this->reco_pes = reco_pes;
    }

    bool IsEmpty(){
      return (flash_channels->empty() || exp_phs->empty() || reco_pes->empty());
    }
};

class LikelihoodComputer{

  public:

    float E_reco;
    float x_reco;
    TGraphErrors* g_logms;

    float GetLikelihoodMatch(const ClusterTPC& tpc_cluster,
        const ClusterPDS& pds_cluster){
      // float delta_time = pds_cluster.time_pds-tpc_cluster.time_tpc;
      float delta_time = tpc_cluster.time_tpc;
      E_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, delta_time, tpc_cluster.charge);
      float reco_pe, exp_ph;
      x_reco = delta_time*drift_velocity;
      float vertex_coor[3] = {x_reco, tpc_cluster.y_reco, tpc_cluster.z_reco};
      int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);

      float log_likelihood = 0.;

      for (size_t idx_opdet=0; idx_opdet<pds_cluster.reco_pes->size(); idx_opdet++){
        float voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];
        exp_ph = E_reco*LY_times_PDE*voxel_vis;
        if(exp_ph==0) exp_ph = E_reco*LY_times_PDE*1.e-15;
        float P_hit_mu = f_reco_prob->Eval(exp_ph);

        reco_pe = pds_cluster.reco_pes->at(idx_opdet);

        if (reco_pe>0.0){
          if (reco_pe > trend_thr){
            f_lognormal->SetParameters(f_logms_trend->Eval(reco_pe), f_sigmas_trend->Eval(reco_pe));
          } else {
            f_lognormal->SetParameters(g_logms->Eval(reco_pe), g_sigmas->Eval(reco_pe));
          }
          // log_likelihood += log(P_hit_mu*f_lognormal->Eval(reco_pe));
          log_likelihood += log(P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(reco_pe, exp_ph)));
        } else {
          log_likelihood += log(1. - P_hit_mu);
        }
      }

      return log_likelihood;
    } // GetLikelihoodMatch

    // LikelihoodComputer constructor
    LikelihoodComputer(TString visibility_file_name,
        size_t n_opdet,
        float drift_velocity,
        float LY_times_PDE,
        TF1* f_reco_prob,
        TF1* f_lognormal,
        TF1* f_logms_trend,
        TF1* f_sigmas_trend,
        TGraphErrors* g_logms,
        TGraphErrors* g_sigmas,
        double trend_thr,
        float calib_c, 
        float calib_slope, 
        float corr_lambda,
        TH2D* h2_exp_reco = nullptr)
      : visibility_file_name(visibility_file_name),
      n_opdet(n_opdet),
      drift_velocity(drift_velocity),
      LY_times_PDE(LY_times_PDE),
      f_reco_prob(f_reco_prob),
      f_lognormal(f_lognormal),
      f_logms_trend(f_logms_trend),
      f_sigmas_trend(f_sigmas_trend),
      g_logms(g_logms),
      g_sigmas(g_sigmas),
      trend_thr(trend_thr),
      calib_c(calib_c), 
      calib_slope(calib_slope), 
      corr_lambda(corr_lambda) {
        this->h2_exp_reco = h2_exp_reco; // Optional, can be nullptr
        setptivatemembers();
      } // LikelihoodComputer constructor

  private:
    // take in input
    TString visibility_file_name;
    size_t n_opdet;
    float drift_velocity;
    float LY_times_PDE;
    std::vector<std::vector<float>> opDet_visMapDirect;
    TF1* f_reco_prob;
    TF1* f_lognormal;
    TF1* f_logms_trend;
    TF1* f_sigmas_trend;
    double trend_thr;
    TGraphErrors* g_sigmas;
    float calib_c, calib_slope, corr_lambda;
    TH2D* h2_exp_reco = nullptr; // Optional, can be nullptr

    // set inside the class
    TH1D* hgrid[3] = {nullptr};
    std::vector<int> cryo_to_tpc;

    void setptivatemembers(){
      TFile* visibility_file = TFile::Open(visibility_file_name, "READ");
      if (!visibility_file || visibility_file->IsZombie()) {
        std::cerr << "Error opening visibility file: " << visibility_file_name << std::endl;
        return;
      }
      hgrid[0] = (TH1D*)visibility_file->Get("photovisAr/hgrid0")->Clone("hclone_hgrid0");
      hgrid[1] = (TH1D*)visibility_file->Get("photovisAr/hgrid1")->Clone("hclone_hgrid1");
      hgrid[2] = (TH1D*)visibility_file->Get("photovisAr/hgrid2")->Clone("hclone_hgrid2");

      TTree* tDimensions = (TTree*)visibility_file->Get("photovisAr/tDimensions");
      Double_t coor_dim[3] = {0.};
      tDimensions->SetBranchAddress("dimension", coor_dim);
      tDimensions->GetEntry(0);
      float tpc_min[3] = {float(coor_dim[0]), float(coor_dim[1]), float(coor_dim[2])}; // x,y,z
      tDimensions->GetEntry(1);
      float tpc_max[3] = {float(coor_dim[0]), float(coor_dim[1]), float(coor_dim[2])}; // x,y,z
      cryo_to_tpc = GetCryoToTPCMap(hgrid, tpc_min, tpc_max);

      TTree* photoVisMap = (TTree*)visibility_file->Get("photovisAr/photoVisMap");
      const size_t n_entriesmap = photoVisMap->GetEntries();
      float opDet_visDirect[n_opdet];
      photoVisMap->SetBranchAddress("opDet_visDirect", &opDet_visDirect);
      opDet_visMapDirect.resize(n_entriesmap, std::vector<float>(n_opdet, 0.));

      for (size_t VisMapEntry = 0; VisMapEntry < n_entriesmap; ++VisMapEntry) {
        photoVisMap->GetEntry(VisMapEntry);
        std::memcpy(opDet_visMapDirect[VisMapEntry].data(), &opDet_visDirect[0], n_opdet * sizeof(float));
      }

      return;
    } // 

}; // LikelihoodComputer


void fm_ana(){
  TFile* calib_file = TFile::Open("./calibrator.root", "READ");
  TTree* calib_tree = static_cast<TTree*>(calib_file->Get("calib_tree"));
  Float_t calib_c = 0.;         Float_t calib_slope = 0.;
  Float_t drift_velocity = 0.0; Float_t corr_lambda = 0.0;
  calib_tree->SetBranchAddress("calib_c", &calib_c); 
  calib_tree->SetBranchAddress("calib_slope", &calib_slope);
  calib_tree->SetBranchAddress("drift_velocity", &drift_velocity);
  calib_tree->SetBranchAddress("corr_lambda", &corr_lambda);
  calib_tree->GetEntry(0);
  // calib_file->Close();

  TFile* distribution_file = TFile::Open("fm_distributions.root", "READ");
  TTree* tpc_pds_tree = static_cast<TTree*>(distribution_file->Get("tpc_pds_tree"));
  int ifile, iev;
  float charge, time_tpc, time_pds, y_reco, z_reco, e_reco;
  float x_true, y_true, z_true, e_true;
  std::vector<size_t>* flash_channels = nullptr;
  std::vector<float>* flash_pes = nullptr;
  std::vector<float>* flash_times = nullptr;
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
  tpc_pds_tree->SetBranchAddress("flash_channels", &flash_channels);
  tpc_pds_tree->SetBranchAddress("flash_pes", &flash_pes);
  tpc_pds_tree->SetBranchAddress("flash_times", &flash_times);
  tpc_pds_tree->SetBranchAddress("x_true", &x_true);
  tpc_pds_tree->SetBranchAddress("y_true", &y_true);
  tpc_pds_tree->SetBranchAddress("z_true", &z_true);
  tpc_pds_tree->SetBranchAddress("e_true", &e_true);

  TTree* pds_tree = static_cast<TTree*>(distribution_file->Get("pds_tree"));
  int n_opdet = 0; float LY_times_PDE = 0.;
  pds_tree->SetBranchAddress("n_opdet", &n_opdet);
  pds_tree->SetBranchAddress("LY_times_PDE", &LY_times_PDE);
  pds_tree->GetEntry(0); // Read the first entry to get the parameters
  TH2D* h2_exp_reco = static_cast<TH2D*>(distribution_file->Get("h2_exp_reco_norm"));

  TFile* parametrizer_file = TFile::Open("fm_parametrizer.root", "READ");
  TF1* f_reco_prob = static_cast<TF1*>(parametrizer_file->Get("f_reco_prob"));
  TF1* f_lognormal = static_cast<TF1*>(parametrizer_file->Get("f_lognormal"));
  TF1* f_logms_trend = static_cast<TF1*>(parametrizer_file->Get("f_logms_trend"));
  TF1* f_sigmas_trend = static_cast<TF1*>(parametrizer_file->Get("f_sigmas_trend"));
  // TH2D* h2_exp_reco = static_cast<TH2D*>(parametrizer_file->Get("h2_exp_reco"));

  TGraphErrors* g_logms = static_cast<TGraphErrors*>(parametrizer_file->Get("g_logms")); 
  TGraphErrors* g_sigmas = static_cast<TGraphErrors*>(parametrizer_file->Get("g_sigmas"));

  // --- HISTOS -----------------------------------------------------------------
  TH1D* h_LL = new TH1D("h_LL",Form("%s;%s;%s","h_LL","LL","counts"),
                        200, 0., 600.);
  TH1D* h_LL_fake = new TH1D("h_LL_fake",Form("%s;%s;%s","h_LL_fake","LL_fake","counts"),
                             200, 0., 600.);

  TH1D* h_dx = new TH1D("h_dx",Form("%s;%s;%s","h_dx","dX","Counts"), 200, -100, 100);

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
    visibility_file_name,      // Visibility file name
    n_opdet,                // Number of optical detectors
    drift_velocity,         // Drift velocity
    LY_times_PDE,          // Light yield times photo detector efficiency
    f_reco_prob,           // Reconstruction probability function
    f_lognormal,           // Lognormal function for extrapolation
    f_logms_trend,         // Trend function for logm
    f_sigmas_trend,        // Trend function for sigma
    g_logms,               // Graph for logm values
    g_sigmas,              // Graph for sigma values
    trend_thr,             // Threshold for trend
    calib_c,               // Calibration constant
    calib_slope,           // Calibration slope
    corr_lambda,            // Correction lambda value
    h2_exp_reco
  );


  // --- LOOP OVER TPC-PDS CLUSTERS ---------------------------------------------
  float ntry = 0; float nmismatch = 0;
  std::vector<ClusterTPC> fake_tpc_clusters;
  for (Long64_t entry = 0; entry < tpc_pds_tree->GetEntries(); entry++) {
    tpc_pds_tree->GetEntry(entry);

    if (fake_tpc_clusters.size() < n_combinations) {
      ClusterTPC t_tpc_cluster = ClusterTPC(charge, time_tpc, y_reco, z_reco);
      ClusterPDS t_pds_cluster = ClusterPDS(time_pds, flash_channels, exp_phs, reco_pes);
      if (t_pds_cluster.IsEmpty()) {
        std::cerr << "This Error: PDS cluster vectors are empty!" << std::endl;
        continue; // Skip empty clusters
      }
      fake_tpc_clusters.push_back(t_tpc_cluster);
      continue;
    }



    ClusterPDS true_pds_cluster = ClusterPDS(time_pds, flash_channels, exp_phs, reco_pes);
    if (true_pds_cluster.IsEmpty()) {
      // std::cerr << "Error: PDS cluster vectors are empty!" << std::endl;
      continue; // Skip empty clusters
    }
    ClusterTPC true_tpc_cluster = ClusterTPC(charge, time_tpc, y_reco, z_reco);

    float true_loglikelihood = likelihood_computer.GetLikelihoodMatch(true_tpc_cluster,     true_pds_cluster);

    h2_xreco_xtrue->Fill(x_true, likelihood_computer.x_reco);
    h2_Ereco_Etrue->Fill(e_true, likelihood_computer.E_reco);
    h2_Ydiff_Zdiff->Fill(y_reco - y_true, z_reco - z_true);
    h_dx->Fill(likelihood_computer.x_reco - x_true);

    // --- Compute likelihood for each combination of TPC and PDS clusters
    for (size_t i = 0; i < fake_tpc_clusters.size(); i++) {
      float fake_loglikelihood = likelihood_computer.GetLikelihoodMatch(fake_tpc_clusters[i], true_pds_cluster);
      // float fake_loglikelihood = -5; 
      h_LL->Fill(-true_loglikelihood);
      h_LL_fake->Fill(-fake_loglikelihood);
      ntry++;
      if (fake_loglikelihood>true_loglikelihood) nmismatch++;
    }

  }
  std::cout << "Mismatch: " << nmismatch << "/" << ntry << "\t" << nmismatch/ntry*100 << std::endl;

  // --- WRITE OUTPUT ----------------------------------------------------------
  TFile* out_file = TFile::Open("fm_ana_output.root", "RECREATE");
  out_file->cd();
  h_LL->Write();
  h_LL_fake->Write();
  h_dx->Write();
  h2_xreco_xtrue->Write();
  h2_Ereco_Etrue->Write();
  h2_Ydiff_Zdiff->Write();
  h2_exp_reco->Write();
  likelihood_computer.g_logms->Write();
  out_file->Close();
  return;
}
