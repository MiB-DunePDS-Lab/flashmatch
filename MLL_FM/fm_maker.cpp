
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "Utils.hpp"
#include <flash_matcher.hpp>

struct Row{
  // Event
  int event_id, cluster_id, flash_id, group_id, label;
  // MLL
  float nll, nll_weighted;
  float reco_term_mean, reco_term_std, reco_term_max, reco_term_min;
  float noreco_term_mean, noreco_term_std, noreco_term_max, noreco_term_min;
  float exp_ph_sum, nhit_expected;
  float time_diff;
  // TPC
  float charge, max_charge, time_tpc, y_reco, z_reco;
  int charge_nhits;
  // PDS
  float total_pe, max_pe, time_pds;
  int nhits;
  float flash_reco_y, flash_reco_z;

  // Adding ..............
  float exp_reco_ratio, totalpe_nhits_ratio, exp_close_totalpe_ratio;
  float nhits_expnhits_ratio, close_nhits_exp_ratio;
  float dt_nearest_flash = 1e9;
  int n_close_flashes = 0;
  float close_totalpe = 0;
  int close_nhits = 0;


  // Matching
  float e_reco, flash_cluster_dist;
  // True
  float x_true, y_true, z_true, e_true;
  float purity;
};

float time_window = 20;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void fm_maker(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLConfigs mll_conf = load_ana_config("./configs/ana_config.json");
  std::string ana_file_name      = mll_conf.ana_file_name;
  std::string visibility_dir     = mll_conf.visibility_dir;
  std::string sample_config_file = mll_conf.sample_config_file;
  float light_yield              = mll_conf.light_yield;
  float arapuca_pde              = mll_conf.arapuca_pde;
  float q_cut_high               = mll_conf.q_cut_high;
  float LY_times_PDE             = light_yield * arapuca_pde;
  
  SampleConfigs sample_conf = load_sample_config("./configs/"+sample_config_file);
  std::string input_dir       = sample_conf.input_dir;
  std::string geom_identifier = sample_conf.geom_identifier;
  double trend_thr            = sample_conf.trend_thr;
  float q_cut_low             = sample_conf.q_cut_low;
  
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

  TFile* parametrizer_file  = TFile::Open((input_dir+"MLL_Parametrizer_"+geom_identifier+".root").c_str(), "READ");
  TF1* f_RecoExpDistr       = static_cast<TF1*>(parametrizer_file->Get("f_RecoExpDistr"));
  TF1* f_par1_trend         = static_cast<TF1*>(parametrizer_file->Get("f_par1_trend"));
  TF1* f_par2_trend         = static_cast<TF1*>(parametrizer_file->Get("f_par2_trend"));
  TH2D* h2_exp_reco         = static_cast<TH2D*>(parametrizer_file->Get("h2_exp_reco"));
  TGraphErrors* g_par1      = static_cast<TGraphErrors*>(parametrizer_file->Get("g_par1")); 
  TGraphErrors* g_par2      = static_cast<TGraphErrors*>(parametrizer_file->Get("g_par2"));

  TEfficiency* he_hit_prob = static_cast<TEfficiency*>(parametrizer_file->Get("he_hit_prob"));

  ana_file_name = input_dir+ana_file_name;
  TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
  TTree* tree = static_cast<TTree*>(ana_file->Get("solarnuana/SolarNuAnaTree"));
  std::vector<size_t> MaxChargeIndxs = take_max_charge_indices(tree, "Event", "Charge");
  printf("Total entries in tree: %lld\n", tree->GetEntries());
  printf("Entries passing max charge selection: %lu\n", MaxChargeIndxs.size());

  TTreeReader treeReader(tree);
  // Cluster variables --------------------------------------------------------
  TTreeReaderValue<float> Charge(treeReader, "Charge");
  TTreeReaderValue<float> MaxCharge(treeReader, "MaxCharge");
  TTreeReaderValue<int>   ChargeNHits(treeReader, "NHits");
  TTreeReaderValue<float> Time(treeReader, "Time");
  TTreeReaderValue<float> RecoY(treeReader, "RecoY");
  TTreeReaderValue<float> RecoZ(treeReader, "RecoZ");
  // Flash variables ----------------------------------------------------------
  TTreeReaderArray<float> AdjOpFlashTime(treeReader, "AdjOpFlashTime");
  TTreeReaderArray<float> AdjOpFlashPE(treeReader, "AdjOpFlashPE");
  TTreeReaderArray<float> AdjOpFlashMaxPE(treeReader, "AdjOpFlashMaxPE");
  TTreeReaderArray<int>   AdjOpFlashNHits(treeReader, "AdjOpFlashNHits");
  TTreeReaderArray<float> AdjOpFlashY(treeReader, "AdjOpFlashPur");
  TTreeReaderArray<float> AdjOpFlashRecoY(treeReader, "AdjOpFlashRecoY");
  TTreeReaderArray<float> AdjOpFlashRecoZ(treeReader, "AdjOpFlashRecoZ");
  TTreeReaderArray<float> AdjOpFlashR(treeReader, "AdjOpFlashR");
  TTreeReaderArray<float> AdjOpFlashPur(treeReader, "AdjOpFlashPur");
  TTreeReaderArray<float> AdjOpFlashPEperOpDet(treeReader, "AdjOpFlashPEperOpDet");
  TTreeReaderValue<float> MatchedOpFlashTime(treeReader, "MatchedOpFlashTime");
  TTreeReaderValue<float> MatchedOpFlashPE(treeReader, "MatchedOpFlashPE");
  TTreeReaderValue<float> MatchedOpFlashMaxPE(treeReader, "MatchedOpFlashMaxPE");
  TTreeReaderValue<int>   MatchedOpFlashNHits(treeReader, "MatchedOpFlashNHits");
  TTreeReaderValue<float> MatchedOpFlashPur(treeReader, "MatchedOpFlashPur");
  TTreeReaderValue<float> MatchedOpFlashRecoY(treeReader, "MatchedOpFlashRecoY");
  TTreeReaderValue<float> MatchedOpFlashRecoZ(treeReader, "MatchedOpFlashRecoZ");
  TTreeReaderValue<float> MatchedOpFlashR(treeReader, "MatchedOpFlashR");
  TTreeReaderValue<bool>  MatchedOpFlashCorrectly(treeReader, "MatchedOpFlashCorrectly");
  TTreeReaderArray<float> MatchedOpFlashPEperOpDet(treeReader, "MatchedOpFlashPEperOpDet");
  // True variables -----------------------------------------------------------
  TTreeReaderValue<float> x_true(treeReader, "SignalParticleX");
  TTreeReaderValue<float> y_true(treeReader, "SignalParticleY");
  TTreeReaderValue<float> z_true(treeReader, "SignalParticleZ");
  TTreeReaderValue<float> e_true(treeReader, "SignalParticleE");

  
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

  // --- PREPARE OUTPUT -------------------------------------------------------
  TFile* out_file = TFile::Open((input_dir+"MLL_Features_"+geom_identifier+".root").c_str(), "RECREATE");
  TTree* feature_tree = new TTree("feature_tree", "feature_tree");

  Row r;
  feature_tree->Branch("event_id", &r.event_id, "event_id/I");
  // feature_tree->Branch("cluster_id", &r.cluster_id, "cluster_id/I");
  feature_tree->Branch("flash_id", &r.flash_id, "flash_id/I");
  // feature_tree->Branch("group_id", &r.group_id, "group_id/I");
  feature_tree->Branch("my_label", &r.label, "my_label/I");
  feature_tree->Branch("nll", &r.nll, "nll/F");
  feature_tree->Branch("nll_weighted", &r.nll_weighted, "nll_weighted/F");
  feature_tree->Branch("charge", &r.charge, "charge/F");
  feature_tree->Branch("max_charge", &r.max_charge, "max_charge/F");
  feature_tree->Branch("charge_nhits", &r.charge_nhits, "charge_nhits/F");
  feature_tree->Branch("y_reco", &r.y_reco, "y_reco/F");
  feature_tree->Branch("z_reco", &r.z_reco, "z_reco/F");
  feature_tree->Branch("time_tpc", &r.time_tpc, "time_tpc/F");
  feature_tree->Branch("total_pe", &r.total_pe, "total_pe/F");
  feature_tree->Branch("max_pe", &r.max_pe, "max_pe/F");
  feature_tree->Branch("nhits", &r.nhits, "nhits/I");
  feature_tree->Branch("flash_reco_y", &r.flash_reco_y, "flash_reco_y/F");
  feature_tree->Branch("flash_reco_z", &r.flash_reco_z, "flash_reco_z/F");
  feature_tree->Branch("exp_reco_ratio", &r.exp_reco_ratio, "exp_reco_ratio/F");
  feature_tree->Branch("totalpe_nhits_ratio", &r.totalpe_nhits_ratio, "totalpe_nhits_ratio/F");
  feature_tree->Branch("exp_close_totalpe_ratio", &r.exp_close_totalpe_ratio, "exp_close_totalpe_ratio/F");
  feature_tree->Branch("nhits_expnhits_ratio", &r.nhits_expnhits_ratio, "nhits_expnhits_ratio/F");
  feature_tree->Branch("dt_nearest_flash", &r.dt_nearest_flash, "dt_nearest_flash/F");
  feature_tree->Branch("n_close_flashes", &r.n_close_flashes, "n_close_flashes/I");
  feature_tree->Branch("close_totalpe", &r.close_totalpe, "close_totalpe/F");
  feature_tree->Branch("close_nhits", &r.close_nhits, "close_nhits/I");
  feature_tree->Branch("close_nhits_exp_ratio", &r.close_nhits_exp_ratio, "close_nhits_exp_ratio/F");
  feature_tree->Branch("time_pds", &r.time_pds, "time_pds/F");
  feature_tree->Branch("time_diff", &r.time_diff, "time_diff/F");
  feature_tree->Branch("flash_cluster_dist", &r.flash_cluster_dist, "flash_cluster_dist/F");
  feature_tree->Branch("e_reco", &r.e_reco, "e_reco/F");
  feature_tree->Branch("x_true", &r.x_true, "x_true/F");
  feature_tree->Branch("y_true", &r.y_true, "y_true/F");
  feature_tree->Branch("z_true", &r.z_true, "z_true/F");
  feature_tree->Branch("e_true", &r.e_true, "e_true/F");
  feature_tree->Branch("purity", &r.purity, "purity/F");
  feature_tree->Branch("reco_term_mean", &r.reco_term_mean, "reco_term_mean/F");
  feature_tree->Branch("reco_term_std", &r.reco_term_std, "reco_term_std/F");
  feature_tree->Branch("reco_term_max", &r.reco_term_max, "reco_term_max/F");
  feature_tree->Branch("reco_term_min", &r.reco_term_min, "reco_term_min/F");
  feature_tree->Branch("noreco_term_mean", &r.noreco_term_mean, "noreco_term_mean/F");
  feature_tree->Branch("noreco_term_std", &r.noreco_term_std, "noreco_term_std/F");
  feature_tree->Branch("noreco_term_max", &r.noreco_term_max, "noreco_term_max/F");
  feature_tree->Branch("noreco_term_min", &r.noreco_term_min, "noreco_term_min/F");
  feature_tree->Branch("exp_ph_sum", &r.exp_ph_sum, "exp_ph_sum/F");
  feature_tree->Branch("nhit_expected", &r.nhit_expected, "nhit_expected/F");


  r.event_id = 0;
  std::vector<float> dummy_vec, dummy_vec2;
  for (auto& idx_entry : MaxChargeIndxs){
    treeReader.SetEntry(idx_entry);
    ClusterTPC cluster = ClusterTPC(*Charge, *Time, *RecoY, *RecoZ);
    r.charge = *Charge;
    r.max_charge = *MaxCharge;
    r.charge_nhits = *ChargeNHits;
    r.time_tpc = *Time;
    r.y_reco = *RecoY;
    r.z_reco = *RecoZ;
    r.x_true = *x_true;
    r.y_true = *y_true;
    r.z_true = *z_true;
    r.e_true = *e_true;

    if (*MatchedOpFlashCorrectly){
      std::vector<float> pe_per_opdet(MatchedOpFlashPEperOpDet.begin(), MatchedOpFlashPEperOpDet.end());
      ClusterPDS flash   = ClusterPDS(*MatchedOpFlashTime, pe_per_opdet);
      r.flash_id = -1;
      r.nll = likelihood_computer.GetLikelihoodMatch(cluster, flash, dummy_vec, dummy_vec2, 1.);
      r.nhits = *MatchedOpFlashNHits;
      r.purity = *MatchedOpFlashPur;
      r.flash_reco_y = *MatchedOpFlashRecoY;
      r.flash_reco_z = *MatchedOpFlashRecoZ;
      r.flash_cluster_dist = *MatchedOpFlashR;
      r.nll_weighted = r.nll / (r.nhits * r.nhits);
      r.total_pe = *MatchedOpFlashPE;
      r.max_pe = *MatchedOpFlashMaxPE;
      r.time_pds = *MatchedOpFlashTime;
      r.e_reco = likelihood_computer.E_reco;
      r.label = 2;
      r.reco_term_mean = likelihood_computer.reco_term_mean;
      r.reco_term_std = likelihood_computer.reco_term_std;
      r.reco_term_max = likelihood_computer.reco_term_max;
      r.reco_term_min = likelihood_computer.reco_term_min;
      r.noreco_term_mean = likelihood_computer.noreco_term_mean;
      r.noreco_term_std = likelihood_computer.noreco_term_std;
      r.noreco_term_max = likelihood_computer.noreco_term_max;
      r.noreco_term_min = likelihood_computer.noreco_term_min;
      r.exp_ph_sum = likelihood_computer.exp_ph_sum;
      r.nhit_expected = likelihood_computer.nhit_expected;

      r.exp_reco_ratio = r.exp_ph_sum / (r.total_pe+1.e-6);
      r.totalpe_nhits_ratio = r.total_pe / (r.charge_nhits+1.e-6);
      r.nhits_expnhits_ratio = r.nhits / (r.nhit_expected+1.e-6);
      r.time_diff = *Time - *MatchedOpFlashTime;
      
      r.dt_nearest_flash = 1e9;
      r.n_close_flashes = 0;
      r.close_totalpe = 0;
      r.close_nhits = 0;
      for (size_t ii = 0; ii < AdjOpFlashTime.GetSize(); ii++){
        float dt = std::abs(*MatchedOpFlashTime - AdjOpFlashTime.At(ii));
        if (dt < r.dt_nearest_flash) r.dt_nearest_flash = dt;
        if (dt < time_window){
          r.n_close_flashes++;
          r.close_totalpe += AdjOpFlashPE.At(ii);
          r.close_nhits += AdjOpFlashNHits.At(ii);
        }
      }
      r.close_nhits_exp_ratio = r.close_nhits / (r.nhit_expected+1.e-6);
      
      feature_tree->Fill();
    }
    
    for (size_t idx_flash = 0; idx_flash < AdjOpFlashTime.GetSize(); idx_flash++){
      if (AdjOpFlashPE.At(idx_flash) == *MatchedOpFlashPE && AdjOpFlashTime.At(idx_flash) == *MatchedOpFlashTime){
        continue; // Skip the correctly matched flash since it's already processed
      }

      std::vector<float> pe_per_opdet(AdjOpFlashPEperOpDet.begin() + idx_flash*geom.n_opdet, AdjOpFlashPEperOpDet.begin() + (idx_flash+1)*geom.n_opdet);
      ClusterPDS flash   = ClusterPDS(AdjOpFlashTime.At(idx_flash), pe_per_opdet);
      r.flash_id = idx_flash;
      r.nll = likelihood_computer.GetLikelihoodMatch(cluster, flash, dummy_vec, dummy_vec2, 1.);
      r.nhits = AdjOpFlashNHits[idx_flash];
      r.purity = AdjOpFlashY[idx_flash];
      r.flash_reco_y = AdjOpFlashRecoY[idx_flash];
      r.flash_reco_z = AdjOpFlashRecoZ[idx_flash];
      r.flash_cluster_dist = AdjOpFlashR[idx_flash];
      r.nll_weighted = r.nll / (r.nhits * r.nhits);
      r.total_pe = AdjOpFlashPE[idx_flash];
      r.max_pe = AdjOpFlashMaxPE[idx_flash];
      r.time_pds = AdjOpFlashTime[idx_flash];
      r.e_reco = likelihood_computer.E_reco;
      r.label = AdjOpFlashPur[idx_flash] > 0. ? 1 : 0;
      r.reco_term_mean = likelihood_computer.reco_term_mean;
      r.reco_term_std = likelihood_computer.reco_term_std;
      r.reco_term_max = likelihood_computer.reco_term_max;
      r.reco_term_min = likelihood_computer.reco_term_min;
      r.noreco_term_mean = likelihood_computer.noreco_term_mean;
      r.noreco_term_std = likelihood_computer.noreco_term_std;
      r.noreco_term_max = likelihood_computer.noreco_term_max;
      r.noreco_term_min = likelihood_computer.noreco_term_min;
      r.exp_ph_sum = likelihood_computer.exp_ph_sum;
      r.nhit_expected = likelihood_computer.nhit_expected;

      r.exp_reco_ratio = r.exp_ph_sum / (r.total_pe+1.e-6);
      r.totalpe_nhits_ratio = r.total_pe / (r.charge_nhits+1.e-6);
      r.nhits_expnhits_ratio = r.nhits / (r.nhit_expected+1.e-6);
      r.time_diff = *Time - AdjOpFlashTime.At(idx_flash);
      
      r.dt_nearest_flash = 1e9;
      r.n_close_flashes = 0;
      r.close_totalpe = 0;
      r.close_nhits = 0;
      for (size_t ii = 0; ii < AdjOpFlashTime.GetSize(); ii++){
        float dt = std::abs(*MatchedOpFlashTime - AdjOpFlashTime.At(ii));
        if (dt < r.dt_nearest_flash) r.dt_nearest_flash = dt;
        if (dt < time_window){
          r.n_close_flashes++;
          r.close_totalpe += AdjOpFlashPE.At(ii);
          r.close_nhits += AdjOpFlashNHits.At(ii);
        }
      }
      r.close_nhits_exp_ratio = r.close_nhits / (r.nhit_expected+1.e-6);

      feature_tree->Fill();
    } // loop over flashes
    r.event_id++;
  } // loop over SolarNuAnaTree entries
 
  out_file->cd();
  feature_tree->Write("", TObject::kOverwrite);
  out_file->Close();
  ana_file->Close();
  return;
}
