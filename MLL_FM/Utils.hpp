#ifndef UTILS_HPP
#define UTILS_HPP

// GENERAL
#include <cstddef>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <nlohmann/json.hpp>

// ROOT
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TMath.h>

using json = nlohmann::json;


// --- PURE UTILITIES ---------------------------------------------------------
template <typename T>
std::vector<T> treeBranch_to_vector(TTree* tree, const std::string& branch_name) {
  std::vector<T> vec;
  tree->Draw(branch_name.c_str(), "", "goff");
  int n = tree->GetSelectedRows();
  double* array = tree->GetV1();

  vec.reserve(n);
  for (int i = 0; i < n; ++i) {
    vec.push_back(static_cast<T>(array[i]));
  }

  return vec;
}

// function to take the indices of the cluster with the max charge
inline std::vector<size_t> take_max_charge_indices(TTree* tree,
                                            const std::string& event_branch_name,
                                            const std::string& charge_branch_name) {
  std::vector<int> Events = treeBranch_to_vector<int>(tree, event_branch_name);
  std::vector<float> Charges = treeBranch_to_vector<float>(tree, charge_branch_name);
  std::vector<size_t> MaxChargeIndxs;
  int old_event = -1;
  for (size_t i=0; i<Events.size(); i++) {
    if (old_event==Events[i]) continue;
    old_event = Events[i];

    size_t max_charge_index = i;
    size_t j=i+1;
    while(j<Events.size() && Events[j] == Events[i]){
      if (Charges[j]>Charges[i]) max_charge_index = j;
      j++;
    }
    MaxChargeIndxs.push_back(max_charge_index);
  }
  return MaxChargeIndxs;
}


inline TGraphErrors* th2d_to_tgraph_mpv(TH2D* h2, const std::string& name){
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name.c_str()); g->SetTitle(h2->GetTitle());
  for(int idx_x=0; idx_x<=h2->GetNbinsX()+1; idx_x++){
    TH1D* h1_proj = h2->ProjectionY("h1_proj", idx_x, idx_x);
    if(h1_proj->GetEntries() > 0){
      double mpv = h1_proj->GetXaxis()->GetBinCenter(h1_proj->GetMaximumBin());
      double err = h1_proj->GetBinWidth(0); // CHECK THIS
      g->SetPoint(idx_x, h2->GetXaxis()->GetBinCenter(idx_x), mpv);
      g->SetPointError(idx_x, 0., err);
    }
  }
  return g; 
}

template <typename Tuple, std::size_t... Is>
double get_by_index(const Tuple& tuple, std::size_t index, std::index_sequence<Is...>){
  double result = 0.;
  ((Is == index ? (result = static_cast<double>(std::get<Is>(tuple))) : 0.), ...);
  return result;
}

template <typename... T>
double get_maximum_from_tuples(const std::vector<std::tuple<T...>>& tuples,
                               size_t column){
  auto seq = std::index_sequence_for<T...>{};

  auto max_it = std::max_element(tuples.begin(), tuples.end(),
                                 [column, seq](const auto& a, const auto& b){
                                   return get_by_index(a, column, seq) < get_by_index(b, column, seq);
                                 });

  return get_by_index(*max_it, column, seq);
}

template <typename... T>
double get_minimum_from_tuples(const std::vector<std::tuple<T...>>& tuples,
                               size_t column){
  auto seq = std::index_sequence_for<T...>{};

  auto min_it = std::min_element(tuples.begin(), tuples.end(),
                                 [column, seq](const auto& a, const auto& b){
                                   return get_by_index(a, column, seq) < get_by_index(b, column, seq);
                                 });

  return get_by_index(*min_it, column, seq);
}


inline float give_me_Ereco(float calib_c, float calib_slope, float corr_lambda,
                    float dt, float charge){
  
  float q_corr = charge * exp(dt * corr_lambda);
  float E_reco = (q_corr - calib_c) / calib_slope;
  
  return E_reco;
}


// --- HANDLE VISIBILITIES ----------------------------------------------------

// From the hgridX histograms created by the PhotonVisibilityExpoort_module.cc, it returns
// a vector where the z-th + y-th * hgrid[2]->GetNbinsX() + x-th * hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()
// entry value is the number of the photoVisMap entry corresponding to the x,y,z coordinates.
template <typename T>
std::vector<int> GetCryoToTPCMap(TH1D* hgrid[3], T tpc_min[3], T tpc_max[3]){
  size_t n_total_entries = hgrid[0]->GetNbinsX() * hgrid[1]->GetNbinsX() * hgrid[2]->GetNbinsX();
  std::vector<int> total_to_tpc(n_total_entries, -1);
  int tpc_entry = 0;
  for (int ix=0; ix<hgrid[0]->GetNbinsX(); ix++){
    for (int iy=0; iy<hgrid[1]->GetNbinsX(); iy++){
      for (int iz=0; iz<hgrid[2]->GetNbinsX(); iz++){
        double x = hgrid[0]->GetBinCenter(ix);
        double y = hgrid[1]->GetBinCenter(iy);
        double z = hgrid[2]->GetBinCenter(iz);
        
        if (x > tpc_min[0] && x < tpc_max[0] && 
            y > tpc_min[1] && y < tpc_max[1] && 
            z > tpc_min[2] && z < tpc_max[2]) {
          total_to_tpc[iz + iy*hgrid[2]->GetNbinsX() + ix*hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()] = tpc_entry;
          tpc_entry++;
        }
        
      }
    }
  }
  return total_to_tpc;
}

template <typename T>
bool isInFiducialVolume(T vertex_coor[3], T vol_min[3], T vol_max[3], T x_cut){
  return (vertex_coor[0] > vol_min[0] && vertex_coor[0] < vol_max[0] &&
          vertex_coor[1] > vol_min[1] && vertex_coor[1] < vol_max[1] &&
          vertex_coor[2] > vol_min[2] && vertex_coor[2] < vol_max[2] &&
          abs(vertex_coor[0]) > x_cut);
}

template <typename T>
int GetTPCIndex(T vertex_coor[3], TH1D* hgrid[3], std::vector<int>& total_to_tpc){
  int bin_x = hgrid[0]->FindBin(vertex_coor[0]);
  int bin_y = hgrid[1]->FindBin(vertex_coor[1]);
  int bin_z = hgrid[2]->FindBin(vertex_coor[2]);

  int tpc_index = total_to_tpc[bin_z + bin_y*hgrid[2]->GetNbinsX()
                               + bin_x*hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()];

  return tpc_index; // Inside TPC
}


// --- FIT FUNCTIONS ----------------------------------------------------------
// Sigmoid function
inline double sigmoid(double* x, double* par){
  return 1 / (1 + exp(-(x[0]-par[0])/par[1]));
}

// Positive error function (0 to 1)
inline double positive_erf(double* x, double* par){
  return 0.5*(1.+TMath::Erf((x[0]-par[0])/par[1]));
}

// Sigmoid times sigmoid times positive_erf 
inline double sigmoid_sigmoid_erf(double* x, double* par){
  return sigmoid(x, &par[0]) * sigmoid(x, &par[2]) * positive_erf(x, &par[4]);
}

// Gamma distribution
inline double gamma_dist(double* x, double* par){
  if (x[0] <= 0) return 0.;
  double k = par[0];
  double theta = par[1];
  return (pow(x[0], k-1) * exp(-x[0]/theta)) / (pow(theta, k) * TMath::Gamma(k));
}

// Weibull distribution
inline double weibull_dist(double* x, double* par){
  if (x[0] <= 0) return 0.;
  double k = par[0];
  double lambda = par[1];
  return (k / lambda) * pow(x[0] / lambda, k - 1) * exp(-pow(x[0] / lambda, k));
}

// Log-logistic distribution
inline double log_logistic_dist(double* x, double* par){
  if (x[0] <= 0) return 0.;
  double alpha = par[0];
  double beta = par[1];
  return (beta / alpha) * pow(x[0]/alpha, beta-1) / pow(1 + pow(x[0]/alpha, beta), 2);
}

struct MLLConfigs{
  std::string ana_file_name;
  std::string visibility_dir;
  std::string sample_config_file;
  double fit_Qcorr_Etrue_low;
  double fit_Qcorr_Etrue_up;
  float pe_up;
  float light_yield;
  float arapuca_pde;
  float min_visibility;
  float yz_fiducial_cut;
  float x_fiducial_cut;
  // bool verbose;
  float q_cut_low;                
  float q_cut_high;  
};

inline MLLConfigs load_ana_config(const std::string &filename){
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }

  json j;
  file >> j;
  MLLConfigs config;
  config.ana_file_name        = j.at("ana_file_name").get<std::string>();
  config.visibility_dir       = j.at("visibility_dir").get<std::string>();
  config.sample_config_file   = j.at("sample_config_file").get<std::string>();
  config.fit_Qcorr_Etrue_low  = j.at("fit_Qcorr_Etrue_low").get<double>();
  config.fit_Qcorr_Etrue_up   = j.at("fit_Qcorr_Etrue_up").get<double>();
  config.pe_up                = j.at("pe_up").get<float>();
  config.light_yield          = j.at("light_yield").get<float>();
  config.arapuca_pde          = j.at("arapuca_pde").get<float>();
  config.min_visibility       = j.at("min_visibility").get<float>();
  config.yz_fiducial_cut      = j.at("yz_fiducial_cut").get<float>();
  config.x_fiducial_cut       = j.at("x_fiducial_cut").get<float>();
  // config.verbose              = j.at("verbose").get<bool>();
  config.q_cut_low            = j.at("q_cut_low").get<float>();
  config.q_cut_high           = j.at("q_cut_high").get<float>();
  return config;
}

struct SampleConfigs{
  std::string input_dir;
  std::string geom_identifier;
  std::string distribution;
  double fit_trend_low;
  double fit_trend_up;
  double trend_thr;
  double q_cut_low;
  bool apply_cut;                 
};

inline SampleConfigs load_sample_config(const std::string &filename){
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }

  json j;
  file >> j;
  SampleConfigs config;
  config.input_dir         = j.at("input_dir").get<std::string>();
  config.geom_identifier   = j.at("geom_identifier").get<std::string>();
  config.distribution      = j.at("distribution").get<std::string>();
  config.fit_trend_low     = j.at("fit_trend_low").get<double>();
  config.fit_trend_up      = j.at("fit_trend_up").get<double>();
  config.trend_thr         = j.at("trend_thr").get<double>();
  config.q_cut_low         = j.at("q_cut_low").get<double>();
  config.apply_cut         = j.at("apply_cut").get<bool>();
  return config;
}

#endif // UTILS_HPP
