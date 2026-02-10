#ifndef UTILS_HPP
#define UTILS_HPP

#include "TGraphErrors.h"
#ifndef NOPDET
#define NOPDET 480
#endif // !NOPDET

// GENERAL
#include <cstddef>
#include <vector>
#include <iostream>
#include <numeric>
#include <fstream>
#include <nlohmann/json.hpp>

// ROOT
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include "TEfficiency.h"
#include "TH2.h"
#include "TTreeReaderArray.h"

using json = nlohmann::json;
namespace fs = std::filesystem;


// --- PURE UTILITIES ---------------------------------------------------------
#include <vector>
#include <string>
#include <filesystem>

inline std::vector<std::string> get_list_of_files_in_folder(std::string folder_name, std::string ends_with){
  std::vector<std::string> result;

  if (!fs::exists(folder_name) || !fs::is_directory(folder_name))
    return result;

  for (const auto& entry : fs::directory_iterator(folder_name))
  {
    if (!entry.is_regular_file())
      continue;

    std::string filename = entry.path().filename().string();

    if (filename.size() >= ends_with.size() &&
      filename.compare(filename.size() - ends_with.size(),
                       ends_with.size(),
                       ends_with) == 0)
    {
      result.push_back(filename);
    }
  }

  return result;
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


// --- FLASH MATCHER UTILS ----------------------------------------------------

class Event {
  public:
    static constexpr size_t n_opdet = NOPDET;
    double true_energy;
    double vertex_coor[3];
    double mu[n_opdet];

    Event& operator=(const Event& other) {
      if (this != &other) {
        true_energy = other.true_energy;
        std::copy(std::begin(other.vertex_coor), std::end(other.vertex_coor), vertex_coor);
        std::copy(std::begin(other.mu), std::end(other.mu), mu);
      }
      return *this;
    }

    bool operator==(const Event& other) {
      return (true_energy == other.true_energy &&
              std::equal(std::begin(vertex_coor), std::end(vertex_coor), std::begin(other.vertex_coor)) &&
              std::equal(std::begin(mu), std::end(mu), std::begin(other.mu)));
    }

  Event() : true_energy(0), vertex_coor{0, 0, 0}, mu{0} {}
  
  Event(const Double_t true_energy,
        const TTreeReaderArray<Double_t>& vertex_coor,
        const TTreeReaderArray<Double_t>& mu) {
    this->true_energy = true_energy;
    for (size_t i = 0; i < 3; ++i) {
      this->vertex_coor[i] = vertex_coor[i];
    }
    for (size_t i = 0; i < n_opdet; ++i) {
      this->mu[i] = mu[i];
    }
  }
};

class Flash {
  public:
    static constexpr size_t n_opdet = NOPDET;
    double reco_pe[n_opdet];
    double hit_t[n_opdet];
    bool hit[n_opdet];

    Flash& operator=(const Flash& other) {
      if (this != &other) {
        std::copy(std::begin(other.reco_pe), std::end(other.reco_pe), reco_pe);
        std::copy(std::begin(other.hit_t), std::end(other.hit_t), hit_t);
        std::copy(std::begin(other.hit), std::end(other.hit), hit);
      }
      return *this;
    }

    bool operator==(const Flash& other) {
      for (size_t i = 0; i < n_opdet; ++i) {
        if (reco_pe[i] != other.reco_pe[i] || hit_t[i] != other.hit_t[i] || hit[i] != other.hit[i]) {
          return false;
        }
      }
      return true;
    }

    Flash() : reco_pe{0}, hit_t{0}, hit{false} {}
    Flash(const TTreeReaderArray<Double_t>& reco_pe,
          const TTreeReaderArray<Double_t>& hit_t,
          const TTreeReaderArray<Bool_t>& hit) {
      for (size_t i = 0; i < n_opdet; ++i) {
        this->reco_pe[i] = reco_pe[i];
        this->hit_t[i] = hit_t[i];
        this->hit[i] = hit[i];
      }
    }
};

inline std::vector<double> compute_loglikelihood_terms(Event& event,
                                                       Flash& flash,
                                                       TEfficiency* he_ProbHit_ExpPe,
                                                       TH2D* h2_exp_reco) {
  std::vector<double> terms(event.n_opdet);
  double term = 0.;
  for(size_t idx_opdet=0; idx_opdet<event.n_opdet; idx_opdet++){
    double P_hit_mu = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(event.mu[idx_opdet]));
    if(flash.hit[idx_opdet]){
      term = P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(flash.reco_pe[idx_opdet], event.mu[idx_opdet]));
    } else {
      term = 1. - P_hit_mu;
    }
    terms[idx_opdet] = -log(term);
  }
  return terms;
}

inline std::vector<double> compute_loglikelihood_terms(Event& event,
                                                       Flash& flash,
                                                       TF1*  f_ProbHit_ExpPe,
                                                       TH2D* h2_exp_reco) {
  std::vector<double> terms(event.n_opdet);
  double term = 0.;
  for(size_t idx_opdet=0; idx_opdet<event.n_opdet; idx_opdet++){
    double P_hit_mu = f_ProbHit_ExpPe->Eval(event.mu[idx_opdet]);
    if(flash.hit[idx_opdet]){
      term = P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(flash.reco_pe[idx_opdet], event.mu[idx_opdet]));
    } else {
      term = 1. - P_hit_mu;
    }
    terms[idx_opdet] = -log(term);
  }
  return terms;
}

template<typename T>
bool flash_matcher(Event& event_a, Flash& flash_a, Event& event_b, Flash& flash_b, T he_ProbHit_ExpPe, TH2D* h2_exp_reco) {
  std::vector<double> ea_fa_vector = compute_loglikelihood_terms(event_a, flash_a, he_ProbHit_ExpPe, h2_exp_reco); 
  double ea_fa = std::accumulate(ea_fa_vector.begin(), ea_fa_vector.end(), 0.0);
  std::vector<double> eb_fb_vector = compute_loglikelihood_terms(event_b, flash_b, he_ProbHit_ExpPe, h2_exp_reco);
  double eb_fb = std::accumulate(eb_fb_vector.begin(), eb_fb_vector.end(), 0.0);
  std::vector<double> ea_fb_vector = compute_loglikelihood_terms(event_a, flash_b, he_ProbHit_ExpPe, h2_exp_reco);
  double ea_fb = std::accumulate(ea_fb_vector.begin(), ea_fb_vector.end(), 0.0);
  std::vector<double> eb_fa_vector = compute_loglikelihood_terms(event_b, flash_a, he_ProbHit_ExpPe, h2_exp_reco);
  double eb_fa = std::accumulate(eb_fa_vector.begin(), eb_fa_vector.end(), 0.0);
 
  std::vector<double> es_fs = {ea_fa, eb_fb, ea_fb, eb_fa};
  double inf = std::numeric_limits<double>::infinity();

  if ( std::count_if(es_fs.begin(), es_fs.end(), [&](double x) { return x == inf; }) == 1 ) {
    if (ea_fb == inf || eb_fa == inf) {
      return true;
    } else {
      return false;
    }
  }
  else {
    double min = std::min({ea_fa, eb_fb, ea_fb, eb_fa});
    if (std::count_if(es_fs.begin(), es_fs.end(), [&](double x) { return x == min; }) > 1) {
      std::cout << "they are equal to " << min << std::endl;
    }
    if (ea_fa == min || eb_fb == min) {
      return true;
    } else {
      return false;
    }
  }
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

// Sigmoid times positive_erf
inline double sigmoid_erf(double* x, double* par){
  return sigmoid(x, &par[0]) * positive_erf(x, &par[2]);
}

// Sigmoid times sigmoid times positive_erf 
inline double sigmoid_sigmoid_erf(double* x, double* par){
  return sigmoid(x, &par[0]) * sigmoid(x, &par[2]) * positive_erf(x, &par[4]);
}


struct MLLcconfigs{
  std::string input_dir;
  std::string visibility_file_name;
  int max_nfiles;
  double fit_Qcorr_Etrue_low;
  double fit_Qcorr_Etrue_up;
  float pe_low;
  float pe_up;
  float light_yield;
  float arapuca_pde;
  float min_visibility;
  float yz_fiducial_cut;
  float x_fiducial_cut;
  // bool verbose;
  double fit_trend_low;
  double fit_trend_up;
  bool rebin_h2; // Rebin h2_exp_reco_norm
  int rebinning; // Rebinning factor for h2_exp_reco_norm
  double trend_thr;
};

inline MLLcconfigs load_ana_config(const std::string &filename){
  std::ifstream file(filename);
  if (!file) {
    throw std::runtime_error("Could not open config file: " + filename);
  }

  json j;
  file >> j;
  MLLcconfigs config;
  config.input_dir            = j.at("input_dir").get<std::string>();
  config.visibility_file_name = j.at("visibility_file_name").get<std::string>();
  config.fit_Qcorr_Etrue_low  = j.at("fit_Qcorr_Etrue_low").get<double>();
  config.fit_Qcorr_Etrue_up   = j.at("fit_Qcorr_Etrue_up").get<double>();
  config.max_nfiles           = j.at("max_nfiles").get<int>();
  config.pe_low               = j.at("pe_low").get<float>();
  config.pe_up                = j.at("pe_up").get<float>();
  config.light_yield          = j.at("light_yield").get<float>();
  config.arapuca_pde          = j.at("arapuca_pde").get<float>();
  config.min_visibility       = j.at("min_visibility").get<float>();
  config.yz_fiducial_cut      = j.at("yz_fiducial_cut").get<float>();
  config.x_fiducial_cut       = j.at("x_fiducial_cut").get<float>();
  // config.verbose              = j.at("verbose").get<bool>();
  config.fit_trend_low        = j.at("fit_trend_low").get<double>();
  config.fit_trend_up         = j.at("fit_trend_up").get<double>();
  config.rebin_h2             = j.at("rebin_h2").get<bool>();
  config.rebinning            = j.at("rebinning").get<int>();
  config.trend_thr            = j.at("trend_thr").get<double>();

  return config;
}

#endif // UTILS_HPP
