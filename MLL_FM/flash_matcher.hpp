#ifndef FLASH_MATCHER_HPP
#define FLASH_MATCHER_HPP

#include <TFile.h>
#include <cmath>
#include <numeric>

#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include <TMVA/TSpline1.h>
#include "Utils.hpp"

class VertexInfo {
  public:
    float x;
    float y;
    float z;
    float energy;

    VertexInfo& operator=(const VertexInfo& other) {
      if (this != &other) {
        x = other.x;
        y = other.y;
        z = other.z;
        energy = other.energy;
      }
      return *this;
    }

    bool operator==(const VertexInfo& other) {
      return (x == other.x && y == other.y && z == other.z && energy == other.energy);
    }

    VertexInfo() : x(0), y(0), z(0), energy(0) {}

    VertexInfo(const float x, const float y, const float z, const float energy) {
      this->x = x; this->y = y; this->z = z; this->energy = energy;
    }
};

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
    std::vector<float>* reco_pes;

    ClusterPDS& operator=(const ClusterPDS& other) {
      if (this != &other) {
        time_pds = other.time_pds;
        reco_pes = other.reco_pes;
      }
      return *this;
    }

    bool operator==(const ClusterPDS& other) {
      return (time_pds == other.time_pds &&
          reco_pes == other.reco_pes &&
          *reco_pes == *other.reco_pes);
    }

    ClusterPDS() : time_pds(0), reco_pes(nullptr) {
      reco_pes = new std::vector<float>();
    }


    ClusterPDS(const float& time_pds,
        std::vector<float>* reco_pes) {
      this->time_pds = time_pds;
      this->reco_pes = reco_pes;
    }

    bool IsEmpty(){
      return (reco_pes->empty());
    }
};

class LikelihoodComputer{

public:

  float E_reco;
  float x_reco;
  TMVA::TSpline1* g_he = nullptr;
  float xprob_max = 0.;
  
  float GetLikelihoodMatch(const ClusterTPC& tpc_cluster,
                           const ClusterPDS& pds_cluster,
                           std::vector<float>& reco_terms,
                           std::vector<float>& noreco_terms,
                           float x_sign = 1.) {

    if (reco_terms.size() >= 0) reco_terms.clear();
    if (noreco_terms.size() >= 0) noreco_terms.clear();
    
    float delta_time = tpc_cluster.time_tpc - pds_cluster.time_pds ;
    E_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, delta_time, tpc_cluster.charge);
    float reco_pe, exp_ph;
    x_reco = delta_time*drift_velocity+anode_x;
    float vertex_coor[3] = {x_reco, tpc_cluster.y_reco, tpc_cluster.z_reco};
    vertex_coor[0] = std::max(vertex_coor[0], tpc_min[0]+15); vertex_coor[0] = std::min(vertex_coor[0], tpc_max[0]-15);
    vertex_coor[1] = std::max(vertex_coor[1], tpc_min[1]+15); vertex_coor[1] = std::min(vertex_coor[1], tpc_max[1]-15);
    vertex_coor[2] = std::max(vertex_coor[2], tpc_min[2]+15); vertex_coor[2] = std::min(vertex_coor[2], tpc_max[2]-15);
    int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);

    float x_min_g_par1s = g_par1->GetX()[0];

    float log_likelihood = 0.;
    // std::cout << "yyy" << log_likelihood << std::endl;
    float term = 0.;
    float n_hit = std::count_if(pds_cluster.reco_pes->begin(), pds_cluster.reco_pes->end(), [](float pe){ return pe > 0; });

    if (n_opdet != pds_cluster.reco_pes->size()){
      std::cerr << "Error: Number of OPDet does not match the size of reco_pes vector!" << std::endl;
      exit(1);
    }

    // float weight_hit = 1.;
    float weight_hit = 1./(n_hit*n_hit);
    // float weight_unhit = 1.;
    float weight_unhit = 1./(n_hit*n_hit);
    float weight_sum = 1.;
    // float weight_sum = 1./(std::accumulate(pds_cluster.reco_pes->begin(), pds_cluster.reco_pes->end(), 0.));
    for (size_t idx_opdet=0; idx_opdet<pds_cluster.reco_pes->size(); idx_opdet++){
      float voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];

      exp_ph = E_reco*LY_times_PDE*voxel_vis;
      if(exp_ph==0) exp_ph = E_reco*LY_times_PDE*1.e-15;
      float P_hit_mu = (exp_ph<xprob_max) ? g_he->Eval(exp_ph) : 1.; // Avoid weird extrapolation where we
      //                                                             // have no points in the efficiency graph,
      //                                                             // and just assume that the hit probability is 1 for very high expected PE.
      if (P_hit_mu <= 0.) P_hit_mu = 1.e-4;
      if (P_hit_mu >= 1.) P_hit_mu = 1. - 1.e-4;

      reco_pe = pds_cluster.reco_pes->at(idx_opdet);

      if (reco_pe>0.0){
        n_hit++;
        if (reco_pe > trend_thr){
          f_RecoExpDistr->SetParameters(f_par1_trend->Eval(exp_ph), f_par2_trend->Eval(exp_ph));
          // term = log(P_hit_mu*f_RecoExpDistr->Eval(reco_pe))*weight_hit;
          // term = log(P_hit_mu*f_RecoExpDistr->Eval(reco_pe)/f_RecoExpDistr->Eval(exp(f_par1_trend->Eval(exp_ph)-pow(f_par2_trend->Eval(exp_ph),2))))*weight_hit;
          term = log(P_hit_mu*f_RecoExpDistr->Integral(reco_pe-sqrt(reco_pe)*0.5, reco_pe+sqrt(reco_pe)*0.5, 0.001))*weight_hit;
          log_likelihood += term;
          reco_terms.push_back(term);
          // std::cout <<  "t " << term << " " << log_likelihood << std::endl;
        } else {
          f_RecoExpDistr->SetParameters(g_par1->Eval(exp_ph), g_par2->Eval(exp_ph));
          // term = log(P_hit_mu*f_RecoExpDistr->Eval(reco_pe))*weight_hit;
          // term = log(P_hit_mu*f_RecoExpDistr->Eval(reco_pe)/f_RecoExpDistr->Eval(exp(g_par1->Eval(exp_ph)-pow(g_par2->Eval(exp_ph),2))))*weight_hit;
          term = log(P_hit_mu*f_RecoExpDistr->Integral(reco_pe-sqrt(reco_pe)*0.5, reco_pe+sqrt(reco_pe)*0.5, 0.001))*weight_hit;
          log_likelihood += (term);
          reco_terms.push_back(term);
          // std::cout << "d " << term << " " << log_likelihood << std::endl;
        }
      } else {
        term = log(1. - P_hit_mu)*weight_unhit;
        log_likelihood += term;
        noreco_terms.push_back(term);
          // std::cout << "e " << term << " " << log_likelihood << std::endl;
      }
    }

    // std::cout << "xxx" << log_likelihood << std::endl;
    x_reco *= x_sign;
    // std::cout << log_likelihood << std::endl;
    return log_likelihood*weight_sum;
  } // GetLikelihoodMatch

  // LikelihoodComputer constructor
  LikelihoodComputer(TString visibility_file_name,
                     size_t n_opdet,
                     float anode_x,
                     float drift_velocity,
                     float LY_times_PDE,
                     TEfficiency* he_hit_prob,
                     TF1* f_RecoExpDistr,
                     TF1* f_par1_trend,
                     TF1* f_par2_trend,
                     TGraphErrors* g_par1,
                     TGraphErrors* g_par2,
                     double trend_thr,
                     float calib_c, 
                     float calib_slope, 
                     float corr_lambda,
                     TH2D* h2_exp_reco = nullptr)
    : visibility_file_name(visibility_file_name),
    n_opdet(n_opdet),
    anode_x(anode_x),
    drift_velocity(drift_velocity),
    LY_times_PDE(LY_times_PDE),
    he_hit_prob(he_hit_prob),
    f_RecoExpDistr(f_RecoExpDistr),
    f_par1_trend(f_par1_trend),
    f_par2_trend(f_par2_trend),
    g_par1(g_par1),
    g_par2(g_par2),
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
  float anode_x;
  float drift_velocity;
  float LY_times_PDE;
  std::vector<std::vector<float>> opDet_visMapDirect;
  TEfficiency* he_hit_prob;
  TF1* f_RecoExpDistr;
  TF1* f_par1_trend;
  TF1* f_par2_trend;
  double trend_thr;
  TGraphErrors* g_par1;
  TGraphErrors* g_par2;
  float calib_c, calib_slope, corr_lambda;
  TH2D* h2_exp_reco = nullptr; // Optional, can be nullptr

  
  
  // set inside the class
  TH1D* hgrid[3] = {nullptr};
  std::vector<int> cryo_to_tpc;
  float tpc_min[3] = {0., 0., 0.};
  float tpc_max[3] = {0., 0., 0.};

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
    tpc_min[0] = float(coor_dim[0]); tpc_min[1] = float(coor_dim[1]); tpc_min[2] = float(coor_dim[2]); // x,y,z
    tDimensions->GetEntry(1);
    tpc_max[0] = float(coor_dim[0]); tpc_max[1] = float(coor_dim[1]); tpc_max[2] = float(coor_dim[2]); // x,y,z
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

    TGraph* g_he_graph = (TGraph*)he_hit_prob->CreateGraph();
    xprob_max = g_he_graph->GetX()[g_he_graph->GetN()-1];
    g_he = new TMVA::TSpline1("spline", g_he_graph);
    
    return;
  } // 

}; // LikelihoodComputer

#endif // FLASH_MATCHER_HPP
