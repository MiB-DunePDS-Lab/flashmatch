#ifndef FLASH_MATCHER_HPP
#define FLASH_MATCHER_HPP

#include <TFile.h>

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
  TGraphErrors* g_logms;
  TF1* f_reco_prob;

  float GetLikelihoodMatch(const ClusterTPC& tpc_cluster,
                           const ClusterPDS& pds_cluster){
    
    float delta_time = tpc_cluster.time_tpc - pds_cluster.time_pds ;
    E_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, delta_time, tpc_cluster.charge);
    float reco_pe, exp_ph;
    x_reco = delta_time*drift_velocity;
    float vertex_coor[3] = {x_reco, tpc_cluster.y_reco, tpc_cluster.z_reco};
    vertex_coor[0] = std::max(vertex_coor[0], tpc_min[0]); vertex_coor[0] = std::min(vertex_coor[0], tpc_max[0]-3);
    vertex_coor[1] = std::max(vertex_coor[1], tpc_min[1]); vertex_coor[1] = std::min(vertex_coor[1], tpc_max[1]);
    vertex_coor[2] = std::max(vertex_coor[2], tpc_min[2]); vertex_coor[2] = std::min(vertex_coor[2], tpc_max[2]);
    int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);

    float x_min_g_logms = g_logms->GetX()[0];

    float log_likelihood = 0.;

    if (n_opdet != pds_cluster.reco_pes->size()){
      std::cerr << "Error: Number of OPDet does not match the size of reco_pes vector!" << std::endl;
      exit(1);
    }

    for (size_t idx_opdet=0; idx_opdet<pds_cluster.reco_pes->size(); idx_opdet++){
      float voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];

      exp_ph = E_reco*LY_times_PDE*voxel_vis;
      if(exp_ph==0) exp_ph = E_reco*LY_times_PDE*1.e-15;
      float P_hit_mu = f_reco_prob->Eval(exp_ph);

      reco_pe = pds_cluster.reco_pes->at(idx_opdet);

      if (reco_pe>0.0){
        if (reco_pe > trend_thr){
          f_lognormal->SetParameters(f_logms_trend->Eval(exp_ph), f_sigmas_trend->Eval(exp_ph));
          log_likelihood += log(P_hit_mu*f_lognormal->Eval(reco_pe));
        } else {
          f_lognormal->SetParameters(g_logms->Eval(exp_ph), g_sigmas->Eval(exp_ph));
          log_likelihood += log(P_hit_mu*f_lognormal->Eval(reco_pe));
        }
        // log_likelihood += log(P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(reco_pe, exp_ph)));
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

    return;
  } // 

}; // LikelihoodComputer

#endif // FLASH_MATCHER_HPP
