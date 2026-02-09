#include "RtypesCore.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TEntryList.h"
#include "TTree.h"
#include <cstddef>
#include <filesystem>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>

#include "Utils.hpp"


bool get_next_entry(TTree* tree, int& entry, int& entries){
  if (++entry>entries) return false;
  tree->GetEntry(entry);
  return true;
}

bool one_event_one_entry(TTree* tree, int& iev, int& entry){
  int current_iev = iev;
  int current_ent = entry;
  
  if (current_ent+1 < tree->GetEntries()) {
    tree->GetEntry(current_ent+1);
    if (current_iev < iev){
      tree->GetEntry(current_ent);
      return true;
    }
    else return false;
  }
  return true;
}

bool syncro_trees(TTree* mc_tree, TTree* sn_tree,
                  int&   mc_iev,  int&   sn_iev,
                  int&   mc_entry,int&   sn_entry,
                  int&   mc_entries,int& sn_entries){
  if (mc_iev == sn_iev){
    if (one_event_one_entry(sn_tree, sn_iev, sn_entry)) return true;
    else {
      if (!get_next_entry(mc_tree, mc_entry, mc_entries)) return false;
      // The while loop start with mc_iev > sn_iev
      while (mc_iev!=sn_iev || !one_event_one_entry(sn_tree, sn_iev, sn_entry)){
        if (mc_iev >= sn_iev){
          if (!get_next_entry(sn_tree, sn_entry, sn_entries)) return false;
        }
        else {
          if (!get_next_entry(mc_tree, mc_entry, mc_entries)) return false;
        }
      }
      return true;
    }
  }
  else if (mc_iev > sn_iev){
    while (mc_iev > sn_iev){
      if (!get_next_entry(sn_tree, sn_entry, sn_entries)) return false;
    }
    return syncro_trees(mc_tree, sn_tree, mc_iev, sn_iev, mc_entry, sn_entry, mc_entries, sn_entries);
  }
  else {
    while (mc_iev < sn_iev){
      if (!get_next_entry(mc_tree, mc_entry, mc_entries)) return false;
    }
    return syncro_trees(mc_tree, sn_tree, mc_iev, sn_iev, mc_entry, sn_entry, mc_entries, sn_entries);
  }
}

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


template <typename T>
void select_unique_elements_in_vector(std::vector<T>& vec) {
  std::map<T, size_t> count_map;
  for (const auto& elem : vec) {
    count_map[elem]++;
  }

  std::vector<T> unique_elements;
  for (const auto& [key, value] : count_map) {
    if (value == 1) {
      unique_elements.push_back(key);
    }
  }

  vec = unique_elements;
  return;
}


template <typename T>
std::vector<T> vector_intersection(const std::vector<T>& a, const std::vector<T>& b) {
  std::set<T> set_a(a.begin(), a.end());
  std::set<T> result_set;

  for (const T& elem : b) {
    if (set_a.count(elem)) result_set.insert(elem);
  }

  return std::vector<T>(result_set.begin(), result_set.end());
}

template <typename T>
TEntryList* get_entry_list(TTree* tree, const std::string& branch_name, const std::vector<T>& vec) {
  TEntryList* entry_list = new TEntryList();

  T value;
  tree->SetBranchAddress(branch_name.c_str(), &value);
  Long64_t n_entries = tree->GetEntries();

  for (long long i = 0; i < n_entries; ++i) {
    tree->GetEntry(i);
    if (std::find(vec.begin(), vec.end(), value) != vec.end()) {
      entry_list->Enter(i);
    }
  }

  return entry_list;
}

// --- HARD CODE HERE ----------------
const int max_nfiles = 8000; // Maximum number of files to analyze
int n_opdet = 480; // Number of optical detectors
const std::string ana_folder_name =  "./nue_cc/";
const std::string base_ana_file_name = "solar_ana_dune10kt_1x2x6_hist_";
const std::string visibility_file_name = "dunevis_fdhd_1x2x6_test_float.root"; // File with the visibility maps
const float pe_low = 0.;     // Plotting lower limit
const float pe_up  = 300.;   // Plotting upper limit
const float light_yield = 20000;
const float arapuca_pde = 0.03;
float LY_times_PDE = light_yield * arapuca_pde; // Light yield times the arapuca pde
const float min_visibility = 1.e-15;
// Leave zero if not using fiducial cuts
const float fiducial_cut = 40.; // Fiducial cut in cm (y=z=fiducial_cut)
const float x_cut = 40.; // x cut in cm
const bool verbose = true; // Print progress messages



template <typename T>
T weighted_average(const std::vector<T>& values, const std::vector<T>& weights) {

  T weighted_sum = 0.0;
  T total_weight = 0.0;

  for (size_t i = 0; i < values.size(); ++i) {
    weighted_sum += values[i] * weights[i];
    total_weight += weights[i];
  }

  return weighted_sum / total_weight;
}

TH2D* get_normalized_bidistribution(TH2D* h2_bidistr){
  TH2D* h2_bidistr_norm = new TH2D(Form("%s_norm", h2_bidistr->GetName()),
                                   Form("%s_norm;%s;%s", h2_bidistr->GetName(), h2_bidistr->GetXaxis()->GetTitle(), h2_bidistr->GetYaxis()->GetTitle()),
                                   h2_bidistr->GetNbinsX(), h2_bidistr->GetXaxis()->GetXmin(), h2_bidistr->GetXaxis()->GetXmax(),
                                   h2_bidistr->GetNbinsY(), h2_bidistr->GetYaxis()->GetXmin(), h2_bidistr->GetYaxis()->GetXmax());
  for(int idx_y=1; idx_y<=h2_bidistr->GetNbinsY(); idx_y++){
    TH1D* h1_proj = h2_bidistr->ProjectionX("h1_proj", idx_y, idx_y);
    double h1_proj_integral = h1_proj->Integral();
    
    if(h1_proj_integral > 0.){
      h1_proj->Scale(1./(h1_proj_integral*h1_proj->GetBinWidth(1)));

      for(int idx_x=1; idx_x<=h2_bidistr->GetNbinsX(); idx_x++){
        h2_bidistr_norm->SetBinContent(idx_x, idx_y, h1_proj->GetBinContent(idx_x));
      } 
    }
  }
  return h2_bidistr_norm;
}

void fm_distributions(){
  // --- INPUTS ---------------------------------------------------------------
  // Calibration fit + correction lambda and drift velocity
  // The latest could be retrieved from fhicls, actually
  TFile* calib_file = TFile::Open("./calibrator.root", "READ");
  TTree* calib_tree = static_cast<TTree*>(calib_file->Get("calib_tree"));
  Float_t calib_c = 0.;         Float_t calib_slope = 0.;
  Float_t drift_velocity = 0.0; Float_t corr_lambda = 0.0;
  calib_tree->SetBranchAddress("calib_c", &calib_c); 
  calib_tree->SetBranchAddress("calib_slope", &calib_slope);
  calib_tree->SetBranchAddress("drift_velocity", &drift_velocity);
  calib_tree->SetBranchAddress("corr_lambda", &corr_lambda);
  calib_tree->GetEntry(0);
  calib_file->Close();
  
  // --- VISIBILITY STUFF -----------------------------------------------------
  TFile* visibility_file = TFile::Open(visibility_file_name.c_str(), "READ");
  TH1D* hgrid[3] = {nullptr};
  hgrid[0] = (TH1D*)visibility_file->Get("photovisAr/hgrid0");
  hgrid[1] = (TH1D*)visibility_file->Get("photovisAr/hgrid1");
  hgrid[2] = (TH1D*)visibility_file->Get("photovisAr/hgrid2");
  
  TTree* tDimensions = (TTree*)visibility_file->Get("photovisAr/tDimensions");
  Double_t coor_dim[3] = {0.};
  tDimensions->SetBranchAddress("dimension", coor_dim);
  tDimensions->GetEntry(0);
  float tpc_min[3] = {float(coor_dim[0]), float(coor_dim[1]), float(coor_dim[2])}; // x,y,z
  tDimensions->GetEntry(1);
  float tpc_max[3] = {float(coor_dim[0]), float(coor_dim[1]), float(coor_dim[2])}; // x,y,z
  float vol_min[3] = {tpc_min[0]+x_cut, tpc_min[1]+fiducial_cut, tpc_min[2]+fiducial_cut};
  float vol_max[3] = {tpc_max[0]-x_cut, tpc_max[1]-fiducial_cut, tpc_max[2]-fiducial_cut};
  std::vector<int> cryo_to_tpc = GetCryoToTPCMap(hgrid, tpc_min, tpc_max); // Get the mapping from cryostat voxel to TPC voxel
  
  TTree* photoVisMap = (TTree*)visibility_file->Get("photovisAr/photoVisMap");
  const size_t n_entriesmap = photoVisMap->GetEntries();
  float opDet_visDirect[n_opdet];
  photoVisMap->SetBranchAddress("opDet_visDirect", &opDet_visDirect);
  std::vector<std::vector<float>> opDet_visMapDirect(n_entriesmap, std::vector<float>(n_opdet, 0.));

  for (size_t VisMapEntry = 0; VisMapEntry < n_entriesmap; ++VisMapEntry) {
    photoVisMap->GetEntry(VisMapEntry);
    std::memcpy(opDet_visMapDirect[VisMapEntry].data(), &opDet_visDirect[0], n_opdet * sizeof(float));
  }


  // --- PREPARE OUTPUT -------------------------------------------------------
  TFile* out_file = TFile::Open("fm_distributions.root", "RECREATE");
  TTree* tpc_pds_tree = new TTree("tpc_pds_tree", "tpc_pds_tree");
  Int_t ifile, iev;
  float charge, time_tpc, time_pds, y_reco, z_reco, e_reco;
  float x_true, y_true, z_true, e_true;
  std::vector<size_t>* flash_channels = nullptr;
  std::vector<float>* flash_pes = nullptr;
  std::vector<float>* flash_times = nullptr;
  std::vector<float>* reco_pes = nullptr;
  std::vector<float>* exp_phs = nullptr;

  tpc_pds_tree->Branch("ifile", &ifile);
  tpc_pds_tree->Branch("iev", &iev);
  tpc_pds_tree->Branch("charge", &charge);
  tpc_pds_tree->Branch("y_reco", &y_reco);
  tpc_pds_tree->Branch("z_reco", &z_reco);
  tpc_pds_tree->Branch("e_reco", &e_reco);
  tpc_pds_tree->Branch("time_tpc", &time_tpc);
  tpc_pds_tree->Branch("time_pds", &time_pds);
  tpc_pds_tree->Branch("flash_pes", &flash_pes);
  tpc_pds_tree->Branch("flash_channels", &flash_channels);
  tpc_pds_tree->Branch("flash_times", &flash_times);
  tpc_pds_tree->Branch("reco_pes", &reco_pes);
  tpc_pds_tree->Branch("exp_phs", &exp_phs);
  tpc_pds_tree->Branch("x_true", &x_true);
  tpc_pds_tree->Branch("y_true", &y_true);
  tpc_pds_tree->Branch("z_true", &z_true);
  tpc_pds_tree->Branch("e_true", &e_true);
  
  TH1D* h_exp = new TH1D("h_exp",
                         Form("%s;%s;%s","h_exp","Expected #Pe","Counts"),
                         300, pe_low, pe_up);
  
  TH1D* h_expreco = new TH1D("h_expreco",
                             Form("%s;%s;%s","h_expreco","Expected #Pe", "Counts"),
                             300, pe_low, pe_up);

  TH2D* h2_exp_reco = new TH2D("h2_exp_reco",
                               Form("%s;%s;%s", "h2_exp_reco", "Reco #Pe", "Expected #Pe"),
                               300*5, pe_low, pe_up*5, 300*5, pe_low, pe_up*5);

  // --- EXTRA VARIABLES ------------------------------------------------------
  Float_t vertex_coor[3] = {0.};
  
  // --- LOOP OVER ANA FILES ---------------------------------------------------
  std::filesystem::path ana_folder(ana_folder_name);
  int nfile_in_folder = std::distance(std::filesystem::directory_iterator(ana_folder), std::filesystem::directory_iterator{});
  int nfile_to_analyze = std::min(nfile_in_folder, max_nfiles);
  int nfile_analyzed = 0;
  int idx_file = 0;
  while(nfile_analyzed < nfile_to_analyze){
    // --- ANA STUFF -----------------------------------------------------------
    // std::string ana_file_name = std::string(ana_folder_name+base_ana_file_name+idx_file+".root"); 
    std::string ana_file_name = "./nue_cc/SolarNuAna_marley.root";
    // Check whether both files exist
    idx_file++;
    if(!std::filesystem::exists(ana_file_name)) continue;
    ifile = idx_file; nfile_analyzed++;
    /* if (idx_file % 10 == 0) */ std::cout << nfile_analyzed <<"--"<< ana_file_name << std::endl;

    TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
    TTree* solarnu_tree = static_cast<TTree*>(ana_file->Get("solarnuana/SolarNuAnaTree"));
    int sn_iev = 0; int sn_flag = 0;
    float TrueY = 0.; float RecoY = 0.; float RecoZ = 0.;
    solarnu_tree->SetBranchAddress("Event", &sn_iev);
    solarnu_tree->SetBranchAddress("Flag", &sn_flag);
    solarnu_tree->SetBranchAddress("RecoY", &RecoY);
    solarnu_tree->SetBranchAddress("RecoZ", &RecoZ);
    
    TTree* mctruth_tree = static_cast<TTree*>(ana_file->Get("solarnuana/MCTruthTree"));
    int mc_iev = 0; int mc_flag = 0;
    float TrueX = 0.; float TrueZ = 0.; float TrueE = 0.;
    std::vector<int>* EDepElectron = nullptr;
    std::vector<float>* OpHitPes = nullptr;
    std::vector<float>* OpHitChannels = nullptr;
    std::vector<float>* OpHitTimes = nullptr;
    mctruth_tree->SetBranchAddress("Event", &mc_iev);
    mctruth_tree->SetBranchAddress("Flag", &mc_flag);
    mctruth_tree->SetBranchAddress("SignalParticleX", &TrueX);
    mctruth_tree->SetBranchAddress("SignalParticleY", &TrueY);
    mctruth_tree->SetBranchAddress("SignalParticleZ", &TrueZ);
    mctruth_tree->SetBranchAddress("SignalParticleE", &TrueE);
    mctruth_tree->SetBranchAddress("TSignalElectronDepList", &EDepElectron);
    mctruth_tree->SetBranchAddress("OpHitPE", &OpHitPes);
    mctruth_tree->SetBranchAddress("OpHitChannel", &OpHitChannels);
    mctruth_tree->SetBranchAddress("OpHitTime", &OpHitTimes);


    int sn_entries = solarnu_tree->GetEntries(); int this_sn_entry = 0;
    int mc_entries = mctruth_tree->GetEntries(); int this_mc_entry = 0;
    bool out_of_sync = false; bool end_tree = false; int this_sn_iev = 0;
    int sn_entry = 0;

    std::vector<int> mc_events = treeBranch_to_vector<int>(mctruth_tree, "Flag");
    std::vector<int> sn_events = treeBranch_to_vector<int>(solarnu_tree, "Flag");
    select_unique_elements_in_vector(sn_events);
    std::vector<int> common_events = vector_intersection(mc_events, sn_events);
    TEntryList* mc_entry_list = get_entry_list<int>(mctruth_tree, "Flag", common_events);
    TEntryList* sn_entry_list = get_entry_list<int>(solarnu_tree, "Flag", common_events);
    mctruth_tree->SetEntryList(mc_entry_list);
    solarnu_tree->SetEntryList(sn_entry_list);
    sn_entries = sn_entry_list->GetN();
    mc_entries = mc_entry_list->GetN();


    for (int mc_entry = 0; mc_entry < mc_entries; mc_entry++) {
      mctruth_tree->GetEntry(mc_entry_list->GetEntry(mc_entry));
      solarnu_tree->GetEntry(sn_entry_list->GetEntry(sn_entry));
      if (sn_iev != mc_iev) std::cout << mc_entry << " out of sync " << std::endl;

      // syncro_trees(mctruth_tree, solarnu_tree, mc_iev, sn_iev, mc_entry, sn_entry, mc_entries, sn_entries);
      sn_entry++;

      y_reco = RecoY; z_reco = RecoZ;
      x_true = TrueX; y_true = TrueY; z_true = TrueZ; e_true = TrueE;
      

      vertex_coor[0] = TrueX; vertex_coor[1] = TrueY; vertex_coor[2] = TrueZ;
      if(!isInFiducialVolume(vertex_coor, vol_min, vol_max, x_cut)) continue;
      int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);
      if (tpc_index < 0 || tpc_index >= int(cryo_to_tpc.size())) {
        std::cout << "Event " << iev << " is not in the fiducial volume!" << std::endl;
        std::cout << vertex_coor[0] << " " << vertex_coor[1] << " " << vertex_coor[2] << std::endl;
        continue; // Skip this entry if TPC index is invalid
      }


      vertex_coor[0] = TrueX; vertex_coor[1] = RecoY; vertex_coor[2] = RecoZ;
      if(!isInFiducialVolume(vertex_coor, vol_min, vol_max, x_cut)) continue;
      tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);
      if (tpc_index < 0 || tpc_index >= int(cryo_to_tpc.size())) {
        std::cout << "Event " << iev << " is not in the fiducial volume!" << std::endl;
        std::cout << vertex_coor[0] << " " << vertex_coor[1] << " " << vertex_coor[2] << std::endl;
        continue; // Skip this entry if TPC index is invalid
      }

      iev = mc_iev;
      charge = std::accumulate(EDepElectron->begin(), EDepElectron->end(), 0.0);
      time_tpc = (TrueX) / drift_velocity;
      e_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, time_tpc, charge);
      std::vector<float> ophit_pes = std::move(*OpHitPes);
      std::vector<float> ophit_channels = std::move(*OpHitChannels);
      std::vector<float> ophit_times = std::move(*OpHitTimes);

      if (ophit_pes.size() != ophit_channels.size() ||
          ophit_pes.size() != ophit_times.size()) {
        std::cerr << "Error: OpHit vectors have different sizes!" << std::endl;
        continue; // Skip this entry if sizes do not match
      }
      if (ophit_pes.empty()) time_pds = -1e3;
      else                   time_pds = weighted_average(ophit_times, ophit_pes);

      std::unordered_map<size_t, std::vector<size_t>> channel_hit_map;
      for (size_t i = 0; i < ophit_channels.size(); i++) {
        channel_hit_map[size_t(ophit_channels[i])].push_back(i);
      }
        
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        float voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];
        float reco_pe = 0.; float exp_ph = 0.;
        exp_ph = e_reco*LY_times_PDE*voxel_vis; // Expected photons
        if(exp_ph==0.) exp_ph = e_reco*LY_times_PDE*min_visibility; 

        auto it = channel_hit_map.find(idx_opdet);
        // If recontructed ----------------------------------------------------
        if (it != channel_hit_map.end()) {
          for (const auto& idx_hit : it->second) {
            reco_pe += ophit_pes[idx_hit];
          }
          flash_pes->push_back(reco_pe);
          flash_channels->push_back(idx_opdet);
          flash_times->push_back(ophit_times[it->second[0]]);

          h_exp->Fill(exp_ph);
          h_expreco->Fill(exp_ph);
          h2_exp_reco->Fill(reco_pe, exp_ph);
        }
        else {
          h_exp->Fill(exp_ph);
        } // if not reconstructed
        reco_pes->push_back(reco_pe); exp_phs->push_back(exp_ph);
      } // opdet loop
      tpc_pds_tree->Fill();
      flash_channels->clear(); flash_pes->clear(); flash_times->clear(); reco_pes->clear(); exp_phs->clear();
    } // SolarNuAna loop
    ana_file->Close();
  }

  // --- EXTRA OUTPUTS --------------------------------------------------------
  TEfficiency* he_hit_prob = new TEfficiency(*h_expreco,*h_exp);
  he_hit_prob->SetTitle("Hit Probability;Expected #Pe;Reconstruction Probability");
  he_hit_prob->SetName("he_hit_prob");

  TH2D* h2_exp_reco_norm = get_normalized_bidistribution(h2_exp_reco);

  TTree* pds_tree = new TTree("pds_tree", "pds_tree");
  pds_tree->Branch("n_opdet", &n_opdet);
  pds_tree->Branch("LY_times_PDE", &LY_times_PDE);
  pds_tree->Fill();

  // --- WRITE OUTPUT ---------------------------------------------------------
  out_file->cd();
  tpc_pds_tree->Write();
  pds_tree->Write();
  h_exp->Write();
  h_expreco->Write();
  he_hit_prob->Write();
  h2_exp_reco->Write();
  h2_exp_reco_norm->Write();
  out_file->Close();

  return;
}
