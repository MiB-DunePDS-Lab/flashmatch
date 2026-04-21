#include "RtypesCore.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TEntryList.h"
#include "TTree.h"
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "Utils.hpp"


void fm_distributions(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLConfigs mll_conf = load_ana_config("./configs/ana_config.json");
  std::string ana_file_name        = mll_conf.ana_file_name;
  std::string visibility_dir       = mll_conf.visibility_dir;
  std::string sample_config_file   = mll_conf.sample_config_file;
  float pe_up                      = mll_conf.pe_up;
  float light_yield                = mll_conf.light_yield;
  float arapuca_pde                = mll_conf.arapuca_pde;
  float min_visibility             = mll_conf.min_visibility;
  float fiducial_cut               = mll_conf.yz_fiducial_cut;
  float x_cut                      = mll_conf.x_fiducial_cut;
  float LY_times_PDE               = light_yield * arapuca_pde; // Light yield times the arapuca pde

  SampleConfigs sample_conf = load_sample_config("./configs/"+sample_config_file);
  std::string input_dir            = sample_conf.input_dir;
  std::string geom_identifier      = sample_conf.geom_identifier;

  TString visibility_file_name     = TString(visibility_dir+"dunevis_"+geom_identifier+".root");
  
  DuneGeom geom = load_dune_geom("./configs/"+geom_identifier+".json");
  const size_t n_opdet = geom.n_opdet;
  float anode_x = geom.anode_x;

  // --- INPUTS ---------------------------------------------------------------
  // Calibration fit + correction lambda and drift velocity
  // The latest could be retrieved from fhicls, actually
  TFile* calib_file = TFile::Open((input_dir+"MLL_Calibrator_"+geom_identifier+".root").c_str(), "READ");
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
  TFile* visibility_file = TFile::Open(visibility_file_name, "READ");
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
  TFile* out_file = TFile::Open((input_dir+"MLL_Distributions_"+geom_identifier+".root").c_str(), "RECREATE");
  TTree* anode_tree = new TTree("anode_tree", "anode_tree");
  anode_tree->Branch("anode_x", &anode_x);
  anode_tree->Fill();

  TTree* tpc_pds_tree = new TTree("tpc_pds_tree", "tpc_pds_tree");
  Int_t ifile, iev;
  float charge, time_tpc, time_pds, my_x_reco, x_reco, y_reco, z_reco, e_reco;
  float x_true, y_true, z_true, e_true;
  std::vector<float>* reco_pes          = nullptr;
  std::vector<float>* exp_phs           = nullptr;
  std::vector<int>* opdets              = nullptr;

  tpc_pds_tree->Branch("ifile",     &ifile);
  tpc_pds_tree->Branch("iev",       &iev);
  tpc_pds_tree->Branch("charge",    &charge);
  tpc_pds_tree->Branch("my_x_reco", &my_x_reco);
  tpc_pds_tree->Branch("x_reco",    &x_reco);
  tpc_pds_tree->Branch("y_reco",    &y_reco);
  tpc_pds_tree->Branch("z_reco",    &z_reco);
  tpc_pds_tree->Branch("e_reco",    &e_reco);
  tpc_pds_tree->Branch("time_tpc",  &time_tpc);
  tpc_pds_tree->Branch("time_pds",  &time_pds);
  tpc_pds_tree->Branch("reco_pes",  &reco_pes);
  tpc_pds_tree->Branch("exp_phs",   &exp_phs);
  tpc_pds_tree->Branch("opdets",    &opdets);
  tpc_pds_tree->Branch("x_true",    &x_true);
  tpc_pds_tree->Branch("y_true",    &y_true);
  tpc_pds_tree->Branch("z_true",    &z_true);
  tpc_pds_tree->Branch("e_true",    &e_true);
  for(size_t i=0; i<n_opdet; i++){
    opdets->push_back(i);
  }


  // --- EXTRA VARIABLES ------------------------------------------------------
  Float_t vertex_coor[3] = {0.};

  // --- LOOP OVER ANA FILES ---------------------------------------------------
  // --- ANA STUFF -----------------------------------------------------------
  ana_file_name = input_dir+ana_file_name;
  if(!std::filesystem::exists(ana_file_name)) {
    printf("File %s does not exist!\n", ana_file_name.c_str());
    return;
  }

  TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
  TTree* config_tree = static_cast<TTree*>(ana_file->Get("solarnuana/ConfigTree"));
  float OpFlashAlgoPE = 0.;
  config_tree->SetBranchAddress("OpFlashAlgoPE", &OpFlashAlgoPE);
  config_tree->GetEntry(0);
  if (OpFlashAlgoPE <= 0.) OpFlashAlgoPE = 1.5;
  std::vector<double> bin_lower_edges;
  bin_lower_edges.push_back(0.);
  double lower_bound = 0.1; double upper_bound = pe_up*5.;
  double bin_width = 0.8/1.05;
  while (lower_bound < upper_bound) {
    if (lower_bound <= OpFlashAlgoPE) {
      bin_lower_edges.push_back(lower_bound);
      lower_bound += 0.1; // Finer bins below the OpFlashAlgoPE threshold
    }
    else {
      lower_bound += bin_width;
      bin_lower_edges.push_back(lower_bound);
      bin_width *= 1.05; // Increase the bin width by 5% each time
    }
  }

  out_file->cd();
  TH1D* h_exp = new TH1D("h_exp",
                         Form("%s;%s;%s","h_exp","Expected #Pe","Counts"),
                         bin_lower_edges.size()-1, bin_lower_edges.data());

  TH1D* h_expreco = new TH1D("h_expreco",
                             Form("%s;%s;%s","h_expreco","Expected #Pe", "Counts"),
                             bin_lower_edges.size()-1, bin_lower_edges.data());


  TH2D* h2_exp_reco = new TH2D("h2_exp_reco",
                               Form("%s;%s;%s", "h2_exp_reco", "Reco #Pe", "Expected #Pe"),
                               bin_lower_edges.size()-1, bin_lower_edges.data(),
                               bin_lower_edges.size()-1, bin_lower_edges.data());
  
  // TGraph for Poisson hit probability vs expected PE. Same x-axis binning as h_exp and h_expreco
  TGraph* g_poisson_hitprob = new TGraph();
  g_poisson_hitprob->SetName("g_poisson_hitprob");
  g_poisson_hitprob->SetTitle("Poisson Hit Probability;Expected #Pe;Hit Probability");
  for (size_t i=1; i<=h_exp->GetNbinsX(); i++) {
    g_poisson_hitprob->SetPoint(g_poisson_hitprob->GetN(), h_exp->GetBinCenter(i), 1. - TMath::Poisson(0, h_exp->GetBinCenter(i)));
  }

  ana_file->cd();
  TTree* solarnu_tree = static_cast<TTree*>(ana_file->Get("solarnuana/SolarNuAnaTree"));
  std::vector<size_t> MaxChargeIndxs = take_max_charge_indices(solarnu_tree, "Event", "Charge");

  int sn_iev = 0; int sn_flag = 0;
  float TrueX = 0.; float TrueY = 0.; float TrueZ = 0.;
  float TrueE = 0.;
  float MatchedOpFlashRecoX = 0.;
  float MatchedOpFlashPurity = 0.;
  // float MatchedOpFlashRecoY = 0.; float MatchedOpFlashRecoZ = 0.;
  float MatchedOpFlashTime = 0.;
  bool MatchedOpFlashCorrectly = false;
  std::vector<double>* MatchedOpFlashPEperOpDet = nullptr;
  float Charge = 0.;
  float tpc_time = 0.;
  float RecoY = 0.; float RecoZ = 0.;

  solarnu_tree->SetBranchAddress("Event", &sn_iev);
  solarnu_tree->SetBranchAddress("Flag",  &sn_flag);
  solarnu_tree->SetBranchAddress("RecoY", &RecoY);
  solarnu_tree->SetBranchAddress("RecoZ", &RecoZ);

  solarnu_tree->SetBranchAddress("SignalParticleE",          &TrueE);
  solarnu_tree->SetBranchAddress("SignalParticleX",          &TrueX);
  solarnu_tree->SetBranchAddress("SignalParticleY",          &TrueY);
  solarnu_tree->SetBranchAddress("SignalParticleZ",          &TrueZ);
  solarnu_tree->SetBranchAddress("MatchedOpFlashRecoX",      &MatchedOpFlashRecoX);
  // solarnu_tree->SetBranchAddress("MatchedOpFlashRecoY",      &MatchedOpFlashRecoY);
  // solarnu_tree->SetBranchAddress("MatchedOpFlashRecoZ",      &MatchedOpFlashRecoZ);
  solarnu_tree->SetBranchAddress("MatchedOpFlashTime",       &MatchedOpFlashTime);
  solarnu_tree->SetBranchAddress("MatchedOpFlashPur",        &MatchedOpFlashPurity);
  solarnu_tree->SetBranchAddress("MatchedOpFlashCorrectly",  &MatchedOpFlashCorrectly);
  solarnu_tree->SetBranchAddress("MatchedOpFlashPEperOpDet", &MatchedOpFlashPEperOpDet);
  solarnu_tree->SetBranchAddress("Charge",                   &Charge);
  solarnu_tree->SetBranchAddress("Time",                     &tpc_time);

  int sn_entries = solarnu_tree->GetEntries();
  // for (int sn_entry = 0; sn_entry < sn_entries; sn_entry++) {
  // solarnu_tree->GetEntry(sn_entry);

  size_t last_index = MaxChargeIndxs[MaxChargeIndxs.size() - 1];
  for (auto& idx_entry : MaxChargeIndxs){
    solarnu_tree->GetEntry(idx_entry);
    if (idx_entry % 100 == 0) std::cout << idx_entry <<"/"<< last_index << "\r" << std::flush;
    if (!MatchedOpFlashCorrectly) continue;
    iev = sn_iev;

    x_reco = MatchedOpFlashRecoX; y_reco = RecoY; z_reco = RecoZ;
    x_true = TrueX; y_true = TrueY; z_true = TrueZ; e_true = TrueE;
    // if (MatchedOpFlashPurity<0.9) continue;


    vertex_coor[0] = TrueX; vertex_coor[1] = TrueY; vertex_coor[2] = TrueZ;
    if(!isInFiducialVolume(vertex_coor, vol_min, vol_max, x_cut, anode_x)) continue;
    int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);
    if (tpc_index < 0 || tpc_index >= int(cryo_to_tpc.size())) {
      std::cout << "Event " << iev << " is not in the fiducial volume!" << std::endl;
      std::cout << vertex_coor[0] << " " << vertex_coor[1] << " " << vertex_coor[2] << std::endl;
      continue; // Skip this entry if TPC index is invalid
    }

    vertex_coor[0] = TrueX; vertex_coor[1] = RecoY; vertex_coor[2] = RecoZ;
    if(!isInFiducialVolume(vertex_coor, vol_min, vol_max, x_cut, anode_x)) continue;
    tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);
    if (tpc_index < 0 || tpc_index >= int(cryo_to_tpc.size())) {
      std::cout << "Event " << iev << " is not in the fiducial volume!" << std::endl;
      std::cout << vertex_coor[0] << " " << vertex_coor[1] << " " << vertex_coor[2] << std::endl;
      continue; // Skip this entry if TPC index is invalid
    }

    charge = Charge;
    double dt = tpc_time - MatchedOpFlashTime;
    my_x_reco = dt*drift_velocity;
    e_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, dt, charge);
    time_pds = MatchedOpFlashTime;
    time_tpc = tpc_time;


    for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
      float voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];
      if (idx_opdet<geom.opdets_per_plane[0].first or idx_opdet>geom.opdets_per_plane[0].second) {
        voxel_vis = 0.;
      }
      float reco_pe = 0.; float exp_ph = 0.;
      exp_ph = e_reco*LY_times_PDE*voxel_vis; // Expected photons
      if(exp_ph==0.) exp_ph = e_reco*LY_times_PDE*min_visibility; 

      // If recontructed ----------------------------------------------------
      if (MatchedOpFlashPEperOpDet->at(idx_opdet)>0.) {
        reco_pe = MatchedOpFlashPEperOpDet->at(idx_opdet);

        h_exp->Fill(exp_ph);
        h_expreco->Fill(exp_ph);
        h2_exp_reco->Fill(reco_pe, exp_ph);
        // if (exp_ph < 50 && reco_pe > 200) {
        //   std::cout << TrueX << " " << TrueY << " " << TrueZ <<" "<< TrueE <<" "<< e_reco <<" "<< e_reco/TrueE << std::endl;
        // }
      }
      else {
        h_exp->Fill(exp_ph);
      } // if not reconstructed
      reco_pes->push_back(reco_pe); exp_phs->push_back(exp_ph);
    } // opdet loop
    tpc_pds_tree->Fill();
    reco_pes->clear(); exp_phs->clear();
  } // SolarNuAna loop
  ana_file->Close();

  // --- EXTRA OUTPUTS --------------------------------------------------------
  TEfficiency* he_hit_prob = new TEfficiency(*h_expreco,*h_exp);
  he_hit_prob->SetTitle("Hit Probability;Expected #Pe;Reconstruction Probability");
  he_hit_prob->SetName("he_hit_prob");

  // --- WRITE OUTPUT ---------------------------------------------------------
  out_file->cd();
  anode_tree->Write();
  tpc_pds_tree->Write();
  h_exp->Write();
  h_expreco->Write();
  he_hit_prob->Write();
  h2_exp_reco->Write();
  g_poisson_hitprob->Write();
  out_file->Close();

  return;
}
