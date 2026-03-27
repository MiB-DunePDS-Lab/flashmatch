#include "RtypesCore.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TEntryList.h"
#include "TTree.h"
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "Utils.hpp"


void fm_distributions(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLcconfigs f = load_ana_config("./config.json");
  std::string input_dir            = f.input_dir;
  std::string visibility_file_name = f.visibility_file_name;
  int max_nfiles                   = f.max_nfiles;
  float pe_low                     = f.pe_low;
  float pe_up                      = f.pe_up;
  float light_yield                = f.light_yield;
  float arapuca_pde                = f.arapuca_pde;
  const size_t n_opdet             = f.n_opdet;
  float min_visibility             = f.min_visibility;
  float fiducial_cut               = f.yz_fiducial_cut;
  float x_cut                      = f.x_fiducial_cut;
  float LY_times_PDE               = light_yield * arapuca_pde; // Light yield times the arapuca pde


  // --- INPUTS ---------------------------------------------------------------
  // Calibration fit + correction lambda and drift velocity
  // The latest could be retrieved from fhicls, actually
  TFile* calib_file = TFile::Open((input_dir+"MLL_Calibrator.root").c_str(), "READ");
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
  TFile* out_file = TFile::Open((input_dir+"MLL_Distributions.root").c_str(), "RECREATE");
  TTree* tpc_pds_tree = new TTree("tpc_pds_tree", "tpc_pds_tree");
  Int_t ifile, iev;
  float charge, time_tpc, time_pds, my_x_reco, x_reco, y_reco, z_reco, e_reco;
  float x_true, y_true, z_true, e_true;
  std::vector<float>* reco_pes          = nullptr;
  std::vector<float>* exp_phs           = nullptr;
  std::vector<int>* opdets           = nullptr;

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
  std::string sample_dir = input_dir+"files/";
  std::vector<std::string> ana_files = get_list_of_files_in_folder(sample_dir, ".root");
  int nfile_to_analyze = std::min(int(ana_files.size()), max_nfiles);
  int nfile_analyzed = 0; int idx_file = 0;
  while(nfile_analyzed < nfile_to_analyze){
    // --- ANA STUFF -----------------------------------------------------------
    std::string ana_file_name = sample_dir+ana_files[idx_file];
    idx_file++;
    if(!std::filesystem::exists(ana_file_name)) continue;
    ifile = idx_file; nfile_analyzed++;
    if (idx_file % 10 == 0) std::cout <<nfile_analyzed<<"--"<< ana_file_name << "\r" << std::flush;

    TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
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
    solarnu_tree->SetBranchAddress("MatchedOpFlashPur",       &MatchedOpFlashPurity);
    solarnu_tree->SetBranchAddress("MatchedOpFlashCorrectly",  &MatchedOpFlashCorrectly);
    solarnu_tree->SetBranchAddress("MatchedOpFlashPEperOpDet", &MatchedOpFlashPEperOpDet);
    solarnu_tree->SetBranchAddress("Charge",                   &Charge);
    solarnu_tree->SetBranchAddress("Time",                     &tpc_time);
    
    int sn_entries = solarnu_tree->GetEntries();
    // for (int sn_entry = 0; sn_entry < sn_entries; sn_entry++) {
      // solarnu_tree->GetEntry(sn_entry);
    for (auto& idx_entry : MaxChargeIndxs){
      solarnu_tree->GetEntry(idx_entry);
      if (!MatchedOpFlashCorrectly) continue;
      iev = sn_iev;

      x_reco = MatchedOpFlashRecoX; y_reco = RecoY; z_reco = RecoZ;
      x_true = TrueX; y_true = TrueY; z_true = TrueZ; e_true = TrueE;
      // if (MatchedOpFlashPurity<0.9) continue;
   

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

      charge = Charge;
      double dt = tpc_time - MatchedOpFlashTime;
      my_x_reco = dt*drift_velocity;
      e_reco = give_me_Ereco(calib_c, calib_slope, corr_lambda, dt, charge);
      time_pds = MatchedOpFlashTime;
      time_tpc = tpc_time;

        
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        float voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];
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
  }

  // --- EXTRA OUTPUTS --------------------------------------------------------
  TEfficiency* he_hit_prob = new TEfficiency(*h_expreco,*h_exp);
  he_hit_prob->SetTitle("Hit Probability;Expected #Pe;Reconstruction Probability");
  he_hit_prob->SetName("he_hit_prob");

  // --- WRITE OUTPUT ---------------------------------------------------------
  out_file->cd();
  tpc_pds_tree->Write();
  h_exp->Write();
  h_expreco->Write();
  he_hit_prob->Write();
  h2_exp_reco->Write();
  out_file->Close();

  return;
}
