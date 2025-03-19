
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <filesystem>  // C++17

#include <THnSparse.h>
#include <TTree.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>
#include <TTreeReader.h>

#include "TEfficiency.h"
#include "file_list.hpp"
#include "Utils.hpp"



// --- HARD CODE HERE ----------------
int rebinning = 4;

int nfile_to_analyze = 1e3;
int max_nfiles = 1e4; // Maximum number of files to analyze

const size_t n_opdet = 480; // 480
float light_yield = 20000;
float arapuca_pde = 0.03;

double min_visibility = 1.e-60;
double hit_threshold = 1.5; // Will integrate Poisson [0, hit_threshold]

TString buildmu_output = "./out_buildmu.root";
TString visibility_file_name = "./dunevis_fdhd_1x2x6_test.root";
std::string ana_folder_name = "/exp/dune/data/users/fgalizzi/prod_eminus/seed_0/ana/"; // Folder where the ana files
std::string base_ana_file_name = "solar_ana_dune10kt_1x2x6_hist_"; // Base name of the ana files
// -----------------------------------

void compute_likelihood(){
  // --- INPUT ----------------------------------------------------------------
  TFile* exp_reco_file = TFile::Open(buildmu_output, "READ");
  TH2D* h2_exp_reco = nullptr;
  exp_reco_file->GetObject("h2_exp_reco", h2_exp_reco);

  // h2_exp_reco->Rebin2D(rebinning, rebinning);

  // for each bin in the y-axis of h2_exp_reco project on the x-axis and normalize
  // the projection to have area=1
  for(int idx_y=1; idx_y<=h2_exp_reco->GetNbinsY(); idx_y++){
    TH1D* h1_proj = h2_exp_reco->ProjectionX("h1_proj", idx_y, idx_y);
    h1_proj->Scale(1./h1_proj->Integral());
    for(int idx_x=1; idx_x<=h2_exp_reco->GetNbinsX(); idx_x++){
      h2_exp_reco->SetBinContent(idx_y, idx_x, h1_proj->GetBinContent(idx_x));
    }
  }

  TEfficiency* he_ProbHit_ExpPe = nullptr;
  exp_reco_file->GetObject("h_Expected_Pe_clone", he_ProbHit_ExpPe);

  // --- VISIBILITY STUFF -------------------------------------------------------
  TFile* ff = TFile::Open("dunevis_fdhd_1x2x6_test.root", "READ");
  TTree* photoVisMap = (TTree*)ff->Get("photovisAr/photoVisMap");
  TH1D* hgrid[3] = {nullptr};
  hgrid[0] = (TH1D*)ff->Get("photovisAr/hgrid0");
  hgrid[1] = (TH1D*)ff->Get("photovisAr/hgrid1");
  hgrid[2] = (TH1D*)ff->Get("photovisAr/hgrid2");
  
  TTree* tDimensions = (TTree*)ff->Get("photovisAr/tDimensions");
  double coor_dim[3] = {0.};
  tDimensions->SetBranchAddress("dimension", coor_dim);
  tDimensions->GetEntry(0);
  double tpc_min[3] = {coor_dim[0], coor_dim[1], coor_dim[2]}; // x,y,z
  tDimensions->GetEntry(1);
  double tpc_max[3] = {coor_dim[0], coor_dim[1], coor_dim[2]}; // x,y,z
  std::vector<int> cryo_to_tpc = GetCryoToTPCMap(hgrid, tpc_min, tpc_max); // Get the mapping from cryostat voxel to TPC voxel
  std::cout << "cryo_to_tpc size " << cryo_to_tpc.size() << std::endl;

  double opDet_visDirect[n_opdet];
  double opDet_visReflct[n_opdet];
  const size_t n_entriesmap = photoVisMap->GetEntries();
  std::vector<std::vector<float>> opDet_visMapDirect(n_entriesmap, std::vector<float>(n_opdet, 0.)); // Map to store visibility for each opdet
  std::vector<std::vector<float>> opDet_visMapReflct(n_entriesmap, std::vector<float>(n_opdet, 0.)); // Map to store visibility for each opdet
  photoVisMap->SetBranchStatus("*", 0); // Disable all branches
  photoVisMap->SetBranchStatus("opDet_visDirect", 1);
  photoVisMap->SetBranchStatus("opDet_visReflct", 1);
  photoVisMap->SetBranchAddress("opDet_visDirect", opDet_visDirect);
  photoVisMap->SetBranchAddress("opDet_visReflct", opDet_visReflct);
  photoVisMap->Draw(">>elist", "", "entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  photoVisMap->SetEntryList(elist);
  std::cout << "Filling maps..." << std::endl;
  for (int i = 0; i < photoVisMap->GetEntries(); ++i) {
    if (i % 1000 == 0)
        std::cout << i << " / " << n_entriesmap << "\r" << std::flush;
    photoVisMap->GetEntry(i);
    for (size_t j = 0; j < n_opdet; ++j) {
      opDet_visMapDirect[i][j] = static_cast<float>(opDet_visDirect[j]);
      opDet_visMapReflct[i][j] = static_cast<float>(opDet_visReflct[j]);
    }
  }
  std::cout << "Done filling maps..." << std::endl;
  
  // --- PDF -------------------------------------------------------------------
  size_t fake_counter = 0;
  size_t evt_counter = 0;

  // --- LOOP OVER ANA FILES ---------------------------------------------------
  int nfile_analyzed = 0;
  int idx_file = 0;
  std::vector<double> vec_Expected_Ophit_OpDet(n_opdet, 0.); // Vector to store expected OpHits per OpDet
  std::vector<double> vec_Reco_Ophit_OpDet(n_opdet, 0.); // Vector to store reconstructed OpHits per OpDet
  while(nfile_analyzed < nfile_to_analyze && idx_file < max_nfiles){
    // --- ANA STUFF -----------------------------------------------------------
    std::string ana_file_name = std::string(ana_folder_name+base_ana_file_name+idx_file+".root"); 
    // Check whether file exists
    idx_file++;
    if(!std::filesystem::exists(ana_file_name)){
      continue;
    }
    std::cout << nfile_analyzed << "- Analyzing file: " << ana_file_name << std::endl;
    nfile_analyzed++;

    TFile* ana_file = TFile::Open(ana_file_name.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(ana_file->Get("solarnuana/MCTruthTree"));
    TTreeReader treeReader(tree);

    TTreeReaderValue<float> E_true(treeReader, "SignalParticleE");
    TTreeReaderValue<float> x_true(treeReader, "SignalParticleX");
    TTreeReaderValue<float> y_true(treeReader, "SignalParticleY");
    TTreeReaderValue<float> z_true(treeReader, "SignalParticleZ");
    TTreeReaderValue<int> event_true(treeReader, "Event");
    TTreeReaderValue<std::vector<float>> OpHitPes(treeReader, "OpHitPE");
    TTreeReaderValue<std::vector<float>> OpHitChannels(treeReader, "OpHitChannel");
  
    // --- LOOP OVER TREE -----------------------------------------------------
    size_t idx_entry = 0;
    double old_E_true = 0.; int old_tpc_index = -1;
    while (treeReader.Next()) {
      double vertex_coor[3] = {*x_true, *y_true, *z_true};
      int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);
      
      if (idx_entry == 0){
        old_E_true = *E_true;
        old_tpc_index = tpc_index;
        idx_entry++;
        continue;
      }
      
      if (tpc_index < 0 || old_tpc_index < 0) {
        std::cout << "New/old event " << *event_true << " has no valid TPC index!" << std::endl;
        old_E_true = *E_true;
        old_tpc_index = tpc_index;
        continue; // Skip this entry if TPC index is invalid
      }
      
      double exp_ph, exp_ph_fake;
      
      double Fq = 0.;
      double Fq_fake = 0.;


      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        double voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet] + opDet_visMapReflct[tpc_index][idx_opdet];
        double voxel_vis_fake = opDet_visMapDirect[old_tpc_index][idx_opdet] + opDet_visMapReflct[old_tpc_index][idx_opdet];

        exp_ph = (*E_true)*light_yield*voxel_vis*arapuca_pde;
        exp_ph_fake = old_E_true*light_yield*voxel_vis_fake*arapuca_pde;

        // Compute a single term of the Fq summatory
        double term = 0.;
        double term_fake = 0.;
        double P_hit_mu = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(exp_ph));
        double P_hit_mu_fake = he_ProbHit_ExpPe->GetEfficiency(he_ProbHit_ExpPe->FindFixBin(exp_ph_fake));

        auto it = std::find((*OpHitChannels).begin(), (*OpHitChannels).end(), float(idx_opdet));
        size_t idx_hit = std::distance((*OpHitChannels).begin(), it);

        if(idx_hit != (*OpHitChannels).size()){ 
          term = P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin((*OpHitPes)[idx_hit], exp_ph));
          term_fake = P_hit_mu_fake*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin((*OpHitPes)[idx_hit], exp_ph_fake));
        } else {
          term = 1. - P_hit_mu;
          term_fake = 1. - P_hit_mu_fake;
        }

        Fq += -log(term);
        Fq_fake += -log(term_fake);
      } // end loop over opdets
      
      if(Fq_fake < Fq) fake_counter++;
      evt_counter++;
      old_E_true = *E_true;
      old_tpc_index = tpc_index;
    } // end loop over tree
    
    ana_file->Close();
    std::cout << "Fake counter: " << fake_counter << "/" << evt_counter << std::endl;
  } // end loop over ana files

}
