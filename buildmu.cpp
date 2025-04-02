#include "Rtypes.h"
#include "THnBase.h"
#include "TProfile.h"
#include "TString.h"
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <THnSparse.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TArrayF.h>
#include <TEfficiency.h>
#include <TStyle.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <filesystem>  // C++17
#include <vector>

#include "DUNEVisUtils.hpp" 
#include "file_list.hpp"
#include "Utils.hpp"



// --- HARD CODE HERE ----------------
int nfile_to_analyze = 1e4;
int max_nfiles = 1e4; // Maximum number of files to analyze
// const size_t n_opdet = 4; // 480
const size_t n_opdet = 480; // 480
float light_yield = 20000;
float arapuca_pde = 0.03;

double pe_low = 0.;
double pe_up  = 300.;

double fit_low = 11.;
double fit_up  = 100.;

double min_visibility = 1.e-60;

TString visibility_file_name = "dunevis_fdhd_1x2x6_test.root"; // File with the visibility maps
std::string ana_folder_name = "/exp/dune/data/users/fgalizzi/prod_eminus/seed_2000/ana/"; // Folder where the ana files.root
std::string base_ana_file_name = "solar_ana_dune10kt_1x2x6_hist_";
// -----------------------------------

void buildmu(){
  // --- VISIBILITY STUFF -------------------------------------------------------
  TFile* visibility_file = TFile::Open(visibility_file_name, "READ");
  TTree* photoVisMap = (TTree*)visibility_file->Get("photovisAr/photoVisMap");
  TH1D* hgrid[3] = {nullptr};
  hgrid[0] = (TH1D*)visibility_file->Get("photovisAr/hgrid0");
  hgrid[1] = (TH1D*)visibility_file->Get("photovisAr/hgrid1");
  hgrid[2] = (TH1D*)visibility_file->Get("photovisAr/hgrid2");
  
  TTree* tDimensions = (TTree*)visibility_file->Get("photovisAr/tDimensions");
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
  TTreeReader VisMapReader(photoVisMap);
  TTreeReaderArray<double> opDet_visDirect_reader(VisMapReader, "opDet_visDirect");
  TTreeReaderArray<double> opDet_visReflct_reader(VisMapReader, "opDet_visReflct");
  
  std::cout << "Filling maps..." << std::endl;
  int VisMapEntry = 0; int n_VisMapEntries = photoVisMap->GetEntries();
  while (VisMapReader.Next()){
    for (size_t j = 0; j < n_opdet; ++j) {
      opDet_visMapDirect[VisMapEntry][j] = static_cast<float>(opDet_visDirect_reader[j]);
      opDet_visMapReflct[VisMapEntry][j] = static_cast<float>(opDet_visReflct_reader[j]);
    }
    VisMapEntry++;
    if (VisMapEntry % 1000 == 0)
      std::cout << VisMapEntry << " / " << n_VisMapEntries << "\r" << std::flush;
  }
  std::cout << "Done filling maps..." << std::endl;
  
  // --- OUTPUT ----------------------------------------------------------------
  TFile* out_file = TFile::Open("Mu_Expected.root", "RECREATE");
  Int_t ifile = 0;
  Int_t iev = 0;
  Double_t mu[n_opdet] = {0.};
  Double_t reco_pe[n_opdet] = {0.};
  Double_t hit_t[n_opdet] = {0.};
  Bool_t   hit[n_opdet] = {false};
  TTree* out_tree = new TTree("MuRecoTree", "summary tree");
  out_tree->Branch("ifile", &ifile);
  out_tree->Branch("iev", &iev);
  out_tree->Branch("mu", mu, Form("mu[%d]/D", int(n_opdet)));
  out_tree->Branch("reco_pe", reco_pe, Form("reco_pe[%d]/D", int(n_opdet)));
  out_tree->Branch("hit_t", hit_t, Form("hit_t[%d]/D", int(n_opdet)));
  out_tree->Branch("hit", hit, Form("hit[%d]/O", int(n_opdet)));

  // --- HISTOS ----------------------------------------------------------------
  TH1D* h_Expected_Ophit_OpDet = new TH1D("h_Expected_Ophit_OpDet",
                                          Form("%s;%s;%s","h_Expected_Ophit_OpDet","OpDet","OpHit"),
                                          n_opdet, 0., double(n_opdet));
  TH1D* h_Reco_Ophit_OpDet     = new TH1D("h_Reco_Ophit_OpDet",
                                          Form("%s;%s;%s","h_Reco_Ophit_OpDet","OpDet","OpHit"),
                                          n_opdet, 0., double(n_opdet));
  h_Reco_Ophit_OpDet->SetLineColor(kRed);

  TH1D* h_Expected_Pe = new TH1D("h_Expected_Pe",
                                 Form("%s;%s;%s","h_Expected_Pe","Expected #Pe","Counts"),
                                 300, pe_low, pe_up);
  
  TH1D* h_Expected_PeReco = new TH1D("h_Expected_PeReco",
                                     Form("%s;%s;%s","h_Expected_PeReco","Expected #Pe", "Counts"),
                                     300, pe_low, pe_up);

  TH2D* h2_exp_reco = new TH2D("h2_exp_reco",
                               Form("%s;%s;%s", "", "Reco #Pe", "Expected #Pe"),
                               800, pe_low, pe_up, 800, pe_low, pe_up);

  TH2D* h2_exp_reco_large = new TH2D("h2_exp_reco_large",
                               Form("%s;%s;%s", "", "Reco #Pe", "Expected #Pe"),
                               800*5, pe_low, pe_up*5, 800*5, pe_low, pe_up*5);

  TH2D* h2_ExpPe_HitTime = new TH2D("h2_HitPe_HitTime", Form("%s;%s;%s", "", "HitTime [ticks]", "Expected #Pe"),
                                    200, -1.0, 5.5,
                                    200, pe_low, pe_up);

  // TH2D for events where expected_photons > 5 and no detection
  TH2D* hfail_Etrue_OpDet = new TH2D("hfail_Etrue_OpDet",Form("%s;%s;%s","hfail_Etrue_OpDet","OpDet","E_{True} [MeV]"),
                                     n_opdet, 0., double(n_opdet),
                                     90, 10., 30.);

  TH2D* hfail_Xtrue_OpDet = new TH2D("hfail_Xtrue_OpDet",Form("%s;%s;%s","hfail_Xtrue_OpDet","OpDet","X_{True} [cm]"),
                                     n_opdet, 0., double(n_opdet),
                                     90, -0., 400.);

  TH2D* hfail_Ytrue_OpDet = new TH2D("hfail_Ytrue_OpDet",Form("%s;%s;%s","hfail_Ytrue_OpDet","OpDet","Ytrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, -600., 600.);

  TH2D* hfail_Ztrue_OpDet = new TH2D("hfail_Ztrue_OpDet",Form("%s;%s;%s","hfail_Ztrue_OpDet","OpDet","Ztrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, 0., 1400.);
  
  TH2D* h_Ghost = new TH2D("h_Ghost",Form("%s;%s;%s", "Reco_Ghost", "Event", "Reco_Ghost"), 
                           100, 0, 100, 100, 0, 100); //Ghost PE-> Photons we need ignoring.
  
  TH2D* h_Residual = new TH2D("h_Residual",Form("%s;%s;%s", "", "Reco-True/True", "True"),
                              100, -1000, 1000, 100, 0, 1000);  //*100
  
  // --- DELETES ---------------------------------------------------------------
  delete photoVisMap; delete tDimensions;

  // --- LOOP OVER ANA FILES ---------------------------------------------------
  int nfile_analyzed = 0;
  int idx_file = 0;
  std::vector<double> vec_Expected_Ophit_OpDet(n_opdet, 0.); // Vector to store expected OpHits per OpDet
  std::vector<double> vec_Reco_Ophit_OpDet(n_opdet, 0.); // Vector to store reconstructed OpHits per OpDet
  while(nfile_analyzed < nfile_to_analyze && idx_file < max_nfiles){
    // --- ANA STUFF -----------------------------------------------------------
    std::string ana_file_name = std::string(ana_folder_name+base_ana_file_name+idx_file+".root"); 
    ifile = idx_file;
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
    TTreeReaderValue<std::vector<float>> OpHitPes(treeReader, "OpHitPE");
    TTreeReaderValue<std::vector<float>> OpHitChannels(treeReader, "OpHitChannel");
    TTreeReaderValue<std::vector<float>> OpHitTimes(treeReader, "OpHitTime");
  

    // --- LOOP OVER TREE -----------------------------------------------------
    size_t idx_entry = 0;
    while (treeReader.Next()) {
      iev = int(idx_entry);
      idx_entry++;
      double exp_ph;
      double exp_ph_min;
      double vertex_coor[3] = {*x_true, *y_true, *z_true};
      int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);
      if (tpc_index < 0 || tpc_index >= int(cryo_to_tpc.size())) {
        std::cout << "Event " << iev << " has no valid TPC index!" << std::endl;
        std::cout << vertex_coor[0] << " " << vertex_coor[1] << " " << vertex_coor[2] << std::endl;
        continue; // Skip this entry if TPC index is invalid
      }

      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        // Reset variables
        mu[idx_opdet] = 0.; reco_pe[idx_opdet] = 0.; hit_t[idx_opdet] = 0.; hit[idx_opdet] = false;

        double voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet] + opDet_visMapReflct[tpc_index][idx_opdet];

        exp_ph = (*E_true)*light_yield*voxel_vis*arapuca_pde;
        exp_ph_min = (*E_true)*light_yield*min_visibility*arapuca_pde;
        if(exp_ph==0.) exp_ph = exp_ph_min;
        mu[idx_opdet] = exp_ph;
        vec_Expected_Ophit_OpDet[idx_opdet] += exp_ph;
        h_Expected_Pe->Fill(exp_ph);

        // --- IF RECONSTRUCTED ------------------------------------------------
        auto it = std::find((*OpHitChannels).begin(), (*OpHitChannels).end(), float(idx_opdet));
        size_t idx_hit = std::distance((*OpHitChannels).begin(), it);
        if(idx_hit != (*OpHitChannels).size()){
          reco_pe[idx_opdet] = (*OpHitPes)[idx_hit];
          hit_t[idx_opdet] = (*OpHitTimes)[idx_hit];
          hit[idx_opdet] = true;
          h2_exp_reco->Fill((*OpHitPes)[idx_hit], exp_ph);
          h2_exp_reco_large->Fill((*OpHitPes)[idx_hit], exp_ph);
          h2_ExpPe_HitTime->Fill((*OpHitTimes)[idx_hit], exp_ph);
          h_Expected_PeReco->Fill(exp_ph);
          vec_Reco_Ophit_OpDet[idx_opdet] += (*OpHitPes)[idx_hit];
          h_Residual->Fill(((((*OpHitPes)[idx_hit])-exp_ph)/exp_ph)*100, exp_ph); 
          if(exp_ph == exp_ph_min && (*OpHitPes)[idx_hit] > 0.0 ){    //true==0, ophit>0 per opdet??
            h_Ghost->Fill(exp_ph, (*OpHitPes)[idx_hit]); 
          } //Ghost PE 
        } // reconstructed
        // --- IF NOT RECONSTRUCTED --------------------------------------------
        else {
          // h2_exp_reco->Fill(0., exp_ph);
          // h2_exp_reco_large->Fill(0., exp_ph);
          if (exp_ph > 5.){
            hfail_Etrue_OpDet->Fill(idx_opdet, *E_true);
            hfail_Xtrue_OpDet->Fill(idx_opdet, *x_true);
            hfail_Ytrue_OpDet->Fill(idx_opdet, *y_true);
            hfail_Ztrue_OpDet->Fill(idx_opdet, *z_true);
          }
        } // not reconstructed
      } // end loop over opdets
      out_tree->Fill();
    } // end loop over tree
    ana_file->Close();
  } // end loop over ana files
 
  // --- EXTRA PLOTS -----------------------------------------------------------
  TProfile* h_Expected_PeReco_Prof = h2_exp_reco_large->ProfileX();
  TF1* f1 = new TF1("f1", "pol1", fit_low, fit_up);
  h_Expected_PeReco_Prof->Fit("f1", "R");

  TEfficiency* he_Hit_Prob = new TEfficiency(*h_Expected_PeReco,*h_Expected_Pe);
  he_Hit_Prob->SetTitle("Hit Probability;Expected #Pe;Detection Probability");
  he_Hit_Prob->SetName("HitProb_ExpPe");
  
  // Fill these histogram in a single loop to avoid GetBinContent calls
  for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
    h_Expected_Ophit_OpDet->SetBinContent(idx_opdet, vec_Expected_Ophit_OpDet[idx_opdet]);
    h_Reco_Ophit_OpDet->SetBinContent(idx_opdet, vec_Reco_Ophit_OpDet[idx_opdet]);
  }
  


  out_file->cd();
  out_tree->Write();
  h2_exp_reco->Write();
  he_Hit_Prob->Write();
  h2_exp_reco_large->Write();



  
  // --- SAVE -------------------------------------------------------------------
  TFile* out_file_hist = TFile::Open("out_buildmu.root", "RECREATE");
  out_file_hist->cd();
  hfail_Etrue_OpDet->Write();
  hfail_Xtrue_OpDet->Write();
  hfail_Ytrue_OpDet->Write();
  hfail_Ztrue_OpDet->Write();
  h_Expected_Ophit_OpDet->Write();
  h_Reco_Ophit_OpDet->Write();
  h_Expected_Pe->Write();
  he_Hit_Prob->Write();
  h_Expected_PeReco->Write();
  h2_exp_reco->Write();
  h2_exp_reco_large->Write();
  h_Expected_PeReco_Prof->Write();
  h2_ExpPe_HitTime->Write();
  out_file_hist->Close();
  out_file->Close();

  return;
}
