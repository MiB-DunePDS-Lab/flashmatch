
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

#include <THnSparse.h>
#include <TTree.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>

#include "file_list.hpp"



// --- HARD CODE HERE ----------------
int rebinning = 4;

size_t n_opdet = 480; // 480
float light_yield = 27000;
float arapuca_pde = 0.02;

double min_visibility = 1.e-60;
double hit_threshold = 1.5; // Will integrate Poisson [0, hit_threshold]

TString exp_reco_file_name = "exp_reco.root";
TString visibility_file_name = "./dunevis_fdhd_1x2x6_test_photovisAr.root";
std::string ana_folder_name = "ana/"; // Folder where the ana files.root
// -----------------------------------

void compute_likelihood(){

  // --- INPUT ----------------------------------------------------------------
  TFile* exp_reco_file = TFile::Open(exp_reco_file_name, "READ");
  TH2D* h2_exp_reco = nullptr;
  exp_reco_file->GetObject("h2_Mu_Pe", h2_exp_reco);
  // exp_reco_file->GetObject("h2_exp_reco", h2_exp_reco);

  h2_exp_reco->Rebin2D(rebinning, rebinning);

  // for each bin in the y-axis of h2_exp_reco project on the x-axis and normalize
  // the projection to have area=1
  for(int idx_y=1; idx_y<=h2_exp_reco->GetNbinsY(); idx_y++){
    TH1D* h1_proj = h2_exp_reco->ProjectionX("h1_proj", idx_y, idx_y);
    h1_proj->Scale(1./h1_proj->Integral());
    for(int idx_x=1; idx_x<=h2_exp_reco->GetNbinsX(); idx_x++){
      h2_exp_reco->SetBinContent(idx_y, idx_x, h1_proj->GetBinContent(idx_x));
    }
  }


  // --- VISIBILITY STUFF -------------------------------------------------------
  TFile* visibility_file = TFile::Open(visibility_file_name, "READ");
  // Get the pointer for each opdet
  THnSparseT<TArrayF>* h3VisMap_opDet[n_opdet];
  for(int idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
    TString name_h3_opdet = "h3VisMap_opDet"+std::to_string(idx_opdet);
    visibility_file->GetObject(name_h3_opdet, h3VisMap_opDet[idx_opdet]);
    // h3VisMap_opDet[idx_opdet] = rebin_visibility_map(h3VisMap_opDet[idx_opdet], 5, 5, 5);
  }
  
  // --- PDF -------------------------------------------------------------------
  TF1* pdf_poisson = new TF1("pdf_poisson", "TMath::Poisson(x, [0])");
  pdf_poisson->SetNpx(10000);


  size_t fake_counter = 0;
  size_t evt_counter = 0;

  // --- LOOP OVER ANA FILES ---------------------------------------------------
  for(const auto &file_name : file_list){
    // --- ANA STUFF -----------------------------------------------------------
    std::cout << file_name << std::endl;
    TFile* ana_file = TFile::Open((ana_folder_name+file_name).c_str(), "READ");
    TDirectory* dir = (TDirectory*)ana_file->Get("solarnuana");
    TTree* tree = (TTree*)(dir->Get("MCTruthTree"));

    // Set branches
    float E_true, x_true, y_true, z_true;
    std::vector<float>* OpHitPes = nullptr;
    std::vector<float>* OpHitChannels = nullptr;
    std::vector<float>* OpHitTimes = nullptr;
    tree->SetBranchAddress("SignalParticleE", &E_true);
    tree->SetBranchAddress("SignalParticleX", &x_true);
    tree->SetBranchAddress("SignalParticleY", &y_true);
    tree->SetBranchAddress("SignalParticleZ", &z_true);

    tree->SetBranchAddress("OpHitPE", &OpHitPes);
    tree->SetBranchAddress("OpHitChannel", &OpHitChannels);
    tree->SetBranchAddress("OpHitTime", &OpHitTimes);
  
    // --- LOOP OVER TREE -----------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t idx_entry = 0; idx_entry < nEntries; ++idx_entry) {
      tree->GetEntry(idx_entry);
      double exp_ph, exp_ph_fake;
      double Fq = 0.;
      double Fq_fake = 0.;


      float x_fake, y_fake, z_fake;
      x_fake = 100.;
      y_fake = 100.;
      z_fake = 100.;

      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        int idx_bin[3];
        idx_bin[0] = h3VisMap_opDet[idx_opdet]->GetAxis(0)->FindBin(x_true);
        idx_bin[1] = h3VisMap_opDet[idx_opdet]->GetAxis(1)->FindBin(y_true);
        idx_bin[2] = h3VisMap_opDet[idx_opdet]->GetAxis(2)->FindBin(z_true);
        double voxel_vis = h3VisMap_opDet[idx_opdet]->GetBinContent(idx_bin);

        idx_bin[0] = h3VisMap_opDet[idx_opdet]->GetAxis(0)->FindBin(x_fake);
        idx_bin[1] = h3VisMap_opDet[idx_opdet]->GetAxis(1)->FindBin(y_fake);
        idx_bin[2] = h3VisMap_opDet[idx_opdet]->GetAxis(2)->FindBin(z_fake);
        double voxel_vis_fake = h3VisMap_opDet[idx_opdet]->GetBinContent(idx_bin);
        

        exp_ph = E_true*light_yield*voxel_vis*arapuca_pde;
        exp_ph_fake = E_true*light_yield*voxel_vis_fake*arapuca_pde;
        if(exp_ph==0.) exp_ph = E_true*light_yield*min_visibility*arapuca_pde;
        if(exp_ph_fake==0.) exp_ph_fake = E_true*light_yield*min_visibility*arapuca_pde;

        // Compute a single term of the Fq summatory
        double term = 0.;
        double term_fake = 0.;
        pdf_poisson->SetParameter(0, exp_ph);
        double P_hit_mu = 1. - pdf_poisson->Integral(0., hit_threshold, 1e-5);

        pdf_poisson->SetParameter(0, exp_ph_fake);
        double P_hit_mu_fake = 1. - pdf_poisson->Integral(0., hit_threshold, 1e-5);

        auto it = std::find((*OpHitChannels).begin(), (*OpHitChannels).end(), float(idx_opdet));
        size_t idx_hit = std::distance((*OpHitChannels).begin(), it);
 
        

        if(idx_hit != (*OpHitChannels).size()){ 
          term = P_hit_mu*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin((*OpHitPes)[idx_hit], exp_ph));
          term_fake = P_hit_mu_fake*h2_exp_reco->GetBinContent(h2_exp_reco->FindBin((*OpHitPes)[idx_hit], exp_ph));
        } else {
          term = 1. - P_hit_mu;
          term_fake = 1. - P_hit_mu_fake;
        }
        Fq += -log(term);
        Fq_fake += -log(term_fake);


      } // end loop over opdets
      if(Fq_fake < Fq){
        fake_counter++;
      }
      evt_counter++;


    } // end loop over tree
    ana_file->Close();

    std::cout << "Fake counter: " << fake_counter << "/" << evt_counter << std::endl;
  } // end loop over ana files










}
