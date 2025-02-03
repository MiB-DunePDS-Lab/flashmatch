#include "Rtypes.h"
#include "THnBase.h"
#include "TString.h"
#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <THnSparse.h>
#include <TArrayF.h>
#include <TEfficiency.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "file_list.hpp"


// --- HARD CODE HERE ----------------
// size_t n_opdet = 30; // 480
size_t n_opdet = 480; // 480
float light_yield = 27000;
float arapuca_pde = 0.02;

double min_visibility = 1.e-60;
double hit_threshold = 1.5; // Will integrate Poisson [0, hit_threshold]

TString visibility_file_name = "./dunevis_fdhd_1x2x6_test_photovisAr.root";
std::string ana_folder_name = "ana/"; // Folder where the ana files.root
// -----------------------------------

void buildmu(){
  
  // --- VISIBILITY STUFF -------------------------------------------------------
  TFile* visibility_file = TFile::Open(visibility_file_name, "READ");
  // Get the pointer for each opdet
  THnSparseT<TArrayF>* h3VisMap_opDet[n_opdet];
  for(int idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
    TString name_h3_opdet = "h3VisMap_opDet"+std::to_string(idx_opdet);
    visibility_file->GetObject(name_h3_opdet, h3VisMap_opDet[idx_opdet]);
    // h3VisMap_opDet[idx_opdet] = rebin_visibility_map(h3VisMap_opDet[idx_opdet], 5, 5, 5);
  }

  

  // --- HISTOS ----------------------------------------------------------------
  TH1D* h_Expected_Ophit_OpDet = new TH1D("h_Expected_Ophit_OpDet",
                                          Form("%s;%s;%s","h_Expected_Ophit_OpDet","OpDet","OpHit"),
                                          n_opdet, 0., double(n_opdet));
  TH1D* h_Reco_Ophit_OpDet = new TH1D("h_Reco_Ophit_OpDet",
                                          Form("%s;%s;%s","h_Reco_Ophit_OpDet","OpDet","OpHit"),
                                          n_opdet, 0., double(n_opdet));
  h_Reco_Ophit_OpDet->SetLineColor(kRed);

  TH2D* h2_Mu_Pe = new TH2D("h2_Mu_Pe", Form("%s;%s;%s", "", "Pe", "Mu"), 500, 0., 100., 500, 0, 100);

  TH2D* h2_HitTime_HitPe = new TH2D("h2_HitTime_HitPe", Form("%s;%s;%s", "", "HitTime", "HitPe"),
                                    200, -1.5, 1.5,
                                    200, 0, 100.);

  // TH2D for events where expected_photons > 5 and no detection
  TH2D* hfail_Etrue_OpDet = new TH2D("hfail_Etrue_OpDet",Form("%s;%s;%s","hfail_Etrue_OpDet","OpDet","Etrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, 0., 20.);

  TH2D* hfail_Xtrue_OpDet = new TH2D("hfail_Xtrue_OpDet",Form("%s;%s;%s","hfail_Xtrue_OpDet","OpDet","Xtrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, -400., 400.);

  TH2D* hfail_Ytrue_OpDet = new TH2D("hfail_Ytrue_OpDet",Form("%s;%s;%s","hfail_Ytrue_OpDet","OpDet","Ytrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, -600., 600.);

  TH2D* hfail_Ztrue_OpDet = new TH2D("hfail_Ztrue_OpDet",Form("%s;%s;%s","hfail_Ztrue_OpDet","OpDet","Ztrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, 0., 1400.);

  // TEfficiency* he_Ophit_OpDet = nullptr;
  
  // --- PDF -------------------------------------------------------------------
  TF1* pdf_poisson = new TF1("pdf_poisson", "TMath::Poisson(x, [0])");
  pdf_poisson->SetNpx(10000);


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
      double exp_ph;
      double Fq = 0.;

      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        int idx_bin[3];
        idx_bin[0] = h3VisMap_opDet[idx_opdet]->GetAxis(0)->FindBin(x_true);
        idx_bin[1] = h3VisMap_opDet[idx_opdet]->GetAxis(1)->FindBin(y_true);
        idx_bin[2] = h3VisMap_opDet[idx_opdet]->GetAxis(2)->FindBin(z_true);
        double voxel_vis = h3VisMap_opDet[idx_opdet]->GetBinContent(idx_bin);
        

        exp_ph = E_true*light_yield*voxel_vis*arapuca_pde;
        if(exp_ph==0.) exp_ph = E_true*light_yield*min_visibility*arapuca_pde;
        h_Expected_Ophit_OpDet->Fill(idx_opdet, exp_ph);


        // Compute a single term of the Fq summatory
        double term = 0.;
        
        pdf_poisson->SetParameter(0, exp_ph);
        double P_hit_mu = 1. - pdf_poisson->Integral(0., hit_threshold, 1e-5);

        auto it = std::find((*OpHitChannels).begin(), (*OpHitChannels).end(), float(idx_opdet));
        size_t idx_hit = std::distance((*OpHitChannels).begin(), it);
        if(idx_hit != (*OpHitChannels).size()){ 
          term = P_hit_mu*pdf_poisson->Eval((*OpHitPes)[idx_hit]);
            // if(x_true<200) {std::cout << idx_entry << "\t" << idx_opdet << "\t" << E_true << "\t" << x_true << "\t" <<
            //   y_true << "\t" << z_true << "\t" << std::endl;
            // }
          h2_Mu_Pe->Fill((*OpHitPes)[idx_hit], exp_ph);
          h2_HitTime_HitPe->Fill((*OpHitTimes)[idx_hit], (*OpHitPes)[idx_hit]);
        } else {
          term = 1. - P_hit_mu;
          h2_Mu_Pe->Fill(0., exp_ph);
          if (exp_ph > 5.){
            // std::cout << idx_entry << "\t" << idx_opdet << "\t" << E_true << "\t" << x_true << "\t" <<
            //   y_true << "\t" << z_true << "\t" << std::endl;
            hfail_Etrue_OpDet->Fill(idx_opdet, E_true);
            hfail_Xtrue_OpDet->Fill(idx_opdet, x_true);
            hfail_Ytrue_OpDet->Fill(idx_opdet, y_true);
            hfail_Ztrue_OpDet->Fill(idx_opdet, z_true);
          }
        }
        Fq += -log(term);


        // if(abs((*OpHitPes)[idx_hit]-exp_ph)>10 && (*OpHitPes)[idx_hit]<1.){
        //   std::cout << idx_opdet << std::endl; 
        // }

      } // end loop over opdets
     
      for(size_t idx_hit=0; idx_hit<(*OpHitChannels).size(); idx_hit++){
        h_Reco_Ophit_OpDet->Fill((*OpHitChannels)[idx_hit], (*OpHitPes)[idx_hit]);
      } 
    } // end loop over tree

    ana_file->Close();
  } // end loop over ana files


  // for each bin in the x-axis of h2_Mu_Pe project on the y-axis and normalize
  // the projection to have area=1
  for(int idx_x=1; idx_x<=h2_Mu_Pe->GetNbinsX(); idx_x++){
    TH1D* h1_proj = h2_Mu_Pe->ProjectionY("h1_proj", idx_x, idx_x);
    h1_proj->Scale(1./h1_proj->Integral());
    for(int idx_y=1; idx_y<=h2_Mu_Pe->GetNbinsY(); idx_y++){
      h2_Mu_Pe->SetBinContent(idx_x, idx_y, h1_proj->GetBinContent(idx_y));
    }
  }

  


  // --- PLOTTING ----------------------------------------------------------------
  TCanvas* cfail_Etrue_Opdet = new TCanvas("cfail_Etrue_Opdet","cfail_Etrue_Opdet",0,0,800,600);
  cfail_Etrue_Opdet->cd();
  hfail_Etrue_OpDet->Draw("colz");
  cfail_Etrue_Opdet->Modified(); cfail_Etrue_Opdet->Update();

  TCanvas* cfail_Xtrue_Opdet = new TCanvas("cfail_Xtrue_Opdet","cfail_Xtrue_Opdet",0,0,800,600);
  cfail_Xtrue_Opdet->cd();
  hfail_Xtrue_OpDet->Draw("colz");
  cfail_Xtrue_Opdet->Modified(); cfail_Xtrue_Opdet->Update();

  TCanvas* cfail_Ytrue_Opdet = new TCanvas("cfail_Ytrue_Opdet","cfail_Ytrue_Opdet",0,0,800,600);
  cfail_Ytrue_Opdet->cd();
  hfail_Ytrue_OpDet->Draw("colz");
  cfail_Ytrue_Opdet->Modified(); cfail_Ytrue_Opdet->Update();

  TCanvas* cfail_Ztrue_Opdet = new TCanvas("cfail_Ztrue_Opdet","cfail_Ztrue_Opdet",0,0,800,600);
  cfail_Ztrue_Opdet->cd();
  hfail_Ztrue_OpDet->Draw("colz");
  cfail_Ztrue_Opdet->Modified(); cfail_Ztrue_Opdet->Update();

  TCanvas* c_Ophit_OpDet= new TCanvas("c_Ophit_OpDet","c_Ophit_OpDet",0,0,800,600);
  c_Ophit_OpDet->cd();
  h_Expected_Ophit_OpDet->Draw();
  h_Reco_Ophit_OpDet->Draw("same");
  c_Ophit_OpDet->Modified(); c_Ophit_OpDet->Update();

  TCanvas* c_Mu_Pe = new TCanvas("c_Mu_Pe ","c_Mu_Pe ",0,0,800,600);
  c_Mu_Pe ->cd();
  h2_Mu_Pe->Draw("colz");
  c_Mu_Pe ->Modified(); c_Mu_Pe ->Update();

  TCanvas* c_HitTime_HitPE = new TCanvas("c_HitTime_HitPE","c_HitTime_HitPE",0,0,800,600);
  c_HitTime_HitPE->cd();
  h2_HitTime_HitPe->Draw("colz"); 
  c_HitTime_HitPE->Modified(); c_HitTime_HitPE->Update();

  // TCanvas* c_Eff_ExpReco = new TCanvas("c_Eff_ExpReco","c_Eff_ExpReco",0,0,800,600);
  // c_Eff_ExpReco->cd();
  // he_Ophit_OpDet = new TEfficiency(*h_Reco_Ophit_OpDet,*h_Expected_Ophit_OpDet); 
  // c_Eff_ExpReco->Modified(); c_Eff_ExpReco->Update();
  // --- FILE CLOSURE ---------------
  // ana_file->Close(); visibility_file->Close();
  // delete ana_file;   delete visibility_file;
  return;
}

