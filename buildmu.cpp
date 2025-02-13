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

#include "DUNEVisUtils.hpp" 
#include "file_list.hpp"


// --- HARD CODE HERE ----------------
// size_t n_opdet = 8; // 480
size_t n_opdet = 480; // 480
float light_yield = 27000;
float arapuca_pde = 0.03;

double min_visibility = 1.e-60;
double hit_threshold = 1.5; // Will integrate Poisson [0, hit_threshold]

TString visibility_file_name = "./dunevis_fdhd_1x2x6_test_photovisAr.root";
std::string ana_folder_name = "ana/"; // Folder where the ana files.root
// -----------------------------------

void buildmu(){
  gStyle->SetPalette(kSunset);
 
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

  TH1D* h_Expected_Pe = new TH1D("h_Expected_Pe",
                                 Form("%s;%s;%s","h_Expected_Pe","PE","Detection Probability"),
                                 300, 0., 200.);
  
  TH1D* h_Expected_PeReco = new TH1D("h_Expected_PeReco",
                                     Form("%s;%s;%s","h_Expected_PeReco","x","y"),
                                     300, 0., 200.);

  TH2D* h2_exp_reco = new TH2D("h2_exp_reco", Form("%s;%s;%s", "", "Reco #Pe", "Expected #Pe"), 800, 0., 200., 800, 0, 200);

  TH2D* h2_ExpPe_HitTime = new TH2D("h2_HitPe_HitTime", Form("%s;%s;%s", "", "HitTime [ticks]", "Expected #Pe"),
                                    200, -1.5, 4.5,
                                    200, 0, 200.);

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
  
  TH2D* h_Ghost = new TH2D("h_Ghost",Form("%s;%s;%s", "Reco_Ghost", "Event", "Reco_Ghost") 
                           , 100, 0, 100, 100, 0, 100); //Ghost PE-> Photons we need ignoring.
  
  TH2D* h_Residual = new TH2D("h_Residual",Form("%s;%s;%s", "", "Reco-True/True", "True") 
                              , 100, -1000, 1000, 100, 0, 1000);  //*100

  // TEfficiency* he_Ophit_OpDet = nullptr;
  
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
      double exp_ph_min;

      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        int idx_bin[3];
        idx_bin[0] = h3VisMap_opDet[idx_opdet]->GetAxis(0)->FindBin(x_true);
        idx_bin[1] = h3VisMap_opDet[idx_opdet]->GetAxis(1)->FindBin(y_true);
        idx_bin[2] = h3VisMap_opDet[idx_opdet]->GetAxis(2)->FindBin(z_true);
        double voxel_vis = h3VisMap_opDet[idx_opdet]->GetBinContent(idx_bin);
        

        exp_ph = E_true*light_yield*voxel_vis*arapuca_pde;
        exp_ph_min = E_true*light_yield*min_visibility*arapuca_pde;
        if(exp_ph==0.) exp_ph = exp_ph_min;
        h_Expected_Ophit_OpDet->Fill(idx_opdet, exp_ph);
        h_Expected_Pe->Fill(exp_ph);


        auto it = std::find((*OpHitChannels).begin(), (*OpHitChannels).end(), float(idx_opdet));
        size_t idx_hit = std::distance((*OpHitChannels).begin(), it);
        if(idx_hit != (*OpHitChannels).size()){ 
          h2_exp_reco->Fill((*OpHitPes)[idx_hit], exp_ph);
          h2_ExpPe_HitTime->Fill((*OpHitTimes)[idx_hit], exp_ph);
          h_Expected_PeReco->Fill(exp_ph);

          //--------Residual --------------------
          //if (exp_ph > exp_ph_min){
                h_Residual->Fill(((((*OpHitPes)[idx_hit])-exp_ph)/exp_ph)*100, exp_ph); 
             // }
         
          //------Found the Ghost PE-------
           if(exp_ph == exp_ph_min && (*OpHitPes)[idx_hit] > 0.0 ){    //true==0, ophit>0 per opdet??
              h_Ghost->Fill(exp_ph, (*OpHitPes)[idx_hit]); 
              //std::cout << "event_true" << event_true << "--"<<(*OpHitPes)[idx_hit] << std::endl;
             } //Ghost PE
               
        } else {
          h2_exp_reco->Fill(0., exp_ph);
          if (exp_ph > 5.){
            // std::cout << idx_entry << "\t" << idx_opdet << "\t" << E_true << "\t" << x_true << "\t" <<
            //   y_true << "\t" << z_true << "\t" << std::endl;
            hfail_Etrue_OpDet->Fill(idx_opdet, E_true);
            hfail_Xtrue_OpDet->Fill(idx_opdet, x_true);
            hfail_Ytrue_OpDet->Fill(idx_opdet, y_true);
            hfail_Ztrue_OpDet->Fill(idx_opdet, z_true);
          }
        }
      } // end loop over opdets
     
      for(size_t idx_hit=0; idx_hit<(*OpHitChannels).size(); idx_hit++){
        h_Reco_Ophit_OpDet->Fill((*OpHitChannels)[idx_hit], (*OpHitPes)[idx_hit]);
      } 
    } // end loop over tree

    ana_file->Close();
  } // end loop over ana files



  // --- SAVE -------------------------------------------------------------------
  TFile* out_file = TFile::Open("h2_exp_reco.root", "RECREATE");
  h2_exp_reco->Write();
  out_file->Close();

  TFile* out_file2 = TFile::Open("h2_ExpPe_HitTime.root", "RECREATE");
  h2_ExpPe_HitTime->Write();
  out_file2->Close();
  


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
  h2_exp_reco->Draw("colz");
  c_Mu_Pe ->Modified(); c_Mu_Pe ->Update();

  TCanvas* c_HitTime_HitPE = new TCanvas("c_HitTime_HitPE","c_HitTime_HitPE",0,0,800,600);
  c_HitTime_HitPE->cd();
  h2_ExpPe_HitTime->Draw("colz"); 
  c_HitTime_HitPE->Modified(); c_HitTime_HitPE->Update();

  TCanvas* c_Ghost= new TCanvas("c_Ghost","c_Ghost",0,0,800,600);
  c_Ghost->cd();
  h_Ghost->Draw("colz");
  c_Ghost->Modified(); c_Ghost->Update();

  TCanvas* c_Resid= new TCanvas("c_Resid","c_Resid",0,0,800,600);
  c_Resid->cd();
  h_Residual->Draw("colz");
  c_Resid->Modified(); c_Resid->Update();
  
  TEfficiency* he_Hit_Prob = new TEfficiency(*h_Expected_PeReco,*h_Expected_Pe);
  TCanvas* c_Hit_Probability = new TCanvas("c_Hit_Probability","c_Hit_Probability",0,0,800,600);
  c_Hit_Probability->cd();
  he_Hit_Prob->Draw();
  c_Hit_Probability->Modified(); c_Hit_Probability->Update();

  return;
}

