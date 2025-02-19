#include "Rtypes.h"
#include "RtypesCore.h"
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
#include <TArrayF.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <TPRegexp.h>
#include <TObjString.h>

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

double pe_low = 0.;
double pe_up  = 300.;

double fit_low = 11.;
double fit_up  = 100.;

double min_visibility = 1.e-60;
double hit_threshold = 1.5; // Will integrate Poisson [0, hit_threshold]

TString visibility_file_name = "/exp/dune/data/users/dguffant/flash-match/dunevis_fdhd_1x2x6_test_photovisAr_rebin222.root"; 
std::string ana_folder_name = "/exp/dune/data/users/fgalizzi/prod_eminus/ana/";
// -----------------------------------

// --- HELPER FUNCTIONS --------------
std::vector<int> findIndices(const std::vector<float>& vec, const float& target) {
  std::vector<int> indices;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (vec[i] == target) {
      indices.push_back(i);
    }
  }
  return indices;
}

int getRootfileNumber(const std::string& filename) {
  TPRegexp re("_([0-9]+)\\.root$");
  TObjArray* matches = re.MatchS(filename);

  if (matches->GetEntries() > 1) {
    TString numStr = ((TObjString*)matches->At(1))->GetString();
    delete matches;
    return numStr.Atoi();
  }

  delete matches;
  return -1;
}

void buildmu(const TString output_file_path){
  gStyle->SetPalette(kSunset);
 
  // --- VISIBILITY STUFF -------------------------------------------------------
  TFile* visibility_file = TFile::Open(visibility_file_name, "READ");
  // Get the pointer for each opdet
  THnSparseT<TArrayF>* h3VisMap_opDet[n_opdet];
  for(int idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
    TString name_h3_opdet = "h3VisMap_opDet"+std::to_string(idx_opdet)+"_2_2_2";
    visibility_file->GetObject(name_h3_opdet, h3VisMap_opDet[idx_opdet]);
  }

  // --- OUTPUT FILE -----------------------------------------------------------
  TFile* out_file = TFile::Open(output_file_path, "RECREATE");
  Int_t ifile = 0; 
  Int_t iev = 0; 
  Int_t iopdet = 0;
  Double_t mu = 0.0; 
  Double_t reco_pe = 0.0; 
  Double_t hit_t = 0.0; 
  TTree* out_tree = new TTree("t", "summary tree"); 
  out_tree->Branch("ifile", &ifile); 
  out_tree->Branch("ievent", &iev); 
  out_tree->Branch("iopdet", &iopdet); 
  out_tree->Branch("mu", &mu); 
  out_tree->Branch("reco_pe", &reco_pe); 
  out_tree->Branch("hit_t", &hit_t); 
  

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
                               800*3, pe_low, pe_up*3, 800*3, pe_low, pe_up*3);

  TH2D* h2_ExpPe_HitTime = new TH2D("h2_HitPe_HitTime", Form("%s;%s;%s", "", "HitTime [ticks]", "Expected #Pe"),
                                    200, -1.5, 4.5,
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

  // TEfficiency* he_Ophit_OpDet = nullptr;
  
  // --- LOOP OVER ANA FILES ---------------------------------------------------
  for(const auto &file_name : file_list) {
    // --- ANA STUFF -----------------------------------------------------------
    ifile = getRootfileNumber( file_name ); 
    std::cout << file_name << " [" << ifile << "]" << std::endl;
    TFile* ana_file = TFile::Open((ana_folder_name+file_name).c_str(), "READ");
    if (ana_file == nullptr) {
      printf("ERROR: Cannot open file %s/%s\n", 
          ana_folder_name.data(), file_name.data()); 
      ifile++; 
      continue;
    }
    TDirectory* dir = (TDirectory*)ana_file->Get("solarnuana");
    TTree* tree = (TTree*)(dir->Get("MCTruthTree"));

    // Set branches
    float E_true, x_true, y_true, z_true;
    int event_true;
    std::vector<float>* OpHitPes = nullptr;
    std::vector<float>* OpHitChannels = nullptr;
    std::vector<float>* OpHitTimes = nullptr;
    tree->SetBranchAddress("SignalParticleE", &E_true);
    tree->SetBranchAddress("SignalParticleX", &x_true);
    tree->SetBranchAddress("SignalParticleY", &y_true);
    tree->SetBranchAddress("SignalParticleZ", &z_true);
    tree->SetBranchAddress("Event", &event_true);

    tree->SetBranchAddress("OpHitPE", &OpHitPes);
    tree->SetBranchAddress("OpHitChannel", &OpHitChannels);
    tree->SetBranchAddress("OpHitTime", &OpHitTimes);
  
    // --- LOOP OVER TREE -----------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t idx_entry = 0; idx_entry < nEntries; ++idx_entry) {
      iev = idx_entry; 
      tree->GetEntry(idx_entry);
      double exp_ph;
      double exp_ph_min;

      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        iopdet = idx_opdet;
        mu = 0.0; 
        reco_pe = 0.0; 
        hit_t = 0.0; 

        int idx_bin[3];
        idx_bin[0] = h3VisMap_opDet[idx_opdet]->GetAxis(0)->FindBin(x_true);
        idx_bin[1] = h3VisMap_opDet[idx_opdet]->GetAxis(1)->FindBin(y_true);
        idx_bin[2] = h3VisMap_opDet[idx_opdet]->GetAxis(2)->FindBin(z_true);
        double voxel_vis = h3VisMap_opDet[idx_opdet]->GetBinContent(idx_bin);
        
        exp_ph = E_true*light_yield*voxel_vis*arapuca_pde;
        exp_ph_min = E_true*light_yield*min_visibility*arapuca_pde;
        if(exp_ph==0.) exp_ph = exp_ph_min;
          
        mu = exp_ph;

        h_Expected_Ophit_OpDet->SetBinContent(idx_opdet, h_Expected_Ophit_OpDet->GetBinContent(idx_opdet)+exp_ph);
        h_Expected_Pe->Fill(exp_ph);

        // --- IF RECONSTRUCTED ------------------------------------------------
        const auto matches = findIndices( *OpHitChannels, float(idx_opdet) ); 
        if ( matches.empty() ) {
          h2_exp_reco->Fill(0., exp_ph);
          h2_exp_reco_large->Fill(0., exp_ph);
          if (exp_ph > 5.){
            hfail_Etrue_OpDet->Fill(idx_opdet, E_true);
            hfail_Xtrue_OpDet->Fill(idx_opdet, x_true);
            hfail_Ytrue_OpDet->Fill(idx_opdet, y_true);
            hfail_Ztrue_OpDet->Fill(idx_opdet, z_true);
          }
        }
        else {
          for (const auto& idx_hit : matches) {
            reco_pe += (*OpHitPes).at(idx_hit);
            if (idx_hit == matches.front()) {
              hit_t = (*OpHitTimes).at(idx_hit);
            }
          }

          h2_exp_reco->Fill(reco_pe, exp_ph);
          h2_exp_reco_large->Fill(reco_pe, exp_ph);
          h2_ExpPe_HitTime->Fill(hit_t, exp_ph);
          h_Expected_PeReco->Fill(exp_ph);
          h_Reco_Ophit_OpDet->SetBinContent(
              idx_opdet, h_Reco_Ophit_OpDet->GetBinContent(idx_opdet)+exp_ph);
          h_Residual->Fill( (reco_pe/exp_ph)*100, exp_ph); 
          if(exp_ph == exp_ph_min && reco_pe> 0.0 ){    //true==0, ophit>0 per opdet??
            h_Ghost->Fill(exp_ph, reco_pe); 
          } //Ghost PE 
        }

        out_tree->Fill(); 
      } // end loop over opdets
    } // end loop over tree
    ana_file->Close();

    ifile++; 
  } // end loop over ana files


  // --- EXTRA PLOTS -----------------------------------------------------------
  TProfile* h_Expected_PeReco_Prof = h2_exp_reco_large->ProfileX();
  h2_exp_reco_large->Delete();
  TF1* f1 = new TF1("f1", "pol1", fit_low, fit_up);
  h_Expected_PeReco_Prof->Fit("f1", "R");

  TEfficiency* he_Hit_Prob = new TEfficiency(*h_Expected_PeReco,*h_Expected_Pe);
  he_Hit_Prob->SetTitle("Hit Probability;Expected #Pe;Detection Probability");
  

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
  f1->Draw("same");
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
  
  TCanvas* c_Hit_Probability = new TCanvas("c_Hit_Probability","c_Hit_Probability",0,0,800,600);
  c_Hit_Probability->cd();
  he_Hit_Prob->Draw();
  c_Hit_Probability->Modified(); c_Hit_Probability->Update();
  
  // --- SAVE -------------------------------------------------------------------
  out_file->cd(); 
  out_tree->Write(); 
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
  h_Expected_PeReco_Prof->Write();
  h2_ExpPe_HitTime->Write();
  out_file->Close();



  return;
}
