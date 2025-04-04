// ROOT
#include <Rtypes.h>
#include <TProfile.h>
#include <TString.h>
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
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
// GENERAL
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>  // C++17
// UTILS
#include "DUNEVisUtils.hpp" 
#include "Utils.hpp"




// --- HARD CODE HERE ----------------
int nfile_to_analyze = 6000;
int max_nfiles = 8000; // Maximum number of files to analyze
size_t n_opdet = 480; // 480
// size_t n_opdet = 480; // 480
float light_yield = 20000;
float arapuca_pde = 0.03;

double pe_low = 0.;
double pe_up  = 300.;

double fit_low = 11.;
double fit_up  = 100.;

double min_visibility = 1.e-60;

TString visibility_file_name = "./dunevis_fdhd_1x2x6_test.root"; // File with the visibility maps
 std::string ana_folder_name =  "./ana/"; //ana_nocx/ly27k"; // Folder where the ana files.root ly=27000

// The output file name will contain the fiducial cuts
double fiducial_cut = 0.; // Fiducial cut in cm (y=z=fiducial_cut)
double x_cut = 0.; // x cut in cm


std::string debug_folder_name = "./debug_ana/"; 
//std::string debug_folder_name =  "/pnfs/dune/scratch/users/jdelgadg/debug_nocx/";      //nocx/_yield:27k

                               
 std::string base_ana_file_name = "solar_ana_dune10kt_1x2x6_hist_";
 std::string base_debug_file_name = "test_opdet_";
// -----------------------------------


void debug_buildmu_v2(){
  gStyle->SetPalette(kSunset);
 
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
  double vol_min[3] = {tpc_min[0]+x_cut, tpc_min[1]+fiducial_cut, tpc_min[2]+fiducial_cut};
  double vol_max[3] = {tpc_max[0]-x_cut, tpc_max[1]-fiducial_cut, tpc_max[2]-fiducial_cut};
  std::vector<int> cryo_to_tpc = GetCryoToTPCMap(hgrid, tpc_min, tpc_max); // Get the mapping from cryostat voxel to TPC voxel
  std::cout << "cryo_to_tpc size " << cryo_to_tpc.size() << std::endl;

  double opDet_visDirect[n_opdet];
  // double opDet_visReflct[n_opdet];
  const size_t n_entriesmap = photoVisMap->GetEntries();
  std::vector<std::vector<float>> opDet_visMapDirect(n_entriesmap, std::vector<float>(n_opdet, 0.)); // Map to store visibility for each opdet
  // std::vector<std::vector<float>> opDet_visMapReflct(n_entriesmap, std::vector<float>(n_opdet, 0.)); // Map to store visibility for each opdet
  TTreeReader VisMapReader(photoVisMap);
  TTreeReaderArray<double> opDet_visDirect_reader(VisMapReader, "opDet_visDirect");
  // TTreeReaderArray<double> opDet_visReflct_reader(VisMapReader, "opDet_visReflct");
  
  std::cout << "Filling maps..." << std::endl;
  int VisMapEntry = 0; int n_VisMapEntries = photoVisMap->GetEntries();
  while (VisMapReader.Next()){
    for (size_t j = 0; j < n_opdet; ++j) {
      opDet_visMapDirect[VisMapEntry][j] = static_cast<float>(opDet_visDirect_reader[j]);
      // opDet_visMapReflct[VisMapEntry][j] = static_cast<float>(opDet_visReflct_reader[j]);
    }
    VisMapEntry++;
    if (VisMapEntry % 1000 == 0)
      std::cout << VisMapEntry << " / " << n_VisMapEntries << "\r" << std::flush;
  }
  std::cout << "Done filling maps..." << std::endl;
  

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
  
  TH1D* h_PeReco_Expected= new TH1D("h_PeReco_Expected",
                                     Form("%s;%s;%s","h_PeReco_Expected","Expected #Pe", "Counts"),
                                     300, pe_low, pe_up);

  TH2D* h2_reco_exp = new TH2D("h2_reco_exp",
                               Form("%s;%s;%s", "Reco vs Expected","Expected #Pe","Reco #Pe"),
                               800, pe_low, pe_up, 800, pe_low, pe_up);

  TH2D* h2_true_exp = new TH2D("h2_true_exp",
                               Form("%s;%s;%s", "True vs Expected","Expected #Pe","True #Pe"),
                               800, pe_low, pe_up, 800, pe_low, pe_up);

  TH2D* h2_reco_exp_large = new TH2D("h2_reco_exp_large",
                               Form("%s;%s;%s", "Reco vs Expected", "Expected #Pe","Reco #Pe"),
                               800*3, pe_low, pe_up*3, 800*3, pe_low, pe_up*3);

  TH2D* h2_reco_true = new TH2D("h2_reco_true",
                               Form("%s;%s;%s", "Reco vs True","True #Pe","Reco #PE"),
                               800, pe_low, pe_up, 800, pe_low, pe_up);
//------------------------------------------------------------------------------------------------------------------------------------

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
  int nfile_analyzed = 0;
  int idx_file = 0;
  while(nfile_analyzed < nfile_to_analyze && idx_file < max_nfiles){
    // --- ANA STUFF -----------------------------------------------------------
    std::string ana_file_name = std::string(ana_folder_name+base_ana_file_name+idx_file+".root"); 
    std::string debug_file_name = std::string(debug_folder_name+base_debug_file_name+idx_file+".root");
    // Check whether both files exist
    idx_file++;
    if(!std::filesystem::exists(ana_file_name) || !std::filesystem::exists(debug_file_name)){
      continue;
    }
    nfile_analyzed++;
    std::cout <<nfile_analyzed<<"--"<< ana_file_name << std::endl;

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
    TTreeReaderValue<std::vector<float>> OpHitTimes(treeReader, "OpHitTime");
  
    // Debug stuff
    TFile* debug_file = TFile::Open(debug_file_name.c_str(), "READ");
    TDirectory* dir_debug = (TDirectory*)debug_file->Get("simphcount");
    TTree* t_DetectedPhotons = (TTree*)(dir_debug->Get("DetectedPhotons"));
    int EventID, OpChannel;
    t_DetectedPhotons->SetBranchAddress("EventID", &EventID);
    t_DetectedPhotons->SetBranchAddress("OpChannel", &OpChannel);


    // --- LOOP OVER TREE -----------------------------------------------------
    Long64_t nEntries_debug = t_DetectedPhotons->GetEntries();
    Long64_t entry_debug = 0;
    int idx_event = 0;
    while (treeReader.Next()) {
      double exp_ph;
      double exp_ph_min;
      double vertex_coor[3] = {*x_true, *y_true, *z_true};
      if(!isInFiducialVolume(vertex_coor, vol_min, vol_max, x_cut)){
        idx_event++;
        continue;
      }
      int tpc_index = GetTPCIndex(vertex_coor, hgrid, cryo_to_tpc);
      if (tpc_index < 0 || tpc_index >= int(cryo_to_tpc.size())) {
        idx_event++;
        std::cout << "Event " << *event_true << " has no valid TPC index!" << std::endl;
        std::cout << vertex_coor[0] << " " << vertex_coor[1] << " " << vertex_coor[2] << std::endl;
        continue; // Skip this entry if TPC index is invalid
      }

      t_DetectedPhotons->GetEntry(0);
      std::map<int, int> opdet_detphotons;
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        opdet_detphotons[idx_opdet] = 0;
      }
      while(EventID <= idx_event+1 && entry_debug < nEntries_debug){
        if(EventID == idx_event+1 && OpChannel < n_opdet){
          opdet_detphotons[OpChannel]++;
        }
        entry_debug++;
        t_DetectedPhotons->GetEntry(entry_debug);
      }


      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        double voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];
        // if (voxel_vis >= .5) continue; // Skip if visibility is 1

        exp_ph = (*E_true)*light_yield*voxel_vis*arapuca_pde;
        exp_ph_min = (*E_true)*light_yield*min_visibility*arapuca_pde;
        if(exp_ph==0.) exp_ph = exp_ph_min;
        h_Expected_Ophit_OpDet->SetBinContent(idx_opdet, h_Expected_Ophit_OpDet->GetBinContent(idx_opdet)+exp_ph);
        h_Expected_Pe->Fill(exp_ph);
        if (opdet_detphotons[idx_opdet]>0) h2_true_exp->Fill(exp_ph,opdet_detphotons[idx_opdet]);
        else h2_true_exp->Fill(exp_ph,0.0);

        // --- IF RECONSTRUCTED ------------------------------------------------
        auto it = std::find((*OpHitChannels).begin(), (*OpHitChannels).end(), float(idx_opdet));
        size_t idx_hit = std::distance((*OpHitChannels).begin(), it);
        if(idx_hit != (*OpHitChannels).size()){ 
          h2_reco_exp->Fill(exp_ph,(*OpHitPes)[idx_hit]);
          h2_reco_exp_large->Fill(exp_ph,(*OpHitPes)[idx_hit]);
          h2_ExpPe_HitTime->Fill((*OpHitTimes)[idx_hit], exp_ph);
          h_PeReco_Expected->Fill(exp_ph);
          h_Reco_Ophit_OpDet->SetBinContent(idx_opdet, h_Reco_Ophit_OpDet->GetBinContent(idx_opdet)+exp_ph);
          h_Residual->Fill(((((*OpHitPes)[idx_hit])-exp_ph)/exp_ph)*100, exp_ph); 
          //-------------------------------------------------------------------------
          if (opdet_detphotons[idx_opdet]>0) h2_reco_true->Fill(opdet_detphotons[idx_opdet],(*OpHitPes)[idx_hit]);
          //-----------------------------------------------------------------------------------------------------
          if(exp_ph == exp_ph_min && (*OpHitPes)[idx_hit] > 0.0 ){    //true==0, ophit>0 per opdet??
            h_Ghost->Fill(idx_opdet, (*OpHitPes)[idx_hit]); 
          } //Ghost PE 
        } // reconstructed
        // --- IF NOT RECONSTRUCTED --------------------------------------------
        else {
          h2_reco_exp->Fill(exp_ph,0.);
          h2_reco_exp_large->Fill(exp_ph,0.);
          if (exp_ph > 15.){
            hfail_Etrue_OpDet->Fill(idx_opdet, *E_true);
            hfail_Xtrue_OpDet->Fill(idx_opdet, *x_true);
            hfail_Ytrue_OpDet->Fill(idx_opdet, *y_true);
            hfail_Ztrue_OpDet->Fill(idx_opdet, *z_true);
          }
        } // not reconstructed
      } // end loop over opdets
      idx_event++;
    } // end loop over tree
    ana_file->Close();
    debug_file->Close();
  } // end loop over ana files


  // --- EXTRA PLOTS -----------------------------------------------------------
  gStyle->SetOptFit(1111);
  // auto c3 = new TCanvas();
  TProfile* h_PeReco_Expected_Prof = h2_reco_exp_large->ProfileX();
  h2_reco_exp_large->Delete();
  TF1* f1 = new TF1("f1", "pol1", fit_low, fit_up);
  h_PeReco_Expected_Prof->Fit("f1", "R");

  // auto c4 = new TCanvas();
  TProfile* h_PeTrue_Expected_Prof = h2_true_exp->ProfileX();
  TF1* f2 = new TF1("f2", "pol1", fit_low, fit_up);
  h_PeTrue_Expected_Prof->Fit("f2", "R");
  //------------------------------------------------------------------------------------------------------------
   gStyle->SetOptFit(1111);
   // auto c5 = new TCanvas();
   TProfile* h_reco_true_Prof = h2_reco_true->ProfileX();
   TF1* f3 = new TF1("f3", "pol1", fit_low, fit_up);
   h_reco_true_Prof->Fit("f3", "R");
  //-----------------------------------------------------------------------------------------------------------------

  TEfficiency* he_Hit_Prob = new TEfficiency(*h_PeReco_Expected,*h_Expected_Pe);
  he_Hit_Prob->SetTitle("Hit Probability;Expected #Pe;Detection Probability");
 
   //----Plotting-----------------------------------------------------------------------------------------------------
   // TCanvas* Reco_Expected = new TCanvas("Reco_Expected","Reco_Expected",0,0,800,600);
   // h2_reco_exp->Draw();
   //
   // TCanvas* True_expected = new TCanvas("True_expected","True_expected",0,0,800,600);
   // h2_true_exp->Draw();
   //
   // TCanvas* reco_true = new TCanvas("reco_true","reco_true",0,0,800,600);
   // h2_reco_true->Draw();


  // --- SAVE -------------------------------------------------------------------------------------------------------
  std::string out_file_name = "debug_exp_true_reco_FullTPC_norefl.root";
  if (x_cut > 0 || fiducial_cut > 0.0) {
    out_file_name = Form("debug_exp_true_reco_x%icm_yz%icm_fiducial_norefl.root", int(x_cut), int(fiducial_cut));
  }
  TFile* out_file = TFile::Open(out_file_name.c_str(), "RECREATE");

  gStyle->SetOptFit(1111);
  hfail_Etrue_OpDet->Write();
  hfail_Xtrue_OpDet->Write();
  hfail_Ytrue_OpDet->Write();
  hfail_Ztrue_OpDet->Write();
  h_Expected_Ophit_OpDet->Write();
  h_Reco_Ophit_OpDet->Write();
  h_Expected_Pe->Write();
  he_Hit_Prob->Write();
  h_PeReco_Expected->Write();
  h2_reco_exp->Write();
  h2_true_exp->Write();
   //-----------------------------------------------------------------------------------------------------------------
  h2_reco_true->Write();
  //-----------------------------------------------------------------------------------------------------------------
  h_PeReco_Expected_Prof->Write();
  gStyle->SetOptFit(1111);
  h_PeTrue_Expected_Prof->Write();
  gStyle->SetOptFit(1111);
  h_reco_true_Prof->Write();
  h2_ExpPe_HitTime->Write();
  out_file->Close();

  return;
}


