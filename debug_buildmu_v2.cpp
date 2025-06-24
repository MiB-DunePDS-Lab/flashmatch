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
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>
#include <filesystem>  // C++17
// UTILS
#include "DUNEVisUtils.hpp" 
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTreeReaderValue.h"
#include "Utils.hpp"

TGraphErrors* th2d_to_tgraph_mpv(TH2D* h2, const std::string& name){
  TGraphErrors* g = new TGraphErrors();
  g->SetName(name.c_str()); g->SetTitle(h2->GetTitle());
  for(int idx_x=0; idx_x<=h2->GetNbinsX()+1; idx_x++){
    TH1D* h1_proj = h2->ProjectionY("h1_proj", idx_x, idx_x);
    if(h1_proj->GetEntries() > 0){
      double mpv = h1_proj->GetXaxis()->GetBinCenter(h1_proj->GetMaximumBin());
      double err = h1_proj->GetBinWidth(0); // CHECK THIS
      g->SetPoint(idx_x, h2->GetXaxis()->GetBinCenter(idx_x), mpv);
      g->SetPointError(idx_x, 0., err);
    }
  }
  return g; 
}


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
std::string ana_folder_name =  "./ophitfinder_edep_debug/ana/";


// The output file name will contain the fiducial cuts
double fiducial_cut = 30.; // Fiducial cut in cm (y=z=fiducial_cut)
double x_cut = 50.; // x cut in cm


std::string debug_folder_name = "./ophitfinder_edep_debug/debug/"; 
//std::string debug_folder_name =  "/pnfs/dune/scratch/users/jdelgadg/debug_nocx/";      //nocx/_yield:27k

                               
 std::string base_ana_file_name = "solar_ana_dune10kt_1x2x6_hist_";
 std::string base_debug_file_name = "test_opdet_";
// -----------------------------------


void debug_buildmu_v2(){
  gStyle->SetPalette(kSunset);
  gStyle->SetOptFit(1111);
 
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
  TH1D* h_exp_opdet  = new TH1D("h_exp_opdet",
                                Form("%s;%s;%s","h_exp_opdet","OpDet","OpHit"),
                                n_opdet, 0., double(n_opdet));
  TH1D* h_reco_opdet = new TH1D("h_reco_opdet",
                                Form("%s;%s;%s","h_reco_opdet","OpDet","OpHit"),
                                n_opdet, 0., double(n_opdet));
  h_reco_opdet->SetLineColor(kRed);

  TH1D* h_exp = new TH1D("h_exp",
                         Form("%s;%s;%s","h_exp","Expected #Pe","Counts"),
                         300, pe_low, pe_up);
  
  TH1D* h_expreco= new TH1D("h_expreco",
                            Form("%s;%s;%s","h_expreco","Expected #Pe", "Counts"),
                            300, pe_low, pe_up);

  TH1D* h_true = new TH1D("h_true",
                            Form("%s;%s;%s","h_true","True #Pe", "Counts"),
                            300, pe_low, pe_up);

  TH1D* h_truereco = new TH1D("h_truereco",
                            Form("%s;%s;%s","h_truereco","True #Pe", "Counts"),
                            300, pe_low, pe_up);

  TH2D* h2_edep_etrue = new TH2D("h2_edep_etrue",
                                 Form("%s;%s;%s", "h2_edep_etrue", "E_{true} [MeV]", "E_{dep} [MeV]"),
                                 100, 14, 26, 200, 5, 80);

  TH2D* h2_true_exp = new TH2D("h2_true_exp",
                               Form("%s;%s;%s", "True vs Expected","Expected #Pe","True #Pe"),
                               500, pe_low, pe_up, 500, pe_low, pe_up);
  
  TH2D* h2_true_expedep = new TH2D("h2_true_expedep",
                               Form("%s;%s;%s", "True vs Expected using E deposited","Expected #Pe","True #Pe"),
                               500, pe_low, pe_up, 500, pe_low, pe_up);

  TH2D* h2_reco_true = new TH2D("h2_reco_true",
                               Form("%s;%s;%s", "Reco vs True","True #Pe","Reco #PE"),
                               500, pe_low, pe_up, 500, pe_low, pe_up);

  TH2D* h2_reco_exp = new TH2D("h2_reco_exp",
                               Form("%s;%s;%s", "Reco vs Expected","Expected #Pe","Reco #Pe"),
                               500, pe_low, pe_up, 500, pe_low, pe_up);

  TH2D* h2_reco_exp_large = new TH2D("h2_reco_exp_large",
                               Form("%s;%s;%s", "Reco vs Expected", "Expected #Pe","Reco #Pe"),
                               500*3, pe_low, pe_up*3, 500*3, pe_low, pe_up*3);
//------------------------------------------------------------------------------------------------------------------------------------

  TH2D* h2_exp_hittime = new TH2D("h2_exp_hittime ", Form("%s;%s;%s", "", "HitTime [ticks]", "Expected #Pe"),
                                    200, -1.5, 4.5,
                                    200, pe_low, pe_up);

  // TH2D for events where expected_photons > thr and no detection
  TH2D* h2_fail_Etrue_opdet = new TH2D("h2_fail_Etrue_opdet",Form("%s;%s;%s","h2_fail_Etrue_opdet","OpDet","E_{True} [MeV]"),
                                     n_opdet, 0., double(n_opdet),
                                     90, 10., 30.);

  TH2D* h2_fail_Xtrue_opdet = new TH2D("h2_fail_Xtrue_opdet",Form("%s;%s;%s","h2_fail_Xtrue_opdet","OpDet","X_{True} [cm]"),
                                     n_opdet, 0., double(n_opdet),
                                     90, -0., 400.);

  TH2D* h2_fail_Ytrue_opdet = new TH2D("h2_fail_Ytrue_opdet",Form("%s;%s;%s","h2_fail_Ytrue_opdet","OpDet","Ytrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, -600., 600.);

  TH2D* h2_fail_Ztrue_opdet = new TH2D("h2_fail_Ztrue_opdet",Form("%s;%s;%s","h2_fail_Ztrue_opdet","OpDet","Ztrue"),
                                     n_opdet, 0., double(n_opdet),
                                     90, 0., 1400.);

  TH2D* h2_noreco_Etrue_Xtrue = new TH2D("h2_noreco_Etrue_Xtrue", Form("%s;%s;%s", "h2_noreco_Etrue_Xtrue", "E_{True} [MeV]", "X_{True} [cm]"),
                                        90, 10., 30,
                                        90, 0.,  400.);

  TH2D* h2_noreco_Xtrue_Ytrue = new TH2D("h2_noreco_Xtrue_Ytrue", Form("%s;%s;%s", "h2_noreco_Xtrue_Ytrue", "X_{True} [cm]", "Y_{True} [cm]"),
                                        90, 0., 400.,
                                        90, -600., 600.);

  TH2D* h2_noreco_Ytrue_Ztrue = new TH2D("h2_noreco_Ytrue_Ztrue", Form("%s;%s;%s", "h2_noreco_Ytrue_Ztrue", "Y_{True} [cm]", "Z_{True} [cm]"),
                                        90, -600., 600.,
                                        90, 0., 1400.);
  
  
  TH2D* h_ghost = new TH2D("h_ghost",Form("%s;%s;%s", "Reco_Ghost", "Event", "Reco_Ghost"), 
                           100, 0, 100, 100, 0, 100); //Ghost PE-> Photons we need ignoring.
  
  TH2D* h_residual = new TH2D("h_residual",Form("%s;%s;%s", "", "Reco-True/True", "True"),
                              100, -1000, 1000, 100, 0, 1000);  //*100
  
  TH1D* h_dxtrack = new TH1D("h_dxtrack",
                             Form("%s;%s;%s","h_dxtrack","dX [cm]","Counts"), 
                             50, 0, 100);

  TH1D* h_dxtrack_outliers = new TH1D("h_dxtrack_outliers",
                                      Form("%s;%s;%s","h_dxtrack_outliers","dX [cm]","Counts"),
                                      50, 0, 100);
  
  TH2D* h2_charge_etrue = new TH2D("h2_charge_etrue",
                                   Form("%s;%s;%s", "h2_charge_etrue", "E_{true} [MeV]", "charge"),
                                   50, 15, 25,
                                   200, 0, 60000);
  
  TH2D* h2_chargeperenergytrue_drifttime = new TH2D("h2_chargeperenergytrue_drifttime",
                                                Form("%s;%s;%s", "h2_chargeperenergytrue_drifttime", "Drift Time [ticks]", "Charge/Energy [ADC\times ticks/MeV]"), 
                                                60, 0, 4500, 60, 0, 280);

  TH2D* h2_chargeperenergydep_drifttime = new TH2D("h2_chargeperenergydep_drifttime",
                                                Form("%s;%s;%s", "h2_chargeperenergydep_drifttime", "Drift Time [ticks]", "Charge/Energy [ADC\times ticks/MeV]"), 
                                                60, 0, 4500, 60, 0, 280);

  TH2D* h2_dcperde_dt = new TH2D("h2_dcperde_dt",
                                 Form("%s;%s;%s", "h2_dcperde_dt", "Drift Time [Ticks]", "d(#Electron)/dEnergy_{dep} [e-/MeV]"), 
                                 150, 0, 4500, 100, 0, 33000);

  // --- EXTRA VARIABLES -------------------------------------------------------
  std::vector<std::tuple<double, double, double>> true_calib_info; // Charge, time, true energy
  std::vector<std::tuple<std::vector<int>, std::vector<float>, std::vector<float>, double>> true_finecalib_info; // #Electrons, EDepXs, EDeps, ETrue
  double drift_velocity = 0.0; // Drift velocity in cm/tick, will be set later from the config_tree

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
    TTreeReaderValue<std::vector<float>> EDepList(treeReader, "TSignalEDepList");
    TTreeReaderValue<std::vector<float>> EDepX(treeReader, "TSignalXDepList");
    TTreeReaderValue<std::vector<int>> EDepElectron(treeReader, "TSignalElectronDepList");
    
    TTree* solar_tree = static_cast<TTree*>(ana_file->Get("solarnuana/SolarNuAnaTree"));
    TTreeReader solar_treeReader(solar_tree);
    TTreeReaderValue<float> Charge(solar_treeReader, "Charge");

    TTree* config_tree = static_cast<TTree*>(ana_file->Get("solarnuana/ConfigTree"));
    // TTreeReader configReader(config_tree);
    // TTreeReaderValue<float> DetectorSizeX(configReader, "DetectorSizeX");
    // TTreeReaderValue<float>DetectorDriftTime(configReader, "DetectorDriftTime");
    int DetectorSizeX = 6000; // Hardcoded for now, should be read from config_tree
    int DetectorDriftTime = 2; // Hardcoded for now, should be read from config_tree
    config_tree->SetBranchAddress("DetectorSizeX", &DetectorSizeX);
    config_tree->SetBranchAddress("DetectorDriftTime", &DetectorDriftTime);
    config_tree->GetEntry(0);
    drift_velocity = double(DetectorSizeX) / double(DetectorDriftTime);
  
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
    while (treeReader.Next() && solar_treeReader.Next()) {
      double exp_ph, exp_ph_edep;
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

      double EDep = std::accumulate((*EDepList).begin(), (*EDepList).end(), 0.0);
      // double EDep = *std::max_element((*EDepList).begin(), (*EDepList).end());
      h2_edep_etrue->Fill(*E_true, EDep);
      h2_charge_etrue->Fill(*E_true, (*Charge));
      h2_chargeperenergytrue_drifttime->Fill((*x_true)/drift_velocity, (*Charge)/(*E_true));
      h2_chargeperenergydep_drifttime->Fill((*x_true)/drift_velocity, 3*(*Charge)/(std::accumulate((*EDepList).begin(), (*EDepList).end(), 0.0)));
      h_dxtrack->Fill((*std::max_element((*EDepX).begin(), (*EDepX).end()) - 
                      *std::min_element((*EDepX).begin(), (*EDepX).end())));
      if ((*Charge)/(*E_true) > 140.){
        h_dxtrack_outliers->Fill((*std::max_element((*EDepX).begin(), (*EDepX).end()) - 
                                  *std::min_element((*EDepX).begin(), (*EDepX).end())));
      }
      for (size_t idx_edep=0; idx_edep<EDepList->size(); idx_edep++){
        h2_dcperde_dt->Fill((*EDepX)[idx_edep]/drift_velocity, double((*EDepElectron)[idx_edep])/(*EDepList)[idx_edep]);
      }
      true_calib_info.push_back(std::make_tuple(*Charge, (*x_true)/drift_velocity, *E_true));
      true_finecalib_info.push_back(std::make_tuple(*EDepElectron, *EDepX, *EDepList, *E_true));
      // reco_
      // std::cout << EDep << " " << *E_true << std::endl;

      // --- LOOP OVER OPDETS -------------------------------------------------
      for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        double voxel_vis = opDet_visMapDirect[tpc_index][idx_opdet];// + opDet_visMapReflct[tpc_index][idx_opdet];
        // if (voxel_vis >= .5) continue; // Skip if visibility is 1
        
        exp_ph = (*E_true)*light_yield*voxel_vis*arapuca_pde;
        exp_ph_min = (*E_true)*light_yield*min_visibility*arapuca_pde;
        if(exp_ph==0.) exp_ph = exp_ph_min;

        exp_ph_edep = EDep*light_yield*voxel_vis*arapuca_pde;
        if(exp_ph==0.) exp_ph_edep = exp_ph_min;


        h_exp_opdet->SetBinContent(idx_opdet, h_exp_opdet->GetBinContent(idx_opdet)+exp_ph);
        h_exp->Fill(exp_ph);
        h_true->Fill(opdet_detphotons[idx_opdet]);
        h2_true_exp->Fill(exp_ph,opdet_detphotons[idx_opdet]);
        h2_true_expedep->Fill(exp_ph_edep,opdet_detphotons[idx_opdet]);

        const auto matches = findIndices(*OpHitChannels, float(idx_opdet));
        double reco_pe = 0;
        // --- IF RECONSTRUCTED ------------------------------------------------
        if(!matches.empty()){
          for (const auto& idx_hit : matches){
            reco_pe += (*OpHitPes)[idx_hit];
            h2_exp_hittime->Fill((*OpHitTimes)[idx_hit], exp_ph);
          }
          h2_reco_exp->Fill(exp_ph,reco_pe);
          h2_reco_exp_large->Fill(exp_ph,reco_pe);
          h_expreco->Fill(exp_ph);
          h_truereco->Fill(opdet_detphotons[idx_opdet]);
          h_reco_opdet->SetBinContent(idx_opdet, h_reco_opdet->GetBinContent(idx_opdet)+reco_pe);
          h_residual->Fill((((reco_pe)-exp_ph)/exp_ph)*100, exp_ph); 
          h2_reco_true->Fill(opdet_detphotons[idx_opdet],reco_pe);
          if (reco_pe > 6 && opdet_detphotons[idx_opdet] == 0)
            std::cout << ana_file_name << "\n" << "Ghost pes in evt: " << idx_event << " Ch: " << idx_opdet << " Reco: " << reco_pe << std::endl;
          
          if(exp_ph == exp_ph_min && reco_pe > 0.0 ){    //true==0, ophit>0 per opdet??
            h_ghost->Fill(idx_opdet, reco_pe); 
          } //Ghost PE 
        } // reconstructed
        // --- IF NOT RECONSTRUCTED --------------------------------------------
        else {
          h2_reco_exp->Fill(exp_ph,reco_pe);
          h2_reco_exp_large->Fill(exp_ph,reco_pe);
          h2_reco_true->Fill(opdet_detphotons[idx_opdet],reco_pe);
          if (exp_ph > 15.){
            h2_fail_Etrue_opdet->Fill(idx_opdet, *E_true);
            h2_fail_Xtrue_opdet->Fill(idx_opdet, *x_true);
            h2_fail_Ytrue_opdet->Fill(idx_opdet, *y_true);
            h2_fail_Ztrue_opdet->Fill(idx_opdet, *z_true);
          }
          if (opdet_detphotons[idx_opdet] > 30 )
            std::cout << ana_file_name << "\n" << "Event: " << idx_event << " OpDet: " << idx_opdet << " True: " << opdet_detphotons[idx_opdet] << " Exp: " << exp_ph << std::endl;
        } // not reconstructed
      } // end loop over opdets
      // --- IF NO DETECTION -----------------------------------------------
      if ((*OpHitChannels).size() == 0){
        h2_noreco_Etrue_Xtrue->Fill(*E_true, *x_true);
        h2_noreco_Xtrue_Ytrue->Fill(*x_true, *y_true);
        h2_noreco_Ytrue_Ztrue->Fill(*y_true, *z_true);
        // std::cout << "------------------------" << std::endl;
        // for(size_t idx_opdet=0; idx_opdet<n_opdet; idx_opdet++){
        //   std::cout << " " << opdet_detphotons[idx_opdet] << " ";
        // }
        // std::cout << std::endl;
      }
      idx_event++;
    } // end loop over tree
    ana_file->Close();
    debug_file->Close();
  } // end loop over ana files


  // --- EXTRA PLOTS -----------------------------------------------------------
  TProfile* hp_expreco = h2_reco_exp_large->ProfileX();
  h2_reco_exp_large->Delete();
  TF1* f1 = new TF1("f1", "pol1", fit_low, fit_up);
  hp_expreco->Fit("f1", "R");

  TProfile* hp_true_exp = h2_true_exp->ProfileX();
  TF1* f2 = new TF1("f2", "pol1", fit_low, fit_up);
  hp_true_exp->Fit("f2", "R");

  TProfile* h_reco_true_Prof = h2_reco_true->ProfileX();
  TF1* f3 = new TF1("f3", "pol1", fit_low, fit_up);
  h_reco_true_Prof->Fit("f3", "R");
  
  // TProfile* hp_chargeperenergytrue_drifttime = h2_chargeperenergytrue_drifttime->ProfileX();
  // hp_chargeperenergytrue_drifttime->SetName("hp_chargeperenergytrue_drifttime");
  TGraphErrors* g_chargeperenergytrue_drifttime = th2d_to_tgraph_mpv(h2_chargeperenergytrue_drifttime, "g_chargeperenergytrue_drifttime");
  TF1* f4 = new TF1("f4", "[0]*exp(-x/[1])", 200, 4500);
  f4->SetParNames("A", "#tau");
  f4->SetParameters(140, 20000);
  g_chargeperenergytrue_drifttime->Fit("f4", "R");
  double f4_par0 = f4->GetParameter(0); double f4_par1 = f4->GetParameter(1);
  f4 = new TF1("f4", "[0]*exp(-x/[1])", -4500, 4500);
  f4->SetParNames("A", "#tau"); f4->SetParameters(f4_par0, f4_par1);
  double tau_inv = 1./f4->GetParameter(1);

  auto max_it = std::max_element(true_calib_info.begin(), true_calib_info.end(),
                                     [](const auto& a, const auto& b) {
                                         return std::get<0>(a) < std::get<0>(b);
                                     });
  double max_charge = std::get<0>(*max_it);

  auto min_it = std::min_element(true_calib_info.begin(), true_calib_info.end(),
                                     [](const auto& a, const auto& b) {
                                         return std::get<0>(a) < std::get<0>(b);
                                     });
  double min_charge = std::get<0>(*min_it);

  TH2D* h2_corrcharge_etrue = new TH2D("h2_corrcharge_etrue",
                                       Form("%s;%s;%s", "h2_corrcharge_etrue", "E_{true} [MeV]", "Corrected Charge"),
                                       200, 14, 26, 200, min_charge*0.9, max_charge*1.1);

  for (const auto& [charge, drift_time, e_true] : true_calib_info) {
    double corrected_charge = charge * exp(drift_time * tau_inv);
    h2_corrcharge_etrue->Fill(e_true, corrected_charge);
  }

  TProfile* hp_corrcharge_etrue = h2_corrcharge_etrue->ProfileX();
  hp_corrcharge_etrue->SetName("hp_corrcharge_etrue");
  TF1* f5 = new TF1("f5", "pol1", 14, 26);
  hp_corrcharge_etrue->Fit("f5", "R");

  TGraphErrors* g_dcperde_dt = th2d_to_tgraph_mpv(h2_dcperde_dt, "g_dcperde_dt");
  TF1* f6 = new TF1("f6", "[0]*exp(-x/[1])", 0, 4500);
  f6->SetParNames("A", "#tau");
  f6->SetParameters(33000, 20000);
  g_dcperde_dt->Fit("f6", "R");
  double f6_par0 = f6->GetParameter(0); double f6_par1 = f6->GetParameter(1);
  f6 = new TF1("f6", "[0]*exp(-x/[1])", -4500, 4500);
  f6->SetParNames("A", "#tau"); f6->SetParameters(f6_par0, f6_par1);
  double lamda_corr = 1./(f6->GetParameter(1)*drift_velocity); // tau_inv_dc in 1/ticks
  std::cout << "Drift velocity " << drift_velocity << std::endl;

  TH2D* h2_finecorrcharge_etrue = new TH2D("h2_finecorrcharge_etrue", Form("%s;%s;%s", "", "", ""),
                                           200, 14, 26, 200, min_charge*0.9, 11000); // AS IT IS; THERE?S AN EXTRA FACTOR 3!
  double eminus_to_charge = f4_par0/ f6_par0; // Conversion factor from e- to charge (ADC x ticks)
  for (const auto& [electrons, edep_x, edep, e_true] : true_finecalib_info) {
    double corrected_charge = 0.;
    for (size_t idx_edep = 0; idx_edep < edep.size(); idx_edep++) {
      corrected_charge += electrons[idx_edep] * exp(edep_x[idx_edep] * lamda_corr);
    }
    corrected_charge *= eminus_to_charge; // Convert to charge (ADC x ticks)
    h2_finecorrcharge_etrue->Fill(e_true, corrected_charge);
  }


  //-----------------------------------------------------------------------------------------------------------------

  TEfficiency* he_expreco = new TEfficiency(*h_expreco,*h_exp);
  he_expreco->SetTitle("Hit Probability;Expected #Pe;Reco Probability");
  he_expreco->SetName("he_expreco");

  TEfficiency* he_truereco = new TEfficiency(*h_truereco,*h_true);
  he_truereco->SetTitle("True Hit Probability;True #Pe;Reco Probability");
  he_truereco->SetName("he_truereco");


  // --- SAVE -------------------------------------------------------------------------------------------------------
  std::string out_file_name = "debug_wfs.root";
  TFile* out_file = TFile::Open(out_file_name.c_str(), "RECREATE");

  h2_fail_Etrue_opdet->Write();
  h2_fail_Xtrue_opdet->Write();
  h2_fail_Ytrue_opdet->Write();
  h2_fail_Ztrue_opdet->Write();
  h2_noreco_Etrue_Xtrue->Write();
  h2_noreco_Xtrue_Ytrue->Write();
  h2_noreco_Ytrue_Ztrue->Write();
  h_exp_opdet->Write();
  h2_charge_etrue->Write();
  h_reco_opdet->Write();
  h_exp->Write();
  he_expreco->Write();
  he_truereco->Write();
  h_expreco->Write();
  h_true->Write();
  h_truereco->Write();
  h2_edep_etrue->Write();
  h2_reco_exp->Write();
  h2_true_exp->Write();
  h2_chargeperenergytrue_drifttime->Write();
  h2_chargeperenergydep_drifttime->Write();
  // hp_chargeperenergytrue_drifttime->Write();
  g_chargeperenergytrue_drifttime->Write();
  h2_finecorrcharge_etrue->Write();
  h2_dcperde_dt->Write();
  g_dcperde_dt->Write();
  h2_true_expedep->Write();
  h2_reco_true->Write();
  h2_corrcharge_etrue->Write();
  hp_corrcharge_etrue->Write();
  hp_expreco->Write();
  hp_true_exp->Write();
  h_reco_true_Prof->Write();
  h2_exp_hittime->Write();
  h_dxtrack->Write();
  h_dxtrack_outliers->Write();
  out_file->Close();

  return;
}


