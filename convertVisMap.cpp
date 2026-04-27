#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include <cstddef>
#include <string>
#include "TTreeReaderArray.h"



int n_opdet = 184; // Number of optical detectors
const std::string visibility_file_name_ar = "./visibility/raw/dunevis_fdvd_1x8x14_ar.root"; // File with the visibility maps
const std::string visibility_file_name_xe = "./visibility/raw/dunevis_fdvd_1x8x14_xe.root"; // File with the visibility maps
TString output_file_name = "./visibility/processed/dunevis_dunevd10kt.root"; // Output file name
// TString output_file_name = "./visibility/processed/dunevis_dune10kt.root"; // Output file name
const bool verbose = true; // Print progress messages
// const bool add_reflected = false; // Add reflected light to direct light in the output maps

void convertVisMap(){
  TFile* visibility_file_ar = TFile::Open(visibility_file_name_ar.c_str(), "READ");
  TTree* photovisAr = (TTree*)visibility_file_ar->Get("photovisAr/photoVisMap");
  TH1D* hgrid[3] = {nullptr};
  hgrid[0] = (TH1D*)visibility_file_ar->Get("photovisAr/hgrid0")->Clone("hgrid0");
  hgrid[1] = (TH1D*)visibility_file_ar->Get("photovisAr/hgrid1")->Clone("hgrid1");
  hgrid[2] = (TH1D*)visibility_file_ar->Get("photovisAr/hgrid2")->Clone("hgrid2");

  TFile* visibility_file_xe = nullptr;
  TTree* photovisXe = nullptr;
  if (visibility_file_name_xe != "") {
    visibility_file_xe = TFile::Open(visibility_file_name_xe.c_str(), "READ");
    photovisXe = (TTree*)visibility_file_xe->Get("photovisXe/photoVisMap");
  }
  
  TFile* out_file = new TFile(output_file_name, "RECREATE");
  out_file->SetCompressionLevel(0);
  out_file->cd();
  hgrid[0]->Write();
  hgrid[1]->Write();
  hgrid[2]->Write();

  TTree* tDimensions = (TTree*)visibility_file_ar->Get("photovisAr/tDimensions");
  TTree* out_tDimensions = tDimensions->CloneTree();
  out_tDimensions->Write();


  TTree* out_tree = new TTree("photoVisMap", "photoVisMap");
  float vis_dir_ar[n_opdet]; float vis_ref_ar[n_opdet]; float vis_dir_xe[n_opdet]; float vis_ref_xe[n_opdet];
  out_tree->Branch("opDet_visDirect_Ar", vis_dir_ar, Form("opDet_visDirect_Ar[%d]/F", n_opdet));
  out_tree->Branch("opDet_visReflct_Ar", vis_ref_ar, Form("opDet_visReflct_Ar[%d]/F", n_opdet));
  if (photovisXe) {
    out_tree->Branch("opDet_visDirect_Xe", vis_dir_xe, Form("opDet_visDirect_Xe[%d]/F", n_opdet));
    out_tree->Branch("opDet_visReflct_Xe", vis_ref_xe, Form("opDet_visReflct_Xe[%d]/F", n_opdet));
  }

  const size_t n_entriesmap = photovisAr->GetEntries();
  TTreeReader VisMapReader_ar(photovisAr);
  TTreeReaderArray<double> opDet_visDirect_reader_ar(VisMapReader_ar, "opDet_visDirect");
  TTreeReaderArray<double> opDet_visReflct_reader_ar(VisMapReader_ar, "opDet_visReflct");

  TTreeReader VisMapReader_xe = TTreeReader(photovisXe);
  TTreeReaderArray<double> opDet_visDirect_reader_xe(VisMapReader_xe, "opDet_visDirect");
  TTreeReaderArray<double> opDet_visReflct_reader_xe(VisMapReader_xe, "opDet_visReflct");
  if (photovisXe) {
    // VisMapReader_xe = TTreeReader(photovisXe);
    opDet_visDirect_reader_xe = TTreeReaderArray<double>(VisMapReader_xe, "opDet_visDirect");
    opDet_visReflct_reader_xe = TTreeReaderArray<double>(VisMapReader_xe, "opDet_visReflct");
  }

  std::cout << "Filling maps..." << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  int VisMapEntry = 0; int n_VisMapEntries = photovisAr->GetEntries();
  while (VisMapReader_ar.Next() && (!photovisXe || VisMapReader_xe.Next())) {
    for (size_t j = 0; j < n_opdet; ++j) {
      vis_dir_ar[j] = static_cast<float>(opDet_visDirect_reader_ar[j]);
      vis_ref_ar[j] = static_cast<float>(opDet_visReflct_reader_ar[j]);
      if (photovisXe) {
        vis_dir_xe[j] = static_cast<float>(opDet_visDirect_reader_xe[j]);
        vis_ref_xe[j] = static_cast<float>(opDet_visReflct_reader_xe[j]);
      }
    }
    VisMapEntry++;
    out_tree->Fill();
    if (verbose && VisMapEntry % 5000 == 0)
      std::cout << VisMapEntry << " / " << n_VisMapEntries << "\r" << std::flush;
  }
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;
  std::cout << "Done reading maps in " << elapsed_time.count() << " seconds." << std::endl;
  std::cout << "Done filling maps..." << std::endl;
  
  out_tree->Write("", TObject::kOverwrite);
  out_file->Close();
  visibility_file_ar->Close();
  if (visibility_file_xe) visibility_file_xe->Close();

  return;
}
