#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TTreeReader.h"
#include <cstddef>
#include <string>
#include "TTreeReaderArray.h"



int n_opdet = 480; // Number of optical detectors
const std::string visibility_file_name = "dunevis_fdhd_1x2x6_test.root"; // File with the visibility maps
const bool verbose = true; // Print progress messages

void convertVisMap(){
  TFile* visibility_file = TFile::Open(visibility_file_name.c_str(), "READ");
  TTree* photoVisMap = (TTree*)visibility_file->Get("photovisAr/photoVisMap");
  TH1D* hgrid[3] = {nullptr};
  hgrid[0] = (TH1D*)visibility_file->Get("photovisAr/hgrid0")->Clone("hgrid0");
  hgrid[1] = (TH1D*)visibility_file->Get("photovisAr/hgrid1")->Clone("hgrid1");
  hgrid[2] = (TH1D*)visibility_file->Get("photovisAr/hgrid2")->Clone("hgrid2");
  
  TFile* out_file = new TFile("dunevis_fdhd_1x2x6_test_float.root", "RECREATE");
  out_file->SetCompressionLevel(0);

  out_file->mkdir("photovisAr");
  out_file->cd("photovisAr");
  hgrid[0]->Write();
  hgrid[1]->Write();
  hgrid[2]->Write();

  TTree* tDimensions = (TTree*)visibility_file->Get("photovisAr/tDimensions");
  TTree* out_tDimensions = tDimensions->CloneTree();
  out_tDimensions->Write();


  TTree* out_tree = new TTree("photoVisMap", "photoVisMap");
  float floatVisMap[n_opdet];
  out_tree->Branch("opDet_visDirect", floatVisMap, Form("opDet_visDirect[%d]/F", n_opdet));
  

  const size_t n_entriesmap = photoVisMap->GetEntries();
  TTreeReader VisMapReader(photoVisMap);
  TTreeReaderArray<double> opDet_visDirect_reader(VisMapReader, "opDet_visDirect");

  std::cout << "Filling maps..." << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  int VisMapEntry = 0; int n_VisMapEntries = photoVisMap->GetEntries();
  while (VisMapReader.Next()){
    for (size_t j = 0; j < n_opdet; ++j) {
      floatVisMap[j] = static_cast<float>(opDet_visDirect_reader[j]);
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
  
  out_tree->Write();
  out_file->Close();
  visibility_file->Close();

  return;

}
