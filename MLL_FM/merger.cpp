#include <vector>

#include "TFileMerger.h"

#include "Utils.hpp"

std::string sample_dir = "../data/testing_MLL_FM/files/";

void merger(){
  std::vector<std::string> ana_files = get_list_of_files_in_folder(sample_dir, ".root");

  TFileMerger* tfmerger = new TFileMerger();
  // Loop over ana files and merge them into one ROOT TFile
  for (const auto& ana_file : ana_files){
    std::string ana_file_name = sample_dir+ana_file;
    if(!std::filesystem::exists(ana_file_name)) continue;
    tfmerger->AddFile(ana_file_name.c_str());
  }
  tfmerger->OutputFile("SolarNuAnaMerged.root");
  tfmerger->Merge();

  return;
}
