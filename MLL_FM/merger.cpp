#include <vector>
#include <filesystem>

#include "TFileMerger.h"

namespace fs = std::filesystem;

inline std::vector<std::string> get_list_of_files_in_folder(std::string folder_name, std::string ends_with){
  std::vector<std::string> result;

  if (!fs::exists(folder_name) || !fs::is_directory(folder_name))
    return result;

  for (const auto& entry : fs::directory_iterator(folder_name))
  {
    if (!entry.is_regular_file())
      continue;

    std::string filename = entry.path().filename().string();

    if (filename.size() >= ends_with.size() &&
      filename.compare(filename.size() - ends_with.size(),
                       ends_with.size(),
                       ends_with) == 0)
    {
      result.push_back(filename);
    }
  }

  return result;
}

std::string sample_dir = "../data/testing_cheating_FM/no_bkg/files/";

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
