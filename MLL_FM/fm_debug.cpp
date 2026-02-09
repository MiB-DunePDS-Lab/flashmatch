

#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
void fm_debug(){
  TFile* distribution_file = TFile::Open("fm_distributions.root", "READ");
  TTree* tpc_pds_tree = static_cast<TTree*>(distribution_file->Get("tpc_pds_tree"));
  std::vector<float>* reco_pes = nullptr;
  std::vector<float>* exp_phs = nullptr;
  tpc_pds_tree->SetBranchAddress("reco_pes", &reco_pes);
  tpc_pds_tree->SetBranchAddress("exp_phs", &exp_phs);
  // TH2D* h2_exp_reco = static_cast<TH2D*>(distribution_file->Get("h2_exp_reco_norm"));

  TFile* parametrizer_file = TFile::Open("fm_parametrizer.root", "READ");
  TH2D* h2_exp_reco = static_cast<TH2D*>(parametrizer_file->Get("h2_exp_reco"));

  float reco_pe, exp_ph;
  for (Long64_t entry = 0; entry < tpc_pds_tree->GetEntries(); entry++) {
    tpc_pds_tree->GetEntry(entry);
    
    for (size_t i=0; i<reco_pes->size(); i++){
      exp_ph = exp_phs->at(i);
      reco_pe= reco_pes->at(i);
      if (reco_pe>0.0){
        if (h2_exp_reco->GetBinContent(h2_exp_reco->FindBin(reco_pe, exp_ph)) == 0){
          std::cout << "Warning: h2_exp_reco bin content is zero for reco_pe=" << reco_pe
                    << " and exp_ph=" << exp_ph << std::endl;
        }

      }
    }


  }
  
  return;
}
