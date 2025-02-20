/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : rebin_vmap.cc
 * @created     : Wednesday Feb 19, 2025 04:53:35 CST
 */

#include <iostream>
#include "TFile.h"
#include "THnSparse.h"
#include "DUNEVisUtils.hpp"

int rebin_vmap(const TString input_vmap_file, const TString output_vmap_file)
{
  const int N_OPDET = 480; 
  TFile* finput = TFile::Open(input_vmap_file); 

  TFile* fout = new TFile(output_vmap_file, "recreate"); 
  
  for (int iopdet = 0; iopdet < N_OPDET; iopdet++) {
    printf("processing OpDet %i\n", iopdet); 
    TString hname = Form("h3VisMap_opDet%i", iopdet); 
    THnSparseF* hog = finput->Get<THnSparseF>(hname); 

    THnSparseF* hnew = static_cast<THnSparseF*>( rebin_visibility_map(hog, 2, 2, 2) );
    fout->cd(); 
    hnew->Write(); 

    delete hnew;
  }
  return 0;
}

