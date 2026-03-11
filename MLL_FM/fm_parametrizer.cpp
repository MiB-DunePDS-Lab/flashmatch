#include "TGraphErrors.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TEfficiency.h>
#include <TROOT.h>
#include "Utils.hpp"


void fm_parametrizer(){
  gROOT->SetBatch(kTRUE);
  // --- CONFIGS ---------------------------------------------------------------
  MLLcconfigs f = load_ana_config("./config.json");
  std::string input_dir    = f.input_dir;
  std::string distribution = f.distribution;
  double fit_trend_low     = f.fit_trend_low;
  double fit_trend_up      = f.fit_trend_up;
  bool rebin_h2            = f.rebin_h2;
  int rebinning            = f.rebinning;
  
  // --- INPUTS ---------------------------------------------------------------
  TFile* distribution_file = TFile::Open((input_dir+"MLL_Distributions.root").c_str(), "READ");
  TH2D* h2_exp_reco        = static_cast<TH2D*>(distribution_file->Get("h2_exp_reco"));
  TEfficiency* he_hit_prob = static_cast<TEfficiency*>(distribution_file->Get("he_hit_prob"));
 

  // --- FIT ------------------------------------------------------------------
  // Fit reco probability -----------------------------------------------------  
  TF1* f_reco_prob_fit = new TF1("f_hit_prob_fit", sigmoid_sigmoid_erf, 0., 60., 6); // To fit only in a short range
  f_reco_prob_fit->SetParNames("amp1","x0_1","slope1","amp2","x0_2","slope2");
  f_reco_prob_fit->SetNpx(3000);
  f_reco_prob_fit->SetParameters(7, 0.5, 1.0, 7., 0.8, 5.0);
  f_reco_prob_fit->SetParLimits(2, 0.75, 1.05);     
  he_hit_prob->Fit(f_reco_prob_fit, "R");
  
  TF1* f_reco_prob = new TF1("f_reco_prob", sigmoid_sigmoid_erf, 0., 2000., 6); // Extend the function range
  f_reco_prob->SetNpx(3000);
  f_reco_prob->SetParameters(f_reco_prob_fit->GetParameters());

  // Create and fit the reco|expected distributions ---------------------------
  TF1* f_RecoExpDistr = nullptr;
  if (distribution == "lognormal") {
    f_RecoExpDistr = new TF1("f_RecoExpDistr", "ROOT::Math::lognormal_pdf(x,[0],[1])", 0, 2000);
    f_RecoExpDistr->SetParNames("log(m)", "#sigma");
  }
  else {
    f_RecoExpDistr = new TF1("f_RecoExpDistr", gamma_dist, 0., 2000., 2);
    f_RecoExpDistr->SetParNames("k", "theta");
  }
  f_RecoExpDistr->SetNpx(3000);

  std::vector<double> par1s, par2s, exp_phs;
  std::vector<double> err_par1, err_par2, err_exp_phs;
  
  TFile* out_file = TFile::Open((input_dir+"MLL_Parametrizer.root").c_str(), "RECREATE");
  out_file->mkdir("projections");
  out_file->cd("projections");
  
  if (rebin_h2) h2_exp_reco->Rebin2D(rebinning, rebinning);
  for(int idx_y=1; idx_y<=h2_exp_reco->GetNbinsY(); idx_y++){
    std::cout << idx_y << "/" << h2_exp_reco->GetNbinsY() << "\r" << std::flush;
    TH1D* h1_proj = h2_exp_reco->ProjectionX("h1_proj", idx_y, idx_y, "e");
    double h1_proj_integral = h1_proj->Integral();
    
    if(h1_proj_integral > 0.){
      h1_proj->Scale(1./(h1_proj_integral*h1_proj->GetBinWidth(1)));
      
      double mean   = h1_proj->GetBinCenter(h1_proj->GetMaximumBin());
      double stddev = h1_proj->GetStdDev();

      if (distribution == "lognormal") {
        f_RecoExpDistr->SetParameters(log(mean), 0.25);
        f_RecoExpDistr->SetParLimits(0, log(mean-5*stddev), log(mean+5*stddev));
      }
      else {
        f_RecoExpDistr->SetParameters(mean, 1.);
        f_RecoExpDistr->SetParLimits(0, mean-5*stddev, mean+5*stddev);
      }

      TFitResultPtr fit_res = nullptr;
      if (h1_proj->GetEntries() > 400) fit_res = h1_proj->Fit(f_RecoExpDistr, "Q", "", std::max(0.,mean-5*stddev), mean+5*stddev);
      h1_proj->SetName(Form("h1_proj_%i_pe", (int)h2_exp_reco->GetYaxis()->GetBinCenter(idx_y)));
      h1_proj->Write();

      for(int idx_x=1; idx_x<=h2_exp_reco->GetNbinsX(); idx_x++){
        h2_exp_reco->SetBinContent(idx_x, idx_y, h1_proj->GetBinContent(idx_x));
      } 

      
      double exp_ph = h2_exp_reco->GetYaxis()->GetBinCenter(idx_y);
      // if (reco > 400) break; // Stop if reco is greater than 40
      if (exp_ph > 0 && h1_proj->GetEntries() > 400 && fit_res==0){
        exp_phs.push_back(exp_ph);                       err_exp_phs.push_back(0);
        par1s.push_back(f_RecoExpDistr->GetParameter(0));   err_par1.push_back(f_RecoExpDistr->GetParError(0));
        par2s.push_back(f_RecoExpDistr->GetParameter(1));  err_par2.push_back(f_RecoExpDistr->GetParError(1));
      }
    }
  } // handle h2_exp_reco and fit lognormal_pdf


  // --- TGRAPHS --------------------------------------------------------------
  TGraphErrors* g_par1 = new TGraphErrors(par1s.size(), &exp_phs[0], &par1s[0], &err_exp_phs[0], &err_par1[0]);
  g_par1->SetTitle("log(MPV) vs Expected Photons;#Expected Photons;log(MPV)");
  g_par1->SetName("g_par1");

  TGraphErrors* g_par2 = new TGraphErrors(par2s.size(), &exp_phs[0], &par2s[0], &err_exp_phs[0], &err_par2[0]);
  g_par2->SetTitle("#sigma vs Expected Photons;#Expected Photons;#sigma");
  g_par2->SetName("g_par2");
  
  // Fit the parameter trends -------------------------------------------------
  TF1* f_par1_trend = nullptr;
  TF1* f_par2_trend = nullptr;
  if (distribution == "lognormal") {
    f_par1_trend = new TF1("f_par1_trend", "[0]+[1]*log(x+[2])", 0., 2000);
    f_par1_trend->SetParNames("offset","log_slope","linear_slope");
    f_par1_trend->SetParameters(-0.5,1.2,1.0);  
    f_par1_trend->SetParLimits(1, 0.0, 5.0);
    f_par1_trend->SetParLimits(2,-3.0, 10.0);
  
    f_par2_trend = new TF1("f_par2_trend","([0] + [1]*x + [2]*x*x)/(1 + [3]*x + [4]*x*x)", 0., 2000);
    f_par2_trend->SetParameters(0.6, 0.01,0.0001,0.05,0.0002);
    f_par2_trend->SetParLimits(2,0,10);    
    f_par2_trend->SetParLimits(4,0,10);   
  }
  else {
    f_par1_trend = new TF1("f_par1_trend", "[0]+[1]*x", 0., 2000.);
    f_par1_trend->SetParNames("const", "log_slope");
    f_par1_trend->SetParameters(.5, 1.2);
  
    f_par2_trend = new TF1("f_par2_trend", "[2]+[1]*sqrt([0]*x)", 0., 2000.);
    f_par2_trend->SetParNames("b", "A", "const");
    f_par2_trend->SetParameters(2., 1., -0.5);
  }
  f_par1_trend->SetNpx(3000); f_par2_trend->SetNpx(3000);
  g_par1->Fit(f_par1_trend, "", "", fit_trend_low, fit_trend_up); 
  g_par2->Fit(f_par2_trend, "", "", fit_trend_low, fit_trend_up);


  out_file->cd();
  he_hit_prob->Write();
  h2_exp_reco->Write();
  f_reco_prob->Write();
  f_RecoExpDistr->SetParameters(f_par1_trend->Eval(20.), f_par2_trend->Eval(20.));
  f_RecoExpDistr->Write();
  f_par1_trend->Write();
  f_par2_trend->Write();
  g_par1->Write();
  g_par2->Write();
  out_file->Close();
  
  distribution_file->Close();

  return;
}
