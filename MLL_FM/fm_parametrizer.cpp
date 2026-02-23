#include "TGraphErrors.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TEfficiency.h>
#include "Utils.hpp"


void fm_parametrizer(){
  // --- CONFIGS ---------------------------------------------------------------
  MLLcconfigs f = load_ana_config("./config.json");
  std::string input_dir = f.input_dir;
  double fit_trend_low = f.fit_trend_low;
  double fit_trend_up  = f.fit_trend_up;
  bool rebin_h2        = f.rebin_h2;
  int rebinning        = f.rebinning;
  
  // --- INPUTS ---------------------------------------------------------------
  TFile* distribution_file = TFile::Open((input_dir+"MLL_Distributions.root").c_str(), "READ");
  TH2D* h2_exp_reco        = static_cast<TH2D*>(distribution_file->Get("h2_exp_reco"));
  TEfficiency* he_hit_prob = static_cast<TEfficiency*>(distribution_file->Get("he_hit_prob"));
  
  // ========================
  //      Fit Prob
  // ========================
  
  TF1* f_reco_prob_fit = new TF1("f_hit_prob_fit", sigmoid_sigmoid_erf, 0., 60., 6); // To fit only in a short range
  f_reco_prob_fit->SetParNames("amp1","x0_1","slope1","amp2","x0_2","slope2");
  f_reco_prob_fit->SetNpx(3000);
  f_reco_prob_fit->SetParameters(7, 0.5, 1.0, 7., 0.8, 5.0);
  f_reco_prob_fit->SetParLimits(2, 0.75, 1.05);     
  he_hit_prob->Fit(f_reco_prob_fit, "R");
  
  TF1* f_reco_prob = new TF1("f_reco_prob", sigmoid_sigmoid_erf, 0., 2000., 6); // Extend the function range
  f_reco_prob->SetNpx(3000);
  f_reco_prob->SetParameters(f_reco_prob_fit->GetParameters());

  TF1* f_lognormal = new TF1("f_lognormal", "ROOT::Math::lognormal_pdf(x,[0],[1])", 0, 2000);
  f_lognormal->SetParNames("log(m)", "#sigma");
  f_lognormal->SetNpx(3000);
  std::vector<double> logms, sigmas, reco_pes;
  std::vector<double> err_logms, err_sigmas, err_reco_pes;
  
  TFile* out_file = TFile::Open((input_dir+"MLL_Parametrizer.root").c_str(), "RECREATE");
  out_file->mkdir("projections");
  out_file->cd("projections");
  
  if (rebin_h2) h2_exp_reco->Rebin2D(rebinning, rebinning);
  for(int idx_y=1; idx_y<=h2_exp_reco->GetNbinsY(); idx_y++){
    std::cout << idx_y << "/" << h2_exp_reco->GetNbinsY() << "\r" << std::flush;
    TH1D* h1_proj = h2_exp_reco->ProjectionX("h1_proj", idx_y, idx_y);
    double h1_proj_integral = h1_proj->Integral();
    
    if(h1_proj_integral > 0.){
      h1_proj->Scale(1./(h1_proj_integral*h1_proj->GetBinWidth(1)));
      
      double mean   = h1_proj->GetBinCenter(h1_proj->GetMaximumBin());
      double stddev = h1_proj->GetStdDev();

      f_lognormal->SetParameters(log(mean), 0.25);
      f_lognormal->SetParLimits(0, log(mean-5*stddev), log(mean+5*stddev));

      TFitResultPtr fit_res = nullptr;
      if (h1_proj->GetEntries() > 400) fit_res = h1_proj->Fit(f_lognormal, "Q", "", std::max(0.,mean-5*stddev), mean+5*stddev);
      h1_proj->SetName(Form("h1_proj_%i_pe", (int)h2_exp_reco->GetYaxis()->GetBinCenter(idx_y)));
      h1_proj->Write();

      for(int idx_x=1; idx_x<=h2_exp_reco->GetNbinsX(); idx_x++){
        h2_exp_reco->SetBinContent(idx_x, idx_y, h1_proj->GetBinContent(idx_x));
      } 

      
      double reco = h2_exp_reco->GetXaxis()->GetBinCenter(idx_y);
      // if (reco > 400) break; // Stop if reco is greater than 40
      if (reco > 3 && h1_proj->GetEntries() > 400 && fit_res==0){
        reco_pes.push_back(reco);                        err_reco_pes.push_back(0);
        logms.push_back(f_lognormal->GetParameter(0));   err_logms.push_back(f_lognormal->GetParError(0));
        sigmas.push_back(f_lognormal->GetParameter(1));  err_sigmas.push_back(f_lognormal->GetParError(1));
      }
    }
  } // handle h2_exp_reco and fit lognormal_pdf


  // --- TGRAPHS ------------------------------------------------
  TGraphErrors* g_logms = new TGraphErrors(logms.size(), &reco_pes[0], &logms[0], &err_reco_pes[0], &err_logms[0]);
  g_logms->SetTitle("log(MPV) vs reco_pe;#Reco Pe;log(MPV)");
  g_logms->SetName("g_logms");

  TGraphErrors* g_sigmas = new TGraphErrors(sigmas.size(), &reco_pes[0], &sigmas[0], &err_reco_pes[0], &err_sigmas[0]);
  g_sigmas->SetTitle("#sigma vs reco_pe;#Reco Pe;#sigma");
  g_sigmas->SetName("g_sigmas");
  
  // Fit the parameter trends
  
  // ========================
  //      Fit logms
  // ========================
  
  TF1* f_logms_trend = new TF1("f_logms_trend", "[0]+[1]*log(x+[2])", 0.4, 2000);
  f_logms_trend->SetParNames("offset","log_slope","linear_slope");
  f_logms_trend->SetParameters(-0.5,1.2,1.0);  
  f_logms_trend->SetParLimits(1, 0.0, 5.0);
  f_logms_trend->SetParLimits(2,-3.0, 10.0);
  //f_logms_trend->SetNpx(3000);
  //g_logms->Fit(f_logms_trend, "R");
  g_logms->Fit(f_logms_trend, "", "", fit_trend_low, fit_trend_up); 
  
  // ========================
  //      Fit sigma
  // ========================
      
  TF1* f_sigmas_trend = new TF1("f_sigmas_trend","([0] + [1]*x + [2]*x*x)/(1 + [3]*x + [4]*x*x)", 3.0, 2000);
  f_sigmas_trend->SetParameters(0.6, 0.01,0.0001,0.05,0.0002);
  f_sigmas_trend->SetNpx(3000);
  f_sigmas_trend->SetParLimits(4,0,10);   
  f_sigmas_trend->SetParLimits(2,0,10);    
  g_sigmas->Fit(f_sigmas_trend, "", "", fit_trend_low, fit_trend_up);
  
  TTree* parametrizer_tree = new TTree("parametrizer_tree", "parametrizer_tree");
  float log_const, sigma_x0, sigma_lambda, sigma_A, sigma_const;
  parametrizer_tree->Branch("log_const",    &log_const, "log_const/F");
  parametrizer_tree->Branch("sigma_x0",     &sigma_x0, "sigma_x0/F");
  parametrizer_tree->Branch("sigma_lambda", &sigma_lambda, "sigma_lambda/F");
  parametrizer_tree->Branch("sigma_A",      &sigma_A, "sigma_A/F");
  parametrizer_tree->Branch("sigma_const",  &sigma_const, "sigma_const/F");
  log_const    = f_logms_trend->GetParameter(0);
  sigma_x0     = f_sigmas_trend->GetParameter(0);
  sigma_lambda = f_sigmas_trend->GetParameter(1);
  sigma_A      = f_sigmas_trend->GetParameter(2);
  sigma_const  = f_sigmas_trend->GetParameter(3);
  parametrizer_tree->Fill();


  out_file->cd();
  parametrizer_tree->Write();
  he_hit_prob->Write();
  h2_exp_reco->Write();
  f_reco_prob->Write();
  f_lognormal->Write();
  f_logms_trend->Write();
  f_sigmas_trend->Write();
  g_logms->Write();
  g_sigmas->Write();
  out_file->Close();
  
  distribution_file->Close();

  return;
}
