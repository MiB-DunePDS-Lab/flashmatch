#include <vector>

#include <TH1.h>
#include <TMath.h>


// --- HANDLE VISIBILITIES ----------------------------------------------------

// From the hgridX histograms created by the PhotonVisibilityExpoort_module.cc, it returns
// a vector where the z-th + y-th * hgrid[2]->GetNbinsX() + x-th * hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()
// entry value is the number of the photoVisMap entry corresponding to the x,y,z coordinates.
std::vector<int> GetCryoToTPCMap(TH1D* hgrid[3], double tpc_min[3], double tpc_max[3]){
  size_t n_total_entries = hgrid[0]->GetNbinsX() * hgrid[1]->GetNbinsX() * hgrid[2]->GetNbinsX();
  std::vector<int> total_to_tpc(n_total_entries, -1);
  int tpc_entry = 0;
  for (int ix=0; ix<hgrid[0]->GetNbinsX(); ix++){
    for (int iy=0; iy<hgrid[1]->GetNbinsX(); iy++){
      for (int iz=0; iz<hgrid[2]->GetNbinsX(); iz++){
        double x = hgrid[0]->GetBinCenter(ix);
        double y = hgrid[1]->GetBinCenter(iy);
        double z = hgrid[2]->GetBinCenter(iz);
        
        if (x > tpc_min[0] && x < tpc_max[0] && 
            y > tpc_min[1] && y < tpc_max[1] && 
            z > tpc_min[2] && z < tpc_max[2]) {
          total_to_tpc[iz + iy*hgrid[2]->GetNbinsX() + ix*hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()] = tpc_entry;
          tpc_entry++;
        }
        
      }
    }
  }
  return total_to_tpc;
}

bool isInFiducialVolume(double vertex_coor[3], double vol_min[3], double vol_max[3], double x_cut){
  return (vertex_coor[0] > vol_min[0] && vertex_coor[0] < vol_max[0] &&
          vertex_coor[1] > vol_min[1] && vertex_coor[1] < vol_max[1] &&
          vertex_coor[2] > vol_min[2] && vertex_coor[2] < vol_max[2] &&
          abs(vertex_coor[0]) > x_cut);
}

int GetTPCIndex(double vertex_coor[3], TH1D* hgrid[3], std::vector<int>& total_to_tpc){
  int bin_x = hgrid[0]->FindBin(vertex_coor[0]);
  int bin_y = hgrid[1]->FindBin(vertex_coor[1]);
  int bin_z = hgrid[2]->FindBin(vertex_coor[2]);

  int tpc_index = total_to_tpc[bin_z + bin_y*hgrid[2]->GetNbinsX()
                               + bin_x*hgrid[2]->GetNbinsX()*hgrid[1]->GetNbinsX()];

  return tpc_index; // Inside TPC
}


// --- FIT FUNCTIONS ----------------------------------------------------------

// My crystalball function: alpha, n, sigma, mu
// double my_crystalball(double* x, double* par){
//   return ROOT::Math::crystalball_pdf(x[0], par[0], par[1], par[2], par[3]);
// }

// Sigmoid function
double sigmoid(double* x, double* par){
  return 1 / (1 + exp(-(x[0]-par[0])/par[1]));
}

// Positive error function (0 to 1)
double positive_erf(double* x, double* par){
  return 0.5*(1.+TMath::Erf((x[0]-par[0])/par[1]));
}

// Sigmoid times positive_erf
double sigmoid_erf(double* x, double* par){
  return sigmoid(x, &par[0]) * positive_erf(x, &par[2]);
}

// Sigmoid times sigmoid times positive_erf 
double sigmoid_sigmoid_erf(double* x, double* par){
  return sigmoid(x, &par[0]) * sigmoid(x, &par[2]) * positive_erf(x, &par[4]);
}

// Sigmoid times my_crystalball
// double sigmoid_crystalball(double* x, double* par){
//   return par[0] * sigmoid(x, &par[1]) * my_crystalball(x, &par[3]);
// }

double langaufun(double *x, double *par) {
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = 0.;       // Landau maximum location
 
      // Control constants
      double np = 500.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

// langaufun * sigmoid_sigmoid_erf
double langaus2(double* x, double* par){
  return langaufun(x, &par[0]) * sigmoid_sigmoid_erf(x, &par[4]);
}
