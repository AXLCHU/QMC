#ifndef __PATH_GENERATION_HPP
#define __PATH_GENERATION_HPP

#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include "Sampling.h"

using namespace std;

// GBM model

void GBM_paths(std::vector<double>& spot_prices, const double& r, const double& v, const double& T) { 

  double dt = T/static_cast<double>(spot_prices.size());
  double drift = exp(dt*(r-0.5*v*v));
  double vol = sqrt(v*v*dt);

  for (int i=1; i<spot_prices.size(); i++) {
    double gauss_bm = gaussian_box_muller(); //
    spot_prices[i] = spot_prices[i-1] * drift * exp(vol*gauss_bm);
  }
}


void generate_Loss(vector<double>& Loss, vector<double>& spot_prices, const double& r, const double& v, const double& T,
    const double& K) {

    double dt = T / static_cast<double>(spot_prices.size());
    double drift = exp(dt * (r - 0.5 * v * v));
    double vol = sqrt(v * v * dt);
    double N = spot_prices.size();
    double dim = Loss.size();
    
    for (int n = 0; n < dim; n++) {
        
        for (int k = 1; k < N; k++) {

            double gauss_bm = gaussian_box_muller(); //

            spot_prices[k] = spot_prices[k - 1] * drift * exp(vol * gauss_bm);
        }

        Loss[n] = std::max(K - spot_prices[N - 1], 0.0) - exp(r * T) * BSM_Put(spot_prices[0], r, v, T, K);
    }
  
}



// Heston model

void Heston_paths(std::vector<double>& spot_prices, const double& r, const double& T,
                            const double& xi, const double& kappa, const double& theta, const double& rho,
                            const int& num_intervals, const double& X0)
{ 
    double dt = T / static_cast<double>(spot_prices.size());
//    std::vector<double> vol_draws(num_intervals, 0);
//    std::vector<double> spot_draws(num_intervals, 0);
    std::vector<double> vol_path(num_intervals, X0);

    for (int i = 1; i < spot_prices.size(); i++) {
 
        double gauss_bm1 = gaussian_box_muller();
        double gauss_bm2 = gaussian_box_muller();
        gauss_bm2 = rho * gauss_bm1 + sqrt((1 - rho)) * gauss_bm2;

        double v_max = std::max(vol_path[i - 1], 0.0);
        
        vol_path[i] = vol_path[i - 1] + kappa * (theta - vol_path[i - 1]) * dt + xi * sqrt(v_max * dt) * gauss_bm1; // CIR

        spot_prices[i] = spot_prices[i - 1] * exp((r - 0.5 * v_max * v_max) * dt + v_max * sqrt(dt) * gauss_bm2); // GBM
    }
}

#endif