#include<stdlib.h>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<climits>
#include<list>
#include<algorithm>
#include<numeric>
#include<vector>
#include<iostream>
#include<random>
#include<chrono>
#include "VaR_CVaR.h"
#include "BSM.hpp"
#include "Path_generation.hpp"
#include "Low_discrepancy_sequences.h"
#include "Inverse_fcts.h"
#include <boost/random/sobol.hpp>

using namespace std;


int main() {

/*	cout << rand;
	cout << "\n" << expo(rand,1);
	cout << "\n" << gauss(rand);*/

	double alpha = 0.95;
	double dim = 20000;

	std::vector<double> MT(dim, 0);	std::vector<double> RN(dim, 0);
	std::vector<double> gausss(dim, 0);
	std::vector<double> halton1(dim, 0);
	std::vector<double> halton2(dim, 0);
	std::vector<double> X(dim,0); std::vector<double> X1(dim, 0);
	std::vector<double> Y(dim, 0); std::vector<double> Y1(dim, 0);
	std::vector<double> Z(dim, 0);
	std::vector<double> Z1(dim, 0);
	std::vector<double> VaR(dim, 10);
	std::vector<double> CVaR(dim, 10);

	std::normal_distribution<double> dist1{ 0, 1 }; //

	for (int i = 0; i < dim; i++) {

		halton1[i] = Halton_seq(i + 1, 2);
		halton2[i] = Halton_seq(i + 1, 3); //

		//std::random_device seed; // uniformly-distributed integer random number generator
		//std::mt19937 gen{ seed() };
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::mt19937 mt1(seed);
		MT[i] = mt1();
		RN[i] = mt1() / static_cast<double>(mt1.max()); 

		X[i] = gauss(halton1[i]);
		X1[i] = gauss(RN[i]);
//		X[i] = gauss(rand()/ static_cast<double>(RAND_MAX) - 1);

		unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine gen(seed1); // linear congruential engine
		Y[i] = dist1(gen);

		Z[i] = gaussian_box_muller();
		//Z1[i] = gaussian_box_muller_v2(rand1[i]);
	//	gausss[i] = gauss2(halton2[i], halton1[i]);
		
	}


	VaR_CVaR(CVaR, VaR, Z, alpha);

	double variance = 0;
	double sum = accumulate(VaR.begin(), VaR.end(), 0);
	double mean = sum / VaR.size();

	for (int n = 0; n < dim; n++) {
		variance += (VaR[n] - mean) * (VaR[n] - mean) / VaR.size();
	}
	double std = sqrt(variance);

	cout << "\nStandand Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1]; 
	cout << " & std dev de VaR avec Box-Muller  = " << std * 100 << "%";
	cout << "\nStandand Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];

	VaR_CVaR(CVaR, VaR, Y, alpha);

	for (int n = 0; n < dim; n++) {
		variance += (VaR[n] - mean) * (VaR[n] - mean) / VaR.size();
	}
	std = sqrt(variance);

	cout << "\n\nStandand Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR avec Y = " << std * 100 << "%";
	cout << "\nCVaR = " << CVaR[dim - 1];

	//
	VaR_CVaR(CVaR, VaR, X1, alpha);

	for (int n = 0; n < dim; n++) {
		variance += (VaR[n] - mean) * (VaR[n] - mean) / VaR.size();
	}
	std = sqrt(variance);

	cout << "\n\nStandand Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR avec X1 = " << std * 100 << "%";
	cout << "\nCVaR = " << CVaR[dim - 1];

	double S0 = 90;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	unsigned num_intervals = 250;

	cout << "\n\nFor S0 = " << S0 << ", K = " << K << ", r = " << r << " and v = " << v << " : ";
	cout << "\n\nPut option price = " << BSM_Put(S0,r,v,T,K);
	cout << "\nCall option price = " << BSM_Call(S0, r, v, T, K);

	vector<double> spot_prices(num_intervals, S0);
	vector<double> loss_fct(dim);

	for (int n = 0; n < dim; n++) {
		GBM_paths(spot_prices, r, v, T);
		loss_fct[n] = std::max(K - spot_prices[num_intervals - 1],0.0) - exp(r*T) * BSM_Put(S0, r, v, T, K); ///
	}

	VaR_CVaR(CVaR, VaR, loss_fct, alpha);

	cout << "\n\nVaR at "<< alpha*100 << "% for Put option = " << VaR[dim-1];
	cout << "\nCVaR at " << alpha * 100 << "% for Put option  = " << CVaR[dim - 1];
	cout << "\n";

	return 0;
}

