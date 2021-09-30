#include<stdlib.h>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include <string>
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
#include "Sampling.h"
#include <boost/random/sobol.hpp>
#include "MeanVariance.h"
#include <boost/math/distributions/normal.hpp>

using namespace std;


int main() {

	double alpha = 0.95;
	double dim = 30000;

	 vector<double> MT(dim, 0);	vector<double> RN(dim, 0);
	 vector<double> gausss(dim, 0);
	 vector<double> halton1(dim, 0);
	 vector<double> halton2(dim, 0);
	 vector<double> X(dim,0);  vector<double> X1(dim, 0);
	 vector<double> Y(dim, 0);  vector<double> Y1(dim, 0); vector<double> Y2(dim, 0);
	 vector<double> Z(dim, 0);
	 vector<double> Z1(dim, 0);
	 vector<double> VaR(dim, 10); // initial value
	 vector<double> CVaR(dim, 10);

	std::normal_distribution<double> normal{ 0, 1 }; //
	//std::uniform_real_distribution<double> unif(0,1);
	boost::math::normal dist(0, 1);

	double q = quantile(dist, alpha);
//	double q1 = quantile(dist, 1 - alpha);
	double qq = (1 / (1 - alpha)) * ((1 / sqrt(2 * M_PI)) * exp(-q * q * 0.5));

	cout << "\nSND quantile at " << alpha * 100 << "% = " << q;
	cout << "\nSND CVaR at " << alpha * 100 << "% = " << qq;

	gauss2(X1);
	VaR_CVaR(CVaR, VaR, X1, alpha);
	double std = StdDev(VaR);

	cout << "\n\nPseudo-random case :";
	cout << "\n\nVaR at " << alpha * 100 << "% = " << VaR[dim - 1]; cout << " & std dev de VaR = " << std * 100 << "%";
	cout << "\nCVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];


	for (int i = 0; i < dim; i++) {
		
		Z[i] = gaussian_box_muller();

		halton1[i] = Halton_seq(i + 1, 2);
		halton2[i] = Halton_seq(i + 1, 3);

		//std::random_device seed; // uniformly-distributed integer random number generator
		//std::mt19937 gen{ seed() };
		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::mt19937 mt1(seed);
		MT[i] = mt1();
		RN[i] = mt1() / static_cast<double>(mt1.max()); 

//		X[i] = gauss(halton1[i]);
//		X1[i] = gauss(RN[i]);

		unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed1); // linear congruential engine
//		mersenne_twister_engine mt_engine(123);
		Y[i] = normal(generator);

//		Y1[i] = unif(generator);
//		Y2[i] = gauss(Y1[i]);
	}

	VaR_CVaR(CVaR, VaR, Z, alpha);
	std = StdDev(VaR);
	
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1]; 
	cout << " & std dev de VaR avec Box-Muller = " << std * 100 << "%";
	cout << "\nStandard Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];
	
//	gauss_low_discrepency(X, halton1, halton2);
	gauss_low_discrepency(X,dim);
	VaR_CVaR(CVaR, VaR, X, alpha);
	std = StdDev(VaR);

	cout << "\n\nQuasi-random case:";
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR avec Halton = " << std * 100 << "%";
	cout << "\nCVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];

	/*
	VaR_CVaR(CVaR, VaR, X, alpha);
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR avec MT = " << std * 100 << "%";
	cout << "\nCVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];
	*/

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

	for (int n = 0; n < dim; n++) {
		GBM_paths(spot_prices, r, v, T);
		loss_fct[n] = std::max(spot_prices[num_intervals - 1] - K, 0.0) - exp(r * T) * BSM_Call(S0, r, v, T, K); ///
	}

	VaR_CVaR(CVaR, VaR, loss_fct, alpha);

	cout << "\n\nVaR at " << alpha * 100 << "% for Call option = " << VaR[dim - 1];
	cout << "\nCVaR at " << alpha * 100 << "% for Call option  = " << CVaR[dim - 1];



	cout << "\n";

	return 0;
}

