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

using namespace std;

double Halton_seq(int index, int base) {
	double f = 1, r = 0;
	while (index > 0) {
		f = f / base;
		r = r + f * (index % base);
		index = index / base;
	}
	return r;
}


double expo(double rand, double lambda) {
	double X;
	X = -log(rand) / lambda;
	return X;
}

// Rejection method to generate SND
double gauss(double rand) {
	double expo1 = expo(rand, 1);
	double expo2 = expo(rand, 1);
	double Z = 0;

	if (expo2 > pow(expo1 - 1, 2)) {
		if (rand <= 0.5) {
			Z = expo1;
		}
		else Z = -expo1;
	}
	return Z;
}



int main() {

	double rand = Halton_seq(10, 4);

/*	cout << rand;
	cout << "\n" << expo(rand,1);
	cout << "\n" << gauss(rand);*/

	double alpha = 0.95;
	double dim = 10000;

	std::vector<double> X(dim,0);
	std::vector<double> Y(dim, 0);
	std::vector<double> Z(dim, 0);
	std::vector<double> VaR(dim, 10);
	std::vector<double> CVaR(dim, 10);

	std::normal_distribution<double> dist1{ 0, 1 };

	for (int i = 0; i < dim; i++) {
		double rand = Halton_seq(i+1, 2); //
		X[i] = gauss(rand);

		//	std::random_device seed;
		//	std::mt19937 gen{ seed() };

		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		default_random_engine gen(seed);
		Y[i] = dist1(gen);

		Z[i] = gaussian_box_muller();
	}


	VaR_CVaR(CVaR, VaR, Z, alpha);

	double variance = 0;
	double sum = accumulate(VaR.begin(), VaR.end(), 0);
	double mean = sum / VaR.size();

	for (int n = 0; n < dim; n++) {
		variance+=(VaR[n] - mean) * (VaR[n] - mean);
	}
	variance /= dim;
	double std = sqrt(variance);

	cout << "\nStandand Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1]; cout << " & Variance de VaR avec X = " << std;
	cout << "\nStandand Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];

	VaR_CVaR(CVaR, VaR, Y, alpha);

	for (int n = 0; n < dim; n++) {
		variance += (VaR[n] - mean) * (VaR[n] - mean);
	}
	variance /= dim;
	std = sqrt(variance);

	cout << "\n\nVaR = " << VaR[dim - 1]; cout << " & Variance de VaR avec Y = " << std;
	cout << "\nCVaR = " << CVaR[dim - 1];
	cout << "\n";

	double S0 = 100;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	unsigned num_intervals = 250;

	cout << "\nPrice Put option = " << BSM_Put(S0,r,v,T,K);
	cout << "\n";

	vector<double> spot_prices(num_intervals, S0);
	vector<double> loss_fct(dim);

	for (int n = 0; n < dim; n++) {
		GBM_paths(spot_prices, r, v, T);
		loss_fct[n] = std::max(K - spot_prices[num_intervals - 1],0.0) - exp(r*T) * BSM_Put(S0, r, v, T, K);
	}

	VaR_CVaR(CVaR, VaR, loss_fct, alpha);

	cout << "\nVaR at "<< alpha*100 << "% for Put option = " << VaR[dim-1];
	cout << "\nCVaR at " << alpha * 100 << "% for Put option  = " << CVaR[dim - 1];
	cout << "\n";


	return 0;
}

