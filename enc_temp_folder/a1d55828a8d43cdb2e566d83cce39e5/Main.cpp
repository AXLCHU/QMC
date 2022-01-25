#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <climits>
#include <list>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>
#include <chrono>
#include "VaR_CVaR.h"
#include "BSM.hpp"
#include "Path_generation.hpp"
#include "Low_discrepancy_sequences.h"
#include "Sampling.h"
#include "MeanVariance.h"
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/sobol.hpp>
#include <boost/random/detail/sobol_table.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>

using namespace std;


int main() {

	double alpha = 0.99;
	double dim = 10000;

	 vector<double> MT(dim, 0);	vector<double> RN(dim, 0);
	 vector<double> gausss(dim, 0);
	 vector<double> halton1(dim, 0);
	 vector<double> halton2(dim, 0);
	 vector<double> sobol(dim, 0);

	 vector<double> X(dim,0);  vector<double> X1(dim, 0);
	 vector<double> Y(dim, 0);  vector<double> Y1(dim, 0); vector<double> Y2(dim, 0);
	 vector<double> Z(dim, 0);
	 vector<double> Z1(dim, 0);
	 vector<double> VaR(dim, 0);
	 vector<double> CVaR(dim, 0);
	 vector<double> Loss(dim, 0);


	 /*
	 static const std::size_t dimension = dim;

//	 typedef boost::variate_generator<boost::random::sobol&, boost::uniform_01<double> > quasi_random_gen_t;	 // Create a generator
	 typedef boost::variate_generator<boost::random::sobol&, boost::normal_distribution<double> > quasi_random_gen_t;

	 boost::random::sobol engine(dimension); // Initialize the engine to draw randomness out of thin air
	 
//	 quasi_random_gen_t gen(engine, boost::uniform_01<double>()); // Glue the engine and the distribution together
	 quasi_random_gen_t gen(engine, boost::random::normal_distribution<double>());

	 std::vector<double> sample(dimension);

	 std::generate(sample.begin(), sample.end(), gen);	 // At this point you can use std::generate, generate member f-n, etc.
	 engine.generate(sample.begin(), sample.end());

//	 boost::normal_distribution<double> ND(0.0, 1.0);
//	 boost::random::sobol QRNG(2);
//	 boost::random::sobol_engine QRNG(dimension);
	 
//	 boost::variate_generator<boost::random::sobol&, boost::random::normal_distribution<>> var_SND(QRNG, ND);

	 for (int i = 0; i < dim; i++) {
		 sobol[i] = var_SND();
	 }
	 
	 VaR_CVaR(CVaR, VaR, sample, alpha);
	 cout << "\n\SOBOL Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1];
	 cout << "\nSOBOL Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];
*/


	//std::uniform_real_distribution<double> unif(0,1);
	boost::math::normal dist(0, 1);

	double q = quantile(dist, alpha);
	double qq = (1 / (1 - alpha)) * ((1 / sqrt(2 * M_PI)) * exp(-q * q * 0.5));

	cout << "\nStandard Gaussian VaR at " << alpha * 100 << "% = " << q;
	cout << "\nStandard Gaussian CVaR at " << alpha * 100 << "% = " << qq;


	gauss2(X1);
	VaR_CVaR(CVaR, VaR, X1, alpha);

	double std_VaR = StdDev(VaR);
	double std_CVaR = StdDev(CVaR);

	cout << "\n\nPseudo-random case :";
	cout << "\n\nRejection method:";
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1]; cout << " & std dev de VaR = " << std_VaR * 100 << "%" ;
	cout << "\nStandard Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1]; cout << " & std dev de CVaR = " << std_CVaR * 100 << "%";
	cout << "\n\nErrror for VaR = " << q - VaR[dim - 1];
	cout << "\nErrror for CVaR = " << qq - CVaR[dim - 1];

	default_random_engine generator(123); // linear congruential engine
	//mersenne_twister_engine mt_engine(seed);
	std::normal_distribution<double> normal(0.0, 1.0);
/*
	for (int i = 0; i < dim; i++) {
		
		Z[i] = gaussian_box_muller();

		halton1[i] = Halton_seq(i + 1, 2);
		halton2[i] = Halton_seq(i + 1, 3);

		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
		std::mt19937 mt1(seed);

		MT[i] = mt1();
		RN[i] = mt1() / static_cast<double>(mt1.max()); 

		Y[i] = normal(generator);
//		Y[i] = normal(mt1);

	}
*/


//	VaR_CVaR_final_procedure(CVaR, VaR, VaR2, theta2, mu2, X1, alpha);


	//VaR_CVaR(CVaR, VaR, Z, alpha);
	std_VaR = StdDev(VaR);
	std_CVaR = StdDev(CVaR);

	cout << "\n\nBox-Muller: ";
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR = " << std_VaR * 100 << "%";
	cout << "\nStandard Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1] << " & std dev de CVaR = " << std_CVaR * 100 << "%";
	cout << "\n\nErrror for VaR = " << q - VaR[dim - 1];
	cout << "\nErrror for CVaR = " << qq - CVaR[dim - 1];

	/* MERSENNE-TWISTER
	boost::normal_distribution<double> ND(0.0, 1.0);
	mt19937 PRNG;
	boost::variate_generator<mt19937&, boost::random::normal_distribution<>> var_SND(PRNG, ND);
	for (int i = 0; i < dim; i++) { MT[i] = var_SND(); } */

	VaR_CVaR(CVaR, VaR, Y, alpha);
	std_VaR = StdDev(VaR);
	std_CVaR = StdDev(CVaR);

	cout << "\n\nMT19937: ";
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR = " << std_VaR * 100 << "%";
	cout << "\nStandard Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1] << " & std dev de CVaR = " << std_CVaR * 100 << "%";
	cout << "\n\nErrror for VaR = " << q - VaR[dim - 1];
	cout << "\nErrror for CVaR = " << qq - CVaR[dim - 1];

//QRNG
/*
	gauss_low_discrepency(X,dim);
	VaR_CVaR(CVaR, VaR, X, alpha);
	std_VaR = StdDev(VaR);
	std_CVaR = StdDev(CVaR);

	cout << "\n\n\nQuasi-random case:";
	cout << "\n\nWith Halton sequence:";
	cout << "\n\nRejection method:";
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR = " << std_VaR * 100 << "%";
	cout << "\nStandard Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1] << " & std dev de CVaR avec Halton = " << std_CVaR * 100 << "%";
	cout << "\n\nErrror for VaR = " << q - VaR[dim - 1];
	cout << "\nErrror for CVaR = " << qq - CVaR[dim - 1];

	BM_gauss_low_discrepency(X, dim);
	VaR_CVaR(CVaR, VaR, X, alpha);

	cout << "\n\nBox-Muller:";
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1];
	cout << "\nStandard Gaussian CVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];  */
	cout << "\n\n";



	/*
	VaR_CVaR(CVaR, VaR, X, alpha);
	cout << "\n\nStandard Gaussian VaR at " << alpha * 100 << "% = " << VaR[dim - 1] << " & std dev de VaR avec MT = " << std * 100 << "%";
	cout << "\nCVaR at " << alpha * 100 << "% = " << CVaR[dim - 1];
	*/



	double S0 = 100;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	unsigned num_intervals = 250;

	cout << "\n\nFor S0 = " << S0 << ", K = " << K << ", r = " << r << " and v = " << v << " : ";
	cout << "\n\nPut option price = " << BSM_Put(S0,r,v,T,K);
	cout << "\nCall option price = " << BSM_Call(S0, r, v, T, K);

	vector<double> spot_prices(num_intervals, S0);
	double P0 = BSM_Put(S0, r, v, T, K);
	std::vector<double> loss_fct(dim);

	generate_Loss(loss_fct, spot_prices, r, v, T, K);

	double dim2 = dim / 100;
	vector<double> VaR2(dim2, 0); ///
	vector<double> theta2(dim2, 0);
	vector<double> mu2(dim2, 0);
	vector<double> VaaR(dim);
	vector<double> CVaaR(dim, 0); 

	VaR_CVaR_final_procedure(CVaaR, VaaR, VaR2, theta2, mu2, loss_fct, alpha);

	cout << "\n\nVaR at "<< alpha*100 << "% for Put option = " << VaaR[dim-1];
	cout << "\nCVaR at " << alpha * 100 << "% for Put option  = " << CVaaR[dim - 1];
/*
	for (int n = 0; n < dim; n++) {
		GBM_paths(spot_prices, r, v, T);
		loss_fct[n] = std::max(spot_prices[num_intervals - 1] - K, 0.0) - exp(r * T) * BSM_Call(S0, r, v, T, K); ///
	}

	VaR_CVaR(CVaR, VaR, loss_fct, alpha);

	cout << "\n\nVaR at " << alpha * 100 << "% for Call option = " << VaR[dim - 1];
	cout << "\nCVaR at " << alpha * 100 << "% for Call option  = " << CVaR[dim - 1];
*/
	cout << "\n";

	return 0;
}

