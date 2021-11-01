#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include "Low_discrepancy_sequences.h"

double expo(double rand, double lambda) {
	double X = -log(rand) / lambda;
	return X;
}

// Acceptance-rejection method to generate SND

void gauss2(vector<double>& X) {

	double rand1 = 0; double rand2 = 0;
	double expo1 = 0; double expo2 = 0;
	bool cond = true;

	do {
		for (int i = 0; i < X.size(); i++) {

				rand1 = rand() / static_cast<double>(RAND_MAX);
				rand2 = rand() / static_cast<double>(RAND_MAX);
				expo1 = expo(rand1, 1);
				expo2 = expo(rand2, 1);
		
					if (expo2 > pow(expo1 - 1, 2)) { 
						if (rand1 <= 0.5) { X[i] = expo1; }
						else { X[i] = -expo1; }
						cond = true;
					}
					else { cond = false; i -= 1; }
		} 
	} while (cond == false);
}


// Generate gaussian distrib via Halton sequence

void gauss_low_discrepency(vector<double>& X, double& dim) {
	double rand1 = 0; double rand2 = 0;
	double expo1 = 0; double expo2 = 0;
	bool cond = true;

	unsigned seed = 1234;
	std::mt19937 mt1(seed);

	do {

		for (int i = 0; i < X.size(); i++) {

			int MT1 = mt1(); int MTT1 = abs(MT1);
			int MT2 = mt1(); int MTT2 = abs(MT2);

			rand1 = Halton_seq(MTT1, 3);
			rand2 = Halton_seq(MTT2, 3);

			expo1 = expo(rand1, 1); // size RN > size X
			expo2 = expo(rand2, 1);

			if (expo2 > pow(expo1 - 1, 2)) { // Accept
				if (rand1 <= 0.5) { X[i] = expo1; }
				else { X[i] = -expo1; }
				cond = true;
			}
			else { cond = false; i -= 1; }
		}
	} while (cond == false);
}



//Ziggurat


//Box-Mûller

double gaussian_box_muller() {
	double x = 0.0;
	double y = 0.0;
	double euclid_sq = 0.0;

	do {
		x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		euclid_sq = x * x + y * y;
	} while (euclid_sq >= 1.0);

	return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}

#pragma once
