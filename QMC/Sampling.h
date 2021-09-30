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
		
					if (expo2 > pow(expo1 - 1, 2)) { // Accept
						if (rand1 <= 0.5) { X[i] = expo1; }
						else { X[i] = -expo1; }
						cond = true;
					}
					else { cond = false; i -= 1; }
		} 
	} while (cond == false);
}



/*
void gauss_low_discrepency(vector<double>& X, vector<double>& rand1, vector<double>& rand2) {

	double expo1 = 0; double expo2 = 0;
	bool cond = true;

		for (int i = 0; i < X.size(); i++) {
			do {
				rand1[i] = Halton_seq(i, 2);
				rand2[i] = Halton_seq(i, 3);

				for (int j = i; j < rand1.size(); j++) { //

					expo1 = expo(rand1[j], 1); // size RN > size X
					expo2 = expo(rand2[j], 1);

					if (expo2 > pow(expo1 - 1, 2)) {
						if (rand1[j] <= 0.5) { X[i] = expo1; }
						else { X[i] = -expo1; }
						cond = true;
					} 
					else { cond = false; j += 1; } //
				}
			} while (cond == false);
		}
}*/


void gauss_low_discrepency(vector<double>& X, double& dim) {

	vector<double> rand1(dim,0); vector<double> rand2(dim, 0);
	double expo1 = 0; double expo2 = 0;
	bool cond = true;

	for (int i = 0; i < X.size(); i++) {

		rand1[i] = Halton(); //
		rand2[i] = Halton(); 

		expo1 = expo(rand1[i], 1); // size RN > size X
		expo2 = expo(rand2[i], 1);

		if (expo2 > pow(expo1 - 1, 2)) {
			if (rand1[i] <= 0.5) { 
				X[i] = expo1; 
			}
			else { 
				X[i] = -expo1; 
			}
			cond = true;
		}
		else { 
			cond = false; 
		}	// take next Halton nbr index j for same i

	}
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
