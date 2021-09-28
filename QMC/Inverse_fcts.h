#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

double expo(double rand, double lambda) {
	double X = -log(rand) / lambda;
	return X;
}

// Acceptance-rejection method to generate SND
double gauss(double rand) {
	double expo1 = expo(rand, 1);
	double expo2 = expo(rand, 1); //
	double Z = 0;

	if (expo2 > pow(expo1 - 1, 2)) {
		if (rand <= 0.5) {
			Z = expo1;
		}
		else Z = -expo1;
	}
/*	else { Z = ; } */
	return Z;
}

/*
double gauss2(double rand1, double rand2) {

//	double rand1 = rand() / static_cast<double>(RAND_MAX) - 1;
//	double rand2 = rand() / static_cast<double>(RAND_MAX) - 1;


	double targetY = 0.4 * exp(-rand1 * rand1 * 0.5);

	if (rand2 < targetY) {

		return rand1;
	}
	else { 
		//gauss2(); 
		return 0;
	}

}
*/


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

/*
double gaussian_box_muller_v2(const double rand) {
	double x = 0.0;
	double y = 0.0;
	double euclid_sq = 0.0;
	do {
		x = 2.0 * rand;
		y = 2.0 * rand;
		euclid_sq = x * x + y * y;
	} while (euclid_sq >= 1.0);

	return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}
*/

#pragma once
