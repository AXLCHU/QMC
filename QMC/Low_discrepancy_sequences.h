#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

double Halton_seq(int index, int base) { // index & spatial dim 
	double f = 1, r = 0;
	while (index > 0) {
		f = f / base;
		r = r + f * (index % base);
		index = index / base;
	}
	return r;
}

// Random output Halton sequence element
double Halton() {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 mt1(seed);
	return Halton_seq(mt1(), 3);//
}
#pragma once
