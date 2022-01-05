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



// + Sobol, Faure, Hammersley


#pragma once
