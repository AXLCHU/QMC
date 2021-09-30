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


double Mean(std::vector<double>& X) {
	double sum = std::accumulate(X.begin(), X.end(), X[0]);
	double mean = sum / X.size();
	return mean;
}

double StdDev(std::vector<double>& X) {
	double variance = 0;
	for (int n = 0; n < X.size(); n++) {
		variance += (X[n] - Mean(X))* (X[n] - Mean(X)) / X.size();
	}
	double std = sqrt(variance);
	return std;
}
