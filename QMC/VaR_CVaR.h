#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std; 


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


void VaR_CVaR(vector<double>& CVaR, vector<double>& VaR, vector<double>& X, const double& alpha) {

    double gamma;
        for (int k = 1; k < VaR.size(); k++) {

            gamma = 1 / (double)k;
            
            if (X[k] > VaR[k - 1]) {
                VaR[k] = VaR[k - 1] + gamma * alpha / (1 - alpha);
                CVaR[k] = CVaR[k - 1] - gamma * (CVaR[k - 1] - VaR[k - 1] - (1 / (1 - alpha)) * (X[k] - VaR[k - 1]));
            }
            else {
                VaR[k] = VaR[k - 1] - gamma;
                CVaR[k] = CVaR[k - 1] - gamma * (CVaR[k - 1] - VaR[k - 1]);

            }
        }
}
