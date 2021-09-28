#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

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
