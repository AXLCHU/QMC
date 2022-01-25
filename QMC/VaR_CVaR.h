#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

void VaR_CVaR(vector<double>& CVaR, vector<double>& VaR, vector<double>& X, const double& alpha) {

    double gamma;
        for (int k = 1; k < VaR.size(); k++) {

            gamma = 1 / (double)pow(k, 2 / 3);
            
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

// Function to give stepwise constant confidence level
double alpha_n(int n, double alpha, int nb_sim) {
    double M1 = nb_sim / 3;
    if (n <= M1) {
        return 0.93;
    }
    else if (M1 < n <= int(2 * M1)) {
        return 0.96;
    }
    else if (int(2 * M1) < n <= nb_sim) {
        return alpha;
    }
}

double normalCDF(double value)
{
    return 0.5 * erfc(-value * sqrt(0.5));
}


void VaR_CVaR_final_procedure(vector<double>& CVaR, vector<double>& VaR, vector<double>& VaR2, vector<double>& theta2, vector<double>& mu2, vector<double>& X, const double& alpha) {

    double M = VaR2.size();
    //int M = int(dim) / 100; // int
    double gamma;

    for (int k = 1; k < M; k++) {

        gamma = 1 / (double)pow(k, 2 / 3);

        if (X[k] > VaR2[k - 1]) {
            VaR2[k] = VaR2[k - 1] + gamma * alpha_n(k, alpha, M) / (1 - alpha_n(k, alpha, M));
        }
        else { VaR2[k] = VaR2[k - 1] - gamma; }

        if (X[k] - theta2[k - 1] > VaR2[k - 1]) {
            theta2[k] = theta2[k - 1] - gamma * (2 * theta2[k - 1] - X[k]);
        }
        else { theta2[k] = theta2[k - 1]; }

        if (X[k] - mu2[k - 1] > VaR2[k - 1]) {
            mu2[k] = mu2[k - 1] - gamma * (1 / (1 + normalCDF(-mu2[k - 1]) + pow(VaR2[k - 1], 2))) * pow(X[k] - theta2[k - 1] - VaR2[k - 1], 2) * (2 * mu2[k - 1] - X[k]);
        }
        else { mu2[k] = mu2[k - 1]; }
    }

    double theta = theta2[M - 1];
    double mu = mu2[M - 1];
    VaR[0] = VaR2[M-1];

    for (int k = 1; k < X.size(); k++) {

        gamma = 1 / (double)pow(k, 2 / 3);

        if (X[k] >= VaR[k - 1]) {
            VaR[k] = VaR[k - 1] + gamma * alpha / (1 - alpha); // theta frozen
            CVaR[k] = CVaR[k - 1] - gamma * (CVaR[k - 1] - VaR[k - 1] - (1 / (1 - alpha)) * (X[k] - VaR[k - 1])); // mu frozen
        }
        else {
            VaR[k] = VaR[k - 1] - gamma;
            CVaR[k] = CVaR[k - 1] - gamma * (CVaR[k - 1] - VaR[k - 1]);
        }
    }
}
