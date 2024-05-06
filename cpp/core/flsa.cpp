#include <vector>
#include "flsa.hpp"

void flsa(const int n, const std::vector<double>& y, const double lam, const std::vector<int>& c, std::vector<double>& beta) {
    /*
    Inputs:
    - n: Number of data points
    - y: Vector of data
    - lam: Regularization parameter
    - c: Vector of 0 and 1, indicating overlaps of the data point
    - beta: Vector of smoothing values
    */

    // Trivial cases
    if (n == 0) return;
    if (n == 1 || lam == 0) {
        beta = y;
        return;
    }

    /*
    Initialize variables
    x: Knots
    slope_diff: Differences of slopes of the adjacent intervals of the piecewise linear derivative function
    intercept_diff: Differences of intercepts of the adjacent intervals of the piecewise linear derivative function
    */
    std::vector<double> x(2 * n);
    std::vector<double> slope_diff(2 * n);
    std::vector<double> intercept_diff(2 * n);

    /*
    New knots
    */
    std::vector<double> km(n - 1);
    std::vector<double> kp(n - 1);

    int i; //Index of data

    int lmost, rmost;// Index of leftmost and rightmost knots
    int lborder, rborder; // Index of knots that are on the border of the interval, where the value of the derivative is between -lam and lam
    double end_slope, first_intercept, last_intercept; // slope and intercept of the first and last linear function
    double slope_lborder, intercept_lborder, slope_rborder, intercept_rborder;

    // Forward pass
    // First iteration
    if (c[0]) { // When overlapping
        lmost = n;
        rmost = n - 1;
        end_slope = 2;
        first_intercept = -y[0] - y[1];
        last_intercept = first_intercept;
    } else {
        km[0] = -lam + y[0];
        kp[0] = lam + y[0];

        lmost = n - 1;
        rmost = n;

        x[lmost] = km[0];
        x[rmost] = kp[0];
        slope_diff[lmost] = 1;
        intercept_diff[lmost] = -y[0] + lam;
        slope_diff[rmost] = -1;
        intercept_diff[rmost] = y[0] + lam;

        end_slope = 1;
        first_intercept = -y[1] - lam;
        last_intercept = -y[1] + lam;
    }

    // Iteration through the data
    for (i = 1; i < n - 1; i++) {
        if(c[i]){ // When overlapping
        end_slope++;
        first_intercept -= y[i + 1];
        last_intercept -= y[i + 1];
        } else {
        // Find the leftmost knot that is on the left border of the interval between -lam and lam
        slope_lborder = end_slope;
        intercept_lborder = first_intercept;
        for (lborder = lmost; lborder <= rmost; lborder++) {
            if (slope_lborder * x[lborder] + intercept_lborder > -lam) break;
            slope_lborder += slope_diff[lborder];
            intercept_lborder += intercept_diff[lborder];
        }

        // Find the rightmost knot that is on the right border of the interval between -lam and lam
        slope_rborder = end_slope;
        intercept_rborder = last_intercept;
        for (rborder = rmost; rborder >= lborder; rborder--) {
            if (slope_rborder * x[rborder] + intercept_rborder < lam) break;
            slope_rborder -= slope_diff[rborder];
            intercept_rborder -= intercept_diff[rborder];
        }

        // Compute the knots
        km[i] = (-lam - intercept_lborder) / slope_lborder;
        lmost = lborder - 1;
        x[lmost] = km[i];

        kp[i] = (lam - intercept_rborder) / slope_rborder;
        rmost = rborder + 1;
        x[rmost] = kp[i];

        // The first and last differences of the slopes and intercepts are changed
        slope_diff[lmost] = slope_lborder;
        intercept_diff[lmost] = intercept_lborder + lam;
        slope_diff[rmost] = -slope_rborder;
        intercept_diff[rmost] = -intercept_rborder + lam;

        end_slope = 1;
        first_intercept = -y[i + 1] - lam;
        last_intercept = -y[i + 1] + lam;
        }
    }

    // Last knot
    slope_lborder = end_slope;
    intercept_lborder = first_intercept;
    for (lborder = lmost; lborder <= rmost; lborder++) {
        if (slope_lborder * x[lborder] + intercept_lborder > 0) break;
        slope_lborder += slope_diff[lborder];
        intercept_lborder += intercept_diff[lborder];
    }
    beta[n - 1] = -intercept_lborder / slope_lborder;

    // Backward pass
    for (i = n - 2; i >= 0; i--) {
        if (c[i]){
            beta[i] = beta[i + 1];
        } else {
            beta[i] = std::max(km[i], std::min(kp[i], beta[i + 1]));
        }
    }
}