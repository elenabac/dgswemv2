#include "basis_functions.h"

std::vector<double> jacobi_polynomial(int n, int a, int b, const std::vector<double>& x) { 
	size_t n_pts = x.size();

	std::vector<double> P(n_pts);
	
	if (n == 0) {
        for (size_t i = 0; i < n_pts; i++) { P[i] = 1; }
    }
    else if (n == 1) {
        for (size_t i = 0; i < n_pts; i++) { P[i] = (a - b + (a + b + 2)*x[i]) / 2; }
    }
    else {
        Array2D<double> J(3);
        J[0].reserve(n_pts);
        J[1].reserve(n_pts);
        J[2].reserve(n_pts);

        for (size_t i = 0; i < n_pts; i++) {
            J[0].push_back(1);
            J[1].push_back((a - b + (a + b + 2)*x[i]) / 2);
			J[2].push_back(0);
        }

        double a1;
        double a2;
        double a3;
        double a4;

        for (int i = 1; i < n; i++) {
            a1 = 2 * (i + 1)*(i + a + b + 1)*(2 * i + a + b);
            a2 = (2 * i + a + b + 1)*(a*a - b*b);
            a3 = (2 * i + a + b)*(2 * i + a + b + 1)*(2 * i + a + b + 2);
            a4 = 2 * (i + a)*(i + b)*(2 * i + a + b + 2);

            for (size_t j = 0; j < n_pts; j++) {
                J[(i + 1) % 3][j] = ((a2 + a3*x[j])*J[i % 3][j] - a4*J[(i - 1) % 3][j]) / a1;
            }
        }

        for (size_t i = 0; i < n_pts; i++) {
            P[i] = J[n % 3][i];
        }
    }

	return P;
}

std::vector<double> jacobi_polynomial_derivative(int n, int a, int b, const std::vector<double>& x) {
	size_t n_pts = x.size();
	
	std::vector<double> dP(n_pts);

	if (n == 0) {
        for (size_t i = 0; i < n_pts; i++) { dP[i] = 0; }
    }
    else if (n == 1) {
        for (size_t i = 0; i < n_pts; i++) { dP[i] = (a + b + n + 1) / 2.0; }
    }
    else {
        dP = jacobi_polynomial(n - 1, a + 1, b + 1, x);

        for (size_t i = 0; i < n_pts; i++) { dP[i] = dP[i] * (a + b + n + 1) / 2.0; }
    }

	return dP;
}