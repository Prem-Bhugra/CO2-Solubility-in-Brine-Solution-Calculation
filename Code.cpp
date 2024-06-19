/* CLL 121 Thermodynamics Coding Project*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
using namespace std;
double A_f(double P, double Pc, double T, double Tc, double R)
{ // Calculating the dimensionless parameter A
    double omega = 0.224;
    double m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    double a = 0.457236 * (pow(R * Tc, 2) / Pc) * pow(1 + m * (1 - sqrt(T / Tc)), 2);
    return (a * P) / pow(R * T, 2);
}

double B_f(double P, double Pc, double T, double Tc, double R)
{ // Calculating the dimensionless parameter B
    double b = 0.077796 * R * Tc / Pc;
    return (b * P) / (R * T);
}

pair<double, double> root(double r, double s)
{
    double temp = (r * r) + (4 * s);
    double root1 = (r + sqrt(temp)) / 2;
    root1 = (root1 < 0 || temp < 0) ? 0 : root1;
    double root2 = (r - sqrt(temp)) / 2;
    root2 = (root2 < 0 || temp < 0) ? 0 : root2;
    pair<double, double> ans;
    ans.first = root1;
    ans.second = root2;
    return ans;
}

pair<double, double> Bairstow(int n, double &r, double &s, vector<double> &A, vector<double> &B)
{
    B[n] = A[n];
    B[n - 1] = A[n - 1] + (r * B[n]);
    for (int i = n - 2; i >= 0; i--)
    {
        B[i] = A[i] + (r * B[i + 1]) + (s * B[i + 2]);
    }
    double *C = new double[n];
    C[n - 1] = B[n];
    C[n - 2] = B[n - 1] + (r * C[n - 1]);
    for (int i = n - 3; i >= 0; i--)
    {
        C[i] = B[i + 1] + (r * C[i + 1]) + (s * C[i + 2]);
    }
    double delta_s = ((-1 * B[0]) + ((C[0] / C[1]) * B[1])) / (C[1] - ((C[0] / C[1]) * C[2]));
    double delta_r = ((-1 * B[1]) - (C[2] * delta_s)) / C[1];
    r = r + delta_r;
    s = s + delta_s;
    pair<double, double> percentage_errors;
    percentage_errors.first = abs((delta_r / r) * 100);
    percentage_errors.second = abs((delta_s / s) * 100);
    return percentage_errors;
}

vector<double> Solver(double a, double b)
{ // Solving the cubic Peng-Robinson Equation of state using Bairstow's method to find multiple roots
    int degree = 3;
    double r = -1, s = -1;
    vector<double> coefficients(degree + 1);
    coefficients[0] = -1 * (a * b - b * b - b * b * b);
    coefficients[1] = a - 2 * b - 3 * b * b;
    coefficients[2] = -1 * (1 - b);
    coefficients[3] = 1;
    vector<double> B(degree + 1);
    for (int i = 0; i < degree + 1; i++)
    {
        B[i] = 0;
    }
    vector<double> roots(degree);
    int k = 0;
    pair<double, double> errors;
    errors = Bairstow(degree, r, s, coefficients, B);
    while (errors.first > 1e-5 || errors.second > 1e-5)
    {
        errors = Bairstow(degree, r, s, coefficients, B);
    }
    pair<double, double> ans = root(r, s);
    roots[k++] = ans.first;
    roots[k++] = ans.second;
    degree = degree - 2;
    int x = 0;
    while (x < degree + 1)
    {
        coefficients[x] = B[x + 2];
        x++;
    }
    while (x < degree + 3)
    {
        coefficients[x] = 0;
        x++;
    }
    roots[k++] = (-1 * coefficients[0]) / coefficients[1];
    roots[k - 1] = (roots[k - 1] > 0) ? roots[k - 1] : 0;
    return roots;
}

double Z_f(vector<double> roots, double a, double b)
{ // Finding the correct value of compressibility factor from the multiply roots obtained
    double Z_l = min(roots[0], min(roots[1], roots[2]));
    double Z_g = max(roots[0], max(roots[1], roots[2]));
    double sigma = 1 + sqrt(2);
    double epsilon = 1 - sqrt(2);
    double difference = Z_g - Z_l + log((Z_l - b) / (Z_g - b)) - (a / (b * (epsilon - sigma))) * log(((Z_l + sigma * b) / (Z_l + epsilon * b)) * ((Z_g + epsilon * b) / (Z_g + sigma * b)));
    return (difference > 0) ? Z_l : Z_g;
}

double phi_f(double Z, double A, double B)
{ // Calculating the fugacity coefficicent for CO2
    double epsilon = 1 - sqrt(2);
    double sigma = 1 + sqrt(2);
    return exp(Z - 1 - log(Z - B) - (A / (B * (epsilon - sigma))) * log((Z + epsilon * B) / (Z + sigma * B)));
}

double gamma_f(double T, double P, double C)
{ // Calculating the activity coefficient for dissolved CO2
    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;
    c1 = -0.0652869;
    c2 = 1.6790636e-4;
    c3 = 40.838951;
    c4 = 0;
    c5 = 0;
    c6 = -3.9266518e-2;
    c7 = 0;
    c8 = 2.1157167e-2;
    c9 = 6.5486487e-6;
    c10 = 0;
    double lambda = c1 + c2 * T + (c3 / T) + c4 * P + (c5 / P) + c6 * (P / T) + c7 * (T / pow(P, 2)) + c8 * P / (630 - T) + c9 * T * log(P) + c10 * (P / pow(T, 2));
    c1 = -1.144624e-2;
    c2 = 2.8274958e-5;
    c3 = 0;
    c4 = 0;
    c5 = 0;
    c6 = 1.3980876e-2;
    c7 = 0;
    c8 = -1.4349005e-2;
    c9 = 0;
    c10 = 0;
    double epsilon = c1 + c2 * T + (c3 / T) + c4 * P + (c5 / P) + c6 * (P / T) + c7 * (T / pow(P, 2)) + c8 * P / (630 - T) + c9 * T * log(P) + c10 * (P / pow(T, 2));
    return exp(2 * C * lambda + 2 * C * C * epsilon);
}

double rho_f(double P, double T)
{ // Calculating the density of pure water
    T -= 273.15;
    double V, B, A1, A2;
    V = (1 + 18.1597 * 0.001 * T) / (0.9998 + 18.2249 * 0.001 * T - 7.9222 * (1e-6) * T * T - 55.4485 * (1e-9) * pow(T, 3) + 149.7562 * (1e-12) * pow(T, 4) - 393.2952 * (1e-15) * pow(T, 5));
    B = 19654.32 + 147.037 * T - 2.2155 * T * T + 1.0478 * 0.01 * pow(T, 3) - 2.2789 * (1e-5) * pow(T, 4);
    A1 = 3.2891 - 2.391 * (1e-3) * T + 2.8446 * (1e-4) * T * T - 2.82 * (1e-6) * pow(T, 3) + 8.477 * (1e-9) * pow(T, 4);
    A2 = 6.245 * (1e-5) - 3.913 * (1e-6) * T - 3.499 * (1e-8) * T * T + 7.942 * (1e-10) * pow(T, 3) - 3.299 * (1e-12) * pow(T, 4);
    double reciprocal = V - (V * P / (B + A1 * P + A2 * P * P));
    return 1 / reciprocal;
}

double fH2O_f(double P, double Pc, double T, double Tc, double R, double rho)
{ // Calculating the fugacity of H2O
    double a1, a2, a3, a4, a5, a6;
    a1 = -7.8595178;
    a2 = 1.8440825;
    a3 = -11.786649;
    a4 = 22.680741;
    a5 = -15.9618719;
    a6 = 1.8012250;
    double Tr = T / Tc;
    double Ps = Pc * exp((1 / Tr) * (a1 * (1 - Tr) + a2 * pow(1 - Tr, 1.5) + a3 * pow(1 - Tr, 3) + a4 * pow(1 - Tr, 3.5) + a5 * pow(1 - Tr, 4) + a6 * pow(1 - Tr, 7.5)));
    return Ps * exp(18.0152 * (P - Ps) / (rho * R * T));
}

double h_f(double fH2O, double rho, double T, double R)
{ // Calculating the Henry's constant for H2O
    double M = 18.01528;
    double eta = -0.114535;
    double delta_B = -5.279063 + 6.187967 * sqrt(1000 / T);
    return exp((1 - eta) * log(fH2O) + eta * log(R * T * rho / M) + 2 * rho * delta_B);
}

double K_f(double h, double gamma, double phi, double P)
{ // Calculating phase equilibrium constant
    return (h * gamma) / (P * phi);
}

pair<double, double> yH2O_f_and_K_h2o(double K, double P, double T, double fH2O, double R)
{ // Calculating mole fraction of H2O in gas phase
    T -= 273.15;
    double KoH2O = pow(10, -2.209 + 3.097 * 0.01 * T - 1.098 * (1e-4) * T * T + 2.048 * (1e-7) * T * T * T);
    T += 273.15;
    double KH2O = (KoH2O / (fH2O * P)) * exp((P - 1) * 18.18 / (R * T));
    return {(1 - (1 / K)) / ((1 / KH2O) - (1 / K)), KH2O};
}

double y_n_f(double yH2O)
{ // Calculating normalized mole fraction of CO2 in gas phase
    return 1 / (1 + yH2O);
}

double x_f(double y_n, double K)
{ // Calculating mole fraction of CO2 in liquid phase
    return y_n / K;
}

int main()
{
    ifstream f;
    ofstream s;
    f.open("Input.txt");
    s.open("Output.txt");

    int t;
    f >> t; // Number of test cases

    s << setw(15) << left << "Temperature" << setw(12) << "Pressure" << setw(17) << "NaCl Molality" << setw(32) << "Fugacity Coefficient for CO2" << setw(28) << "Henry's Constant for CO2" << setw(40) << "Phase Equilibrium Constant for H2O" << setw(37) << "Mole Fraction of H2O in Gas Phase" << setw(54) << "Mole fraction of CO2 in liquid phase (Calculated)" << endl;

    while (t--)
    {
        double T;
        f >> T; // Kelvin
        double P;
        f >> P; // bar
        double C;
        f >> C; // mole/KgH2O

        double Tc = 31.05; // degree Celsius
        double Pc = 73.83; // bar
        double R = 83.14;  //(bar*cm^3)/(mol*K)

        Tc += 273.15;

        double A = A_f(P, Pc, T, Tc, R);
        double B = B_f(P, Pc, T, Tc, R);
        vector<double> roots = Solver(A, B);
        double Z = Z_f(roots, A, B);
        double phi = phi_f(Z, A, B);
        double gamma = gamma_f(T, P, C);
        double rho = rho_f(P, T);
        double fH2O = fH2O_f(P, 220.55, T, 647.1, R, rho);
        double h = h_f(fH2O, rho, T, R);
        double K = K_f(h, gamma, phi, P);
        pair<double, double> temp = yH2O_f_and_K_h2o(K, P, T, fH2O, R);
        double yH2O = temp.first;
        double Kh2o = temp.second;
        double y_n = y_n_f(yH2O);
        double x = x_f(y_n, K);

        s << left <<setw(15) << T << setw(12) << P << setw(17) << C << setw(32) << phi << setw(28) << h << setw(40) << Kh2o << setw(37) << yH2O << setw(54) << setprecision(5) << x << endl;
    }
}

// Reference:-
/*
Shabani, B., & Vilcáez, J. (2017).
Prediction of CO2-CH4-H2S-N2 gas mixtures
solubility in brine using a non-iterative
fugacity-activity model relevant to CO2-MEOR.
Journal of Petroleum Science and Engineering,
150, 162-179.
https://doi.org/10.1016/j.petrol.2016.12.012
*/