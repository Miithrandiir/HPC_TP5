#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#define NB_POINTS 1000000000
#define REAL_PI 3.141592653589793238462643
#define PRECISION 10

void seq_monte_carlo() {
    //std::vector<std::pair<int, int>> pts;
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    int I = 0;
    long double pi_precision = std::pow(10,-PRECISION);
    long double PI = 0.0;
    for (int i = 1; i <= NB_POINTS; ++i) {
        double x = dist(mt);
        double y = dist(mt);
        if (x * x + y * y < 1) {
            ++I;
        }

        PI = 4.0 * ((long double) I / (long double) i);
        if (std::abs(PI - REAL_PI) < pi_precision) {
            break;
        }

    }

    std::cout << "PI = " << std::scientific << std::setprecision(PRECISION) << PI << std::endl;
}

int main() {
    seq_monte_carlo();
    return 0;
}