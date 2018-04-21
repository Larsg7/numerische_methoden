#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double m = 1;
const double l = 1;
const double g = 9.81;
const double h = 0.01;

const double t_max = 10;

double force(double phi)
{
    return m * g * std::sin(phi);
}

vector<double> get_f(double y[2], double t)
{
    vector<double> result = {y[1], force(y[0])};
    return result;
}

vector<double> euler(double phi_start)
{
    double phi = phi_start;
    double y[] = {phi, 0};
    vector<double> results(t_max / h);

    for (double t = 0; t < t_max; t += h)
    {
        results[t / h] = y[0];
        vector<double> f = get_f(y, t);
        y[0] += f[0] * h;
        y[1] += f[1] * h;
    }

    return results;
}

vector<double> get_k2(double y[2], double t)
{
    vector<double> k1 = get_f(y, t);
    y[0] += 0.5 * k1[0] * h;
    y[1] += 0.5 * k1[1] * h;
    vector<double> result = get_f(y, t + 0.5 * h);
    result[0] += result[0] * h;
    result[1] += result[1] * h;
    return result;
}

vector<double> get_k3(double y[2], double t)
{
    vector<double> k2 = get_k2(y, t);
    y[0] += 0.5 * k2[0] * h;
    y[1] += 0.5 * k2[1] * h;
    vector<double> result = get_f(y, t + 0.5 * h);
    return result;
}

vector<double> get_k4(double y[2], double t)
{
    vector<double> k3 = get_k3(y, t);
    y[0] += k3[0] * h;
    y[1] += k3[1] * h;
    vector<double> result = get_f(y, t + h);
    return result;
}

vector<double> runge_4(double phi_start)
{
    double phi = phi_start;
    double y[] = {phi, 0};
    vector<double> results(t_max / h);

    for (double t = 0; t < t_max; t += h)
    {
        results[t / h] = y[0];
        vector<double> f = get_f(y, t);
        y[0] += f[0] * h;
        y[1] += f[1] * h;
    }

    return results;
}

int main()
{
    double phi_start = 1;

    vector<double> euler_results = euler(phi_start);
}

