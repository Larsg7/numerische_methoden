#include <iostream>
#include <cmath>
#include <boost/format.hpp>

template <typename T>
T phi(bool negative = false)
{
    return negative ? (-1 - std::sqrt(5)) / 2 : (-1 + std::sqrt(5)) / 2;
}

template <typename T>
T iterative(unsigned n)
{
    T this_phi = phi<T>();
    T last_phi = 1;
    T next = 0;
    if (n == 0)
    {
        return 1;
    }
    if (n == 1)
    {
        return this_phi;
    }
    for (int i = 1; i < n; i++)
    {
        next = last_phi - this_phi;
        last_phi = this_phi;
        this_phi = next;
    }
    return next;
}

template <typename T>
T raise(unsigned n)
{
    return std::pow(phi<T>(), n);
}

int main()
{
    int limit = 20;

    for (int i = 0; i < limit; i++)
    {
        std::cout << boost::format("%o\t%.15f\t%.15f\t%.15f\t%.15f") % i % raise<float>(i) % raise<double>(i) % iterative<float>(i) % iterative<double>(i) << std::endl;
    }
}