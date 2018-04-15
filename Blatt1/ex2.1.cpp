#include <iostream>
#include <cmath>
#include <fstream>
#include <boost/format.hpp>

unsigned long fib(unsigned n)
{
    unsigned long i = 1;
    unsigned long j = 1;
    unsigned long next = 0;
    if (n < 1)
    {
        return 0;
    }
    if (n == 1 || n == 2)
    {
        return 1;
    }
    for (unsigned k = 2; k < n; k++)
    {
        next = i + j;
        i = j;
        j = next;
    }
    return next;
}

template<typename T>
T phi()
{
    return (1+std::sqrt(5)) / 2;
}

template<typename T>
T delta(unsigned n)
{
    if (n < 2)
    {
        return 0;
    }
    return static_cast<T>(fib(n)) / static_cast<T>(fib(n-1)) - phi<T>();
}

int main()
{
    std::string float_file_path = "./float.dat";
    std::string double_file_path = "./double.dat";

    std::ofstream float_file (float_file_path, std::ios::out);
    std::ofstream double_file (double_file_path, std::ios::out);

    if (!float_file.is_open() || !double_file.is_open())
    {
        std::cerr << "Could not open files!\n";
        return 1;
    }

    const unsigned limit = 200;

    for (unsigned i = 0; i < limit; i++)
    {
        float_file << boost::format("%o\t%.6f\n") % i % delta<float>(i);
        double_file << boost::format("%o\t%.15f\n") % i % delta<double>(i);
    }
}