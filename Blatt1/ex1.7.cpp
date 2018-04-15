#include <limits>
#include <iostream>

using namespace std;

template<typename T>
T smallestNumber()
{
    T last = 1;
    T next = 1;

    while (next != 0)
    {
        last = next;
        next = last / 2;
        // std::cout << next << std::endl;
    }
    return last;
}

int main()
{
    cout << "Numeric denormalized min float: " << numeric_limits<float>::denorm_min() << endl;
    cout << "Calculated min float: " << smallestNumber<float>() << endl;

    cout << endl;

    cout << "Numeric denormalized min double: " << numeric_limits<double>::denorm_min() << endl;
    cout << "Calculated min double: " << smallestNumber<double>() << endl;
}