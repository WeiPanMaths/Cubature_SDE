#include <iostream>
#include <climits>
#include <string.h>
#include "Utils.h"

#define dlength sizeof(double)*CHAR_BIT

//template<typename InputIterator>
//double bitstring_to_double(InputIterator begin, InputIterator end)
//{
//    uint64_t x = 0;
//    for (; begin != end; ++begin)
//        x = (x << 1) + (*begin - 0);
//
//    double d;
//    memcpy(&d, &x, 8);
//    return d;
//}

int main()
{
    std::cout.precision(15);
    double v = 10.5;
    //   unsigned int bit_max = 52;

    int bit_rep[dlength] = { 0 };
    int lower[dlength] = { 0 };
    int upper[dlength] = { 0 };

    union {
        double value;
        char   array[sizeof(double)];
    };

    value = v;

    for (unsigned int i = 0; i < dlength; ++i) {
        int relativeToByte = i % CHAR_BIT;
        bool isBitSet = (array[sizeof(double) - 1 - i / CHAR_BIT] & (1 << (CHAR_BIT - relativeToByte - 1)))
            == (1 << (CHAR_BIT - relativeToByte - 1));
        // std::cout << (isBitSet ? "1" : "0");
        bit_rep[i] = (isBitSet ? 1 : 0);
        lower[i] = bit_rep[i];
        upper[i] = bit_rep[i];
    }
    for (unsigned int i = 0; i < dlength; ++i)
        std::cout << bit_rep[i];

    std::cout << std::endl << "original " << bitstring_to_double(std::begin(bit_rep), std::end(bit_rep))
        << std::endl;

    unsigned int bit = 35;

    //   for (unsigned int i = 0; i < bit_max - bit; ++i)
    for (unsigned int i = 0; i < bit; ++i)
    {//  		half_bit_rep[12 + i] = 1;
        lower[dlength - bit + i] = 0;
        upper[dlength - bit + i] = 1;
    }

    for (unsigned int i = 0; i < dlength; ++i) {
        std::cout << lower[i];
    }
    std::cout << std::endl;
    for (unsigned int i = 0; i < dlength; ++i) {
        std::cout << upper[i];
    }

    std::cout << std::endl << "lower: " << bitstring_to_double(std::begin(lower), std::end(lower)) << std::endl;
    std::cout << "upper: " << bitstring_to_double(std::begin(upper), std::end(upper));

    return 0;
}

