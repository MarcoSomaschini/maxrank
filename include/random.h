#ifndef _random_h
#define _random_h

#include <time.h>

const unsigned long maxshort = 65536L;
const unsigned long multiplier = 1194211693L;
const unsigned long adder = 12345L;

class randomnumber {
public:
    unsigned long randseed;

public:
    randomnumber(unsigned long s = 0);

    unsigned short random(unsigned long n);

    double frandom(void);
};

randomnumber::randomnumber(unsigned long s) {
    if (s == 0)
        randseed = time(0);
    else
        randseed = s;
}

unsigned short randomnumber::random(unsigned long n) {
    randseed = multiplier * randseed + adder;
    return (unsigned short) ((randseed >> 16) % n);
}

double randomnumber::frandom(void) {
    double res;

    do {
        res = random(maxshort) / double(maxshort);
    } while (res == 0);

    return res;
}

#endif