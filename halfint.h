#ifndef _HALFINT_H
#define _HALFINT_H

#include <stdbool.h>

typedef int half_int;

half_int intToHalfInt(int input);
int halfIntToInt(half_int input);
double halfIntToDouble(half_int input);
half_int intDivBy2ToHalfInt(int input);
bool halfIntIsIntegral(half_int input);
half_int halfIntCeil(half_int input);
half_int halfIntFloor(half_int input);

#endif
