#include <stdio.h>
#include <stdbool.h>
#include "utvpiInterpreter.h"

bool parseFileDelayedSysGen(FILE * constraintFile,
    void * object,
    void (* initializeSystem)(void * object, int n, Parser * parser),
    void (* addConstraint)(void * object, Constraint * constraint, Parser * parser));
