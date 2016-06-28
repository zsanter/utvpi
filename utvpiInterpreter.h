#ifndef _UTVPI_INTERPRETER_H
#define _UTVPI_INTERPRETER_H

#include <stdio.h>
#include <stdbool.h>
#include "constraint.h"

typedef struct Parser Parser;

bool parseFile(FILE * constraintFile,
    void * object,
    void (* initializeSystem)(void * object, int n, Parser * parser), 
    void (* addConstraint)(void * object, Constraint * constraint, Parser * parser));

void parseError(Parser * parser, const char * message);

#endif
