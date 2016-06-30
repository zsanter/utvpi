/* 
 * utvpiInterpreter.h
 * 
 * utvpiInterpreter is used to read UTVPI input files, and transfer information read from those files into UTVPI system solvers.
 * The UTVPI files are human-readable text files, with support for C-style line and block comments. No semicolons are used. The 
 * newline character determines the end of a statement.
 *
 * Modified Backus-Naur Form:
 * <UTVPI_INPUT_FILE> ::= <VARIABLES_LINE><CONSTRAINT_LINE>*
 * <VARIABLES_LINE> ::= <INTEGER> variables<NEWLINE>
 * <CONSTRAINT_LINE> ::= {<SIGN>}<VARIABLE>{<SIGN><VARIABLE>}<LESS_THAN_OR_EQUAL><INTEGER><NEWLINE>
 *
 * The first statement in a UTVPI input file must be an integer followed by the keyword "variables". The integer represents the
 * number of variables within the constraint system. No variable may have an index higher than this value. Additionally, no 
 * variable may have an index less than 1.
 * 
 * All following statements must be constraint definitions. A constraint definition consists of one or two variables, a less-than-
 * or-equal sign, and an integral weight. The first variable may be implicitly or explicitly positive, or explicitly negative. The
 * second variable, if one is present, must be explicitly positive or negative. There may or may not be space between the sign and
 * the variable that it modifies. An integral weight may be implicitly or explicitly positive, or explicitly negative. There must 
 * be no space between a sign and an integer that it modifies. 
 *
 * A variable consists of an "x", either upper- or lower-case, optionally followed by an underscore, then immediately followed by
 * an unsigned integer representing the index of the variable.
 *
 * The following is all compliant with this interpreter:
 *
 *      5 variables
 *      x1+x2<=5
 *      -x1-x3<=-6
 *      -x_4 +X5 <= 8
 *      x2 + x5 <= -7
 *      + x4 - x2 <= 4
 *      - X_3 + x4 <= -2
 *      -x3 <= 4
 * 
 * As each constraint line is parsed, the Constraint struct defined in constraint.h is filled with information. The first variable
 * discovered provides sign[0] and index[0]. The second variable, if one exists, fills sign[1] and index[1]. If there is no second
 * variable, sign[1] == CONSTRAINT_NONE, and index[1] == 0. sign[1] == CONSTRAINT_NONE is the standard way to determine that the
 * Constraint struct describes an absolute constraint, one with only one variable. Obviously, weight holds the weight of the 
 * constraint.
 */

#ifndef _UTVPI_INTERPRETER_H
#define _UTVPI_INTERPRETER_H

#include <stdio.h>
#include <stdbool.h>
#include "constraint.h"

/*
 * struct Parser is defined within the source file. The typedef is here so that a Parser struct pointer may be passed into the 
 * calling system solver's initializeSystem() and addConstraint() functions, where it may be used to call parseError(). The 
 * contents of the Parser struct are opaque to the calling system solver.
 */
typedef struct Parser Parser;

/*
 * parseFile()
 *   A single function call to perform parsing of a UTVPI input file, facilitating the initialization and building of the calling
 *   system solver's representation of the system as the file is parsed
 * constraintFile - FILE pointer to a text file with the formatting specified above
 * object - void pointer, not used by this function, but passed back to the calling system solver as an argument to 
 *   initializeSystem() and addConstraint(). Can be used to store a reference to an already-declared overall graph system struct, 
 *   thereby avoiding the use of global variables
 * initializeSystem - pointer to a function within the calling system solver, used to initialize the solver's representation of 
 *   the system, given the number of variables within the system
 *   object - same void pointer originally passed into parseFile()
 *   n - number of variables present within the system, as stated at the beginning of the UTVPI input file
 *   parser - Parser struct pointer, necessary when calling parseError()
 * addConstraint - pointer to a function within the calling system solver, used to add the given constraint to the solver's 
 *   representation of the system
 *   object - same void pointer originally passed into parseFile()
 *   constraint - pointer to a Constraint struct, storing all information about a single constraint read from the UTVPI input file
 *   parser - Parser struct pointer, necessary when calling parseError()
 * returns true upon successful parse, false if any errors were encountered
 */
bool parseFile(FILE * constraintFile,
    void * object,
    void (* initializeSystem)(void * object, int n, Parser * parser), 
    void (* addConstraint)(void * object, Constraint * constraint, Parser * parser));

/*
 * parseError()
 *   In lieu of throwing an exception, not an option available to us in C, the calling system solver can call this function when a
 *   constraint is passed to it that it can't represent, or for any other reason. parseError() will handle this call the same way 
 *   as it would handle any error generated internally.
 *   parseError() will print the line and column numbers of the location within the input file from which the error originated, 
 *   followed by the text of the line of the input file, and then the text of the message given as input to the function, all to
 *   stderr.
 * parser - pointer to the Parser struct used for tracking information about the parse
 * message - message to be printed with the error output
 */
void parseError(Parser * parser, const char * message);

#endif
