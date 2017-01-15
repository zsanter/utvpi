/**/

#include <stdio.h>
#include <stdbool.h>
#include "utvpiInterpreter.h"

typedef struct Delayer {
  ConstraintRefList * constraints;
  int n;
} Delayer;

void delayedInitializeSystem(void * object, int n, Parser * parser);
void delayedAddConstraint(void * object, Constraint * constraint, Parser * parser);

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
bool parseFileDelayedSysGen(FILE * constraintFile,
    void * object,
    void (* initializeSystem)(void * object, int n, Parser * parser),
    void (* addConstraint)(void * object, Constraint * constraint, Parser * parser)){

  Delayer delayer;
  delayer.constraints = generateConstraintRefList();
  delayer.n = 0;

  bool parseSuccesful = parseFile(constraintFile, delayer, delayedInitializeSystem, delayedAddConstraint);
  if( !parseSuccessful ){
    freeConstraintRefList( delayer.constraints );
    return false;
  }
  else{
    initializeSystem(object, delayer.n, NULL);
    Constraint * constraint = constraintRefListNext( delayer.constraints );
    while( constraint != NULL ){
      addConstraint(object, constraint, NULL);
      constraint = constraintRefListNext( delayer.constraints );
    }
    freeConstraintRefList( constraints );
    return true;
  }
}

void delayedInitializeSystem(void * object, int n, Parser * parser){
  Delayer * delayer = (Delayer *) object;
  delayer->n = n;
}

void delayedAddConstraint(void * object, Constraint * constraint, Parser * parser){
  Delayer * delayer = (Delayer *) object;

  Constraint * heapConstraint = (Constraint *) malloc( sizeof( Constraint ) );
  heapConstraint->sign[0] = constraint->sign[0];
  heapConstraint->sign[1] = constraint->sign[1];
  heapConstraint->index[0] = constraint->index[0];
  heapConstraint->index[1] = constraint->index[1];
  heapConstraint->weight = constraint->weight;

  constraintRefListPrepend(delayer.constraints, heapConstraint);
}
