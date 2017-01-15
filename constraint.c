#include <stdlib.h>
#include <stdbool.h>
#include "constraint.h"

/*
 * fputConstraint() prints the constraint equation corresponding to constraint to output
 *
 * constraint - pointer to a Constraint struct to print as a constraint equation
 * output - FILE pointer to print the constraint equation to
 */
void fputConstraint(Constraint * constraint, FILE * output){
  for(int i = 0; i < VARIABLES_PER_CONSTRAINT; i++){
    switch( constraint->sign[i] ){
    case CONSTRAINT_PLUS:
      fprintf(output, "+x%d ", constraint->index[i]);
      break;
    case CONSTRAINT_MINUS:
      fprintf(output, "-x%d ", constraint->index[i]);
      break;
    }
  }
  fprintf(output, "<= %d\n", constraint->weight);
}

/*
 * generateConstraintRefList() allocates and initializes a ConstraintRefList, returning a pointer to it.
 */
ConstraintRefList * generateConstraintRefList(){
  return (ConstraintRefList *) generateVoidRefList();
}

/*
 * constraintRefListAppend() adds a Constraint pointer to the end of a ConstraintRefList.
 *
 * crl - pointer to a ConstraintRefList to append a Constraint pointer to
 * constraint - pointer to a constraint, which will be added to the end of crl
 */
void constraintRefListAppend(ConstraintRefList * crl, Constraint * constraint){
  voidRefListAppend((VoidRefList *) crl, (void *) constraint);
}

/*
 * constraintRefListPrepend() adds a constraint pointer to the beginning of a ConstraintRefList.
 *
 * crl - pointer to a ConstraintRefList to prepend with a Constraint pointer
 * constraint - pointer to a Constraint, which will be added to the beginning of crl
 */
void constraintRefListPrepend(ConstraintRefList * crl, Constraint * constraint){
  voidRefListPrepend((VoidRefList *) crl, (void *) constraint);
}

/*
 * constraintRefListNext() allows for the ConstraintRefList to be iterated through. Each call to this function returns another
 * element within the ConstraintRefList, until the end of the list is reached, when NULL is returned. New elements should not be
 * added to the ConstraintRefList while this process is ongoing. Reaching the end of the list resets the internal iterator which
 * enables this process. If iteration through the list ends before the list has been fully traversed,
 * constraintRefListIteratorReset() should be called.
 *
 * crl - the ConstraintRefLIst to be iterated through
 */
Constraint * constraintRefListNext(ConstraintRefList * crl){
  return (Constraint *) voidRefListNext((VoidRefList *) crl);
}

/*
 * constraintRefListIteratorReset() resets the ConstraintRefList's internal iterator after a partial traverse of the
 * ConstraintRefLIst has been conducted using constraintRefListNext().
 *
 * crl - the ConstraintRefList
 */
void constraintRefListIteratorReset(ConstraintRefList * crl){
  voidRefListIteratorReset((VoidRefList *) crl);
}


/*
 * freeConstraintRefList() frees the input constraintRefList, as well as referenced Constraints.
 *
 * crl - pointer to the ConstraintRefList to be freed
 */
void freeConstraintRefList(ConstraintRefList * crl){
  freeVoidRefList((VoidRefList *) crl, true);
}
