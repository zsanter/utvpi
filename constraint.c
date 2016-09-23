#include "constraint.h"

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
  ConstraintRefList * newCRL = (ConstraintRefList *) malloc( sizeof(ConstraintRefList) );
  newCRL->first = NULL;
  newCRL->last = NULL;
  return newCRL;
}

/*
 * constraintRefListAppend() adds a Constraint pointer to the end of a ConstraintRefList.
 * crl - pointer to a ConstraintRefList to append a Constraint pointer to
 * constraint - pointer to a constraint, which will be added to the end of crl
 */
void constraintRefListAppend(constraintRefList * crl, Constraint * constraint){
  ConstraintRefListNode * newCRLN = (ConstraintRefListNode *) malloc( sizeof(ConstraintRefListNode) );
  newCRLN->constraint = constraint;
  newCRLN->next = NULL;
  if( crl->first == NULL ){
    crl->first = newCRLN;
  }
  else {
    crl->last->next = newCRLN;
  }
  crl->last = newCRLN;
}

/*
 * constraintRefListPrepend() adds a constraint pointer to the beginning of a ConstraintRefList.
 * crl - pointer to a ConstraintRefList to prepend with a Constraint pointer
 * constraint - pointer to a Constraint, which will be added to the beginning of crl
 */
void constraintRefListPrepend(ConstraintRefList * crl, Constraint * constraint){
  ConstraintRefListNode * newCRLN = (ConstraintRefListNode *) malloc( sizeof(ConstraintRefListNode) );
  newCRLN->constraint = constraint;
  newCRLN->next = crl->first;
  crl->first = newCRLN;
  if( crl->last == NULL ){
    crl->last = newCRLN;
  }
}

/*
 * freeConstraintRefList() frees the input constraintRefList, as well as referenced Constraints.
 * crl - pointer to the ConstraintRefList to be freed
 */
void freeConstraintRefList(ConstraintRefList * crl){
  if( crl != NULL ){
    ConstraintRefListNode * crln = crl->first;
    while( crln != NULL ){
      ConstraintRefListNode * oldCRLN = crln;
      crln = crln->next;
      free( oldCRLN->constraint );
      free( oldCRLN );
    }
    free( crl );
  }
}