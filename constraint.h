/*
 * constraint.h
 *
 * A standard representation of UTVPI constraints, portable across inequality system solvers, input file parsers, and anything 
 * else. It is expected that inequality system solvers will translate Constraints into their own internal representations and can
 * translate their own internal representations back into Constraints for returning proofs of infeasibility, or for any other
 * purpose.
 */

#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

#include <stdio.h>

#define VARIABLES_PER_CONSTRAINT 2

/*
 * CONSTRAINT_PLUS and CONSTRAINT_MINUS represent positive and negative variable signs. CONSTRAINT_NONE represents that a variable
 * is not present, i.e. that the other variable within the Constraint struct is the only variable present within an absolute 
 * constraint.
 */
typedef enum ConstraintSign {
  CONSTRAINT_PLUS,
  CONSTRAINT_MINUS,
  CONSTRAINT_NONE,
} ConstraintSign;

typedef struct Constraint Constraint;
typedef struct ConstraintRefList ConstraintRefList;
typedef struct ConstraintRefListNode ConstraintRefListNode;

/*
 * sign[i] and index[i] represent all information about variable i, and weight gives the weight of the constraint. By convention,
 * absolute constraints, consisting of only a single variable and a weight, have sign[1] == CONSTRAINT_NONE, and index[1] == 0.
 * The single variable should be represented by sign[0] and index[0].
 */
struct Constraint {
  ConstraintSign sign[VARIABLES_PER_CONSTRAINT];
  int index[VARIABLES_PER_CONSTRAINT];
  int weight;
};

/*
 * The ConstraintRefList is a way to reference a list of Constraints
 *
 * first - pointer to the first ConstraintRefListNode within the ConstraintRefList
 * last - pointer to the last ConstraintRefListNode within the ConstraintRefList
 */
struct ConstraintRefList {
  ConstraintRefListNode * first;
  ConstraintRefListNode * last;
};

/*
 * The ConstraintRefListNode is one node within the ConstraintRefList
 *
 * constraint - pointer to a Constraint
 * next - pointer to the next ConstraintRefListNode within the ConstraintRefList
 */
struct ConstraintRefListNode {
  Constraint * constraint;
  ConstraintRefListNode * next;
};

void fputConstraint(Constraint * constraint, FILE * output);
ConstraintRefList * generateConstraintRefList();
void constraintRefListAppend(ConstraintRefList * crl, Constraint * constraint);
void constraintRefListPrepend(ConstraintRefList * crl, Constraint * constraint);
void freeConstraintRefList(ConstraintRefList * crl);

#endif