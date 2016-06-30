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

/*
 * sign[i] and index[i] represent all information about variable i, and weight gives the weight of the constraint. By convention,
 * absolute constraints, consisting of only a single variable and a weight, have sign[1] == CONSTRAINT_NONE, and index[1] == 0.
 * The single variable should be represented by sign[0] and index[0].
 */
typedef struct Constraint {
  ConstraintSign sign[2];
  int index[2];
  int weight;
} Constraint;

#endif