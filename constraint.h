#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

typedef enum ConstraintSign {
  CONSTRAINT_PLUS,
  CONSTRAINT_MINUS,
  CONSTRAINT_NONE,
} ConstraintSign;

typedef struct Constraint {
  ConstraintSign sign[2];
  int index[2];
  int weight;
} Constraint;

#endif