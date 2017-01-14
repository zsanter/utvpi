/*
 * utvpiInterpreterTester.c
 *
 * Call with [executable] [input file]
 * Any improperly-formatted input will be printed to stderr with an explanation by utvpiInterpreter. Any properly-formatted input
 * will be printed to stdout by the functions in this source file, so that it may be compared against the input file.
 */
 
#include <stdio.h>
#include <stdbool.h>
#include "utvpiInterpreter.h"

void main(int argc, char * argv[]);
static void initializeSystem(void * object, int n, Parser * parser);
static void addEdge(void * object, Constraint * constraint, Parser * parser);

void main(int argc, char * argv[]){
  if( argc == 2 ){
    FILE * input = fopen(argv[1], "r");
    if( input == NULL ){
      printf("File \"%s\" not found.\n", argv[1]);
    }
    else{
      if( parseFile(input, NULL, initializeSystem, addEdge) ){
        puts("Successful parse.");
      }
      else {
        puts("Parse unsuccessful.");
      }
    }
  }
  else{
    printf("%s [input file]\n", argv[0]);
  }
}

static void initializeSystem(void * object, int n, Parser * parser){
  printf("%i variables\n", n);
}

static void addEdge(void * object, Constraint * constraint, Parser * parser){
  for(int i = 0; i < 2; i++){
    if( constraint->sign[i] != CONSTRAINT_NONE){
      if( constraint->sign[i] == CONSTRAINT_PLUS){
        printf("+");
      }
      else if( constraint->sign[i] == CONSTRAINT_MINUS){
        printf("-");
      }
      printf("x%i ", constraint->index[i]);
    }
  }
  printf("<= %i\n", constraint->weight);
}
