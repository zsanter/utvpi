#ifndef _UTVPI_INTERPRETER_H
#define _UTVPI_INTERPRETER_H

#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>
#include "constraint.h"

typedef enum TokenType{
  SIGN,
  SIGNED_VARIABLE,
  UNSIGNED_VARIABLE,
  LESS_THAN_OR_EQUAL,
  INTEGER,
  NEWLINE,
  VARIABLES_KEYWORD,
  END_OF_FILE,
  UNDEFINED,
} TokenType;

typedef enum BlockCommentState {
  NOTHING,
  STAR,
  SLASH,
} BlockCommentState;

typedef struct Parser Parser;
typedef struct Token Token;

struct Parser {
  FILE * constraintFile;
  void * object;
  void (* initializeSystem)(void * object, int n, Parser * parser);
  void (* addEdge)(void * object, Constraint * constraint, Parser * parser);
  int systemMaxIndex;
  long int fileLineOffset;
  int line;
  int column;
  bool inMiddleOfBlockComment;
  bool completedSuccessfully;
};

struct Token {
  TokenType type;
  ConstraintSign sign;
  int integerComponent;
};

bool parseFile(FILE * constraintFile,
    void * object,
    void (* initializeSystem)(void * object, int n, Parser * parser), 
    void (* addEdge)(void * object, Constraint * constraint, Parser * parser));
bool parseVertexCount(Parser * parser);
void parseConstraints(Parser * parser);
void afterInitialSign(Parser * parser, Constraint * constraint);
void afterFirstVariable(Parser * parser, Constraint * constraint);
void secondVariable(Parser * parser, Constraint * constraint);
void lessThanOrEqual(Parser * parser, Constraint * constraint);
void weight(Parser * parser, Constraint * constraint);
void newLine(Parser * parser, Constraint * constraint);
bool variableIndexWithinBounds(Parser * parser, Token * token);
bool variableIndecesDiffer(Parser * parser, Constraint * constraint);
void newLineAdjustment(Parser * parser);
Token getToken(Parser * parser);
void afterSign(Parser * parser, Token * token);
void afterX(Parser * parser, Token * token);
void variablesKeyword(Parser * parser, Token * token);
void lessThanOrEqualToken(Parser * parser, Token * token);
void continuingInteger(Parser * parser, Token * token);
int myfgetc(Parser * parser);
void myungetc(int character, Parser * parser);
void parseError(Parser * parser, const char * message);

#endif
