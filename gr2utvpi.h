#ifndef _GR_2_UTVPI_H
#define _GR_2_UTVPI_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>

typedef enum TokenType{
  INTEGER,
  NEWLINE,
  C_KEYWORD,
  P_KEYWORD,
  SP_KEYWORD,
  A_KEYWORD,
  END_OF_FILE,
  UNDEFINED,
} TokenType;

typedef enum Sign {
  PLUS,
  MINUS,
  NONE,
} Sign;

typedef struct Parser Parser;
typedef struct Token Token;
typedef struct Constraint Constraint;

struct Parser {
  FILE * grFile;
  FILE * utvpiFile;
  long int fileLineOffset;
  int line;
  int column;
  bool completedSuccessfully;
};

struct Token {
  TokenType type;
  Sign sign;
  int integerComponent;
};

int main(int argc, char *argv[]);
void parse(Parser * parser);
void copyComment(Parser * parser);
void copyLine(Parser * parser);
void sp(Parser * parser);
void tailNodeID(Parser * parser);
void newLineAdjustment(Parser * parser);
Token getToken(Parser * parser);
void spKeyword(Parser * parser, Token * token);
bool endOfKeyword(Parser * parser);
void continuingInteger(Parser * parser, Token * token);
int myfgetc(Parser * parser);
void myungetc(int character, Parser * parser);
void parseError(Parser * parser, const char * message);
void fprintLineFromFilePointer(Parser * parser, FILE * output);

#endif