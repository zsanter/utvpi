#include <ctype.h>
#include "utvpiInterpreter.h"

typedef enum TokenType {
  SIGN,
  VARIABLE,
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

struct Parser {
  FILE * constraintFile;
  void * object;
  void (* initializeSystem)(void * object, int n, Parser * parser);
  void (* addConstraint)(void * object, Constraint * constraint, Parser * parser);
  int systemMaxIndex;
  long int fileLineOffset;
  int line;
  int column;
  int characterStore;
  bool characterStored;
  bool inMiddleOfBlockComment;
  bool completedSuccessfully;
};

typedef struct Token {
  TokenType type;
  ConstraintSign sign;
  int integerComponent;
} Token;

static bool parseVertexCount(Parser * parser);
static void parseConstraints(Parser * parser);
static void afterInitialSign(Parser * parser, Constraint * constraint);
static void afterFirstVariable(Parser * parser, Constraint * constraint);
static void secondVariable(Parser * parser, Constraint * constraint);
static void lessThanOrEqual(Parser * parser, Constraint * constraint);
static void weight(Parser * parser, Constraint * constraint);
static void newLine(Parser * parser, Constraint * constraint);
static bool variableIndexWithinBounds(Parser * parser, Token * token);
static bool variableIndecesDiffer(Parser * parser, Constraint * constraint);
static void newLineAdjustment(Parser * parser);
static Token getToken(Parser * parser);
static void afterSign(Parser * parser, Token * token);
static void afterX(Parser * parser, Token * token);
static void variablesKeyword(Parser * parser, Token * token);
static void lessThanOrEqualToken(Parser * parser, Token * token);
static void continuingInteger(Parser * parser, Token * token);
static int myfgetc(Parser * parser);
static void myungetc(int character, Parser * parser);

bool parseFile(FILE * constraintFile,
        void * object,
        void (* initializeSystem)(void * object, int n, Parser * parser),
        void (* addConstraint)(void * object, Constraint * constraint, Parser * parser)){

  Parser parser;
  parser.constraintFile = constraintFile;
  parser.object = object;
  parser.initializeSystem = initializeSystem;
  parser.addConstraint = addConstraint;
  parser.systemMaxIndex = 0;
  parser.fileLineOffset = ftell( constraintFile );
  parser.line = 1;
  parser.column = 0;
  parser.characterStore = 0;
  parser.characterStored = false;
  parser.inMiddleOfBlockComment = false;
  parser.completedSuccessfully = true;

  bool foundVertexCount = parseVertexCount(&parser);
  if( !foundVertexCount ){
    parseError(&parser, "No graph definition line found.");
  }
  else if( parser.systemMaxIndex <= 0 ){
    parseError(&parser, "Variable count must be positive.");
  }
  else{
    parser.initializeSystem(parser.object, parser.systemMaxIndex, &parser);
    newLineAdjustment(&parser);
    parseConstraints(&parser);
  }
  return parser.completedSuccessfully;
}

static bool parseVertexCount(Parser * parser){
  Token token;
  do {
    token = getToken(parser);
    if( token.type == NEWLINE ){
      newLineAdjustment(parser);
    }
  } while( token.type == NEWLINE );
  if( token.type == INTEGER ){
    parser->systemMaxIndex = token.integerComponent;
    token = getToken(parser);
    if( token.type == VARIABLES_KEYWORD ){
      token = getToken(parser);
      if( token.type == NEWLINE || token.type == END_OF_FILE ){
        return true;
      }
    }
  }
  return false;
}

static void parseConstraints(Parser * parser){
  Token token;
  do {

    Constraint constraint;
    constraint.sign[0] = CONSTRAINT_NONE;
    constraint.index[0] = 0;
    constraint.sign[1] = CONSTRAINT_NONE;
    constraint.index[1] = 0;
    constraint.weight = 0;

    do {
      token = getToken(parser);
      if( token.type == NEWLINE ){
        newLineAdjustment(parser);
      }
    } while( token.type == NEWLINE );
    switch(token.type){
    case SIGN:
      constraint.sign[0] = token.sign;
      afterInitialSign(parser, &constraint);
      break;
    case VARIABLE:
      constraint.sign[0] = CONSTRAINT_PLUS;
      constraint.index[0] = token.integerComponent;
      if( variableIndexWithinBounds( parser, &token ) ){
        afterFirstVariable(parser, &constraint);
      }
      break;
    case END_OF_FILE:
      break;
    default:
      parseError(parser, "Unexpected token.");
    }
    newLineAdjustment(parser);
  } while( token.type != END_OF_FILE );
}

static void afterInitialSign(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch( token.type ){
  case VARIABLE:
    constraint->index[0] = token.integerComponent;
    if( variableIndexWithinBounds( parser, &token ) ){
      afterFirstVariable(parser, constraint);
    }
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

static void afterFirstVariable(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch( token.type ){
  case SIGN:
    constraint->sign[1] = token.sign;
    secondVariable(parser, constraint);
    break;
  case LESS_THAN_OR_EQUAL:
    weight(parser, constraint);
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

static void secondVariable(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch(token.type){
  case VARIABLE:
    constraint->index[1] = token.integerComponent;
    if( variableIndexWithinBounds( parser, &token )
        && variableIndecesDiffer( parser, constraint ) ){
      lessThanOrEqual(parser, constraint);
    }
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

static void lessThanOrEqual(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch(token.type){
  case LESS_THAN_OR_EQUAL:
    weight(parser, constraint);
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

static void weight(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch(token.type){
  case INTEGER:
    constraint->weight = token.integerComponent;
    newLine(parser, constraint);
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

static void newLine(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch(token.type){
  case NEWLINE:
  case END_OF_FILE:
    parser->addConstraint(parser->object, constraint, parser);
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

static bool variableIndexWithinBounds(Parser * parser, Token * token){
  if( token->integerComponent < 1 || token->integerComponent > parser->systemMaxIndex ){
    parseError(parser, "Variable index outside of defined range.");
    return false;
  }
  return true;
}

static bool variableIndecesDiffer(Parser * parser, Constraint * constraint){
  if( constraint->index[0] == constraint->index[1] ){
    parseError(parser, "Variable indices must differ.");
    return false;
  }
  return true;
}

static void newLineAdjustment(Parser * parser){
  parser->fileLineOffset = ftell( parser->constraintFile );
  parser->line++;
  parser->column = 0;
}

static Token getToken(Parser * parser){

  Token token;
  token.type = UNDEFINED;
  token.sign = CONSTRAINT_NONE;
  token.integerComponent = 0;

  int character;
  if( parser->inMiddleOfBlockComment ){
    parser->inMiddleOfBlockComment = false;
    goto inMiddleOfBlockCommentLabel;
  }
  do {
    character = myfgetc(parser);
    if(character == '/'){
      character = myfgetc(parser);
      switch(character){
      case '/':
        do {
          character = myfgetc(parser);
        } while(character != '\n' && character != EOF);
        break;
      case '*':
        inMiddleOfBlockCommentLabel:
        ;
        BlockCommentState bcs = NOTHING;
        do {
          character = myfgetc(parser);
          if(character == '\n'){
            newLineAdjustment(parser);
          }
          if(character == '*'){
            bcs = STAR;
          }
          else if(bcs == STAR && character == '/'){
            bcs = SLASH;
          }
          else{
            bcs = NOTHING;
          }
        } while( bcs != SLASH && character != EOF );
        if( character == EOF ){
          parseError(parser, "Unexpected EOF or bad block comment.");
        }
        character = myfgetc(parser);
        break;
      default:
        myungetc(character, parser);
        character = '/';
      }
    }
  } while( isspace(character) && character != '\n' );
  switch(character){
  case '+':
    token.sign = CONSTRAINT_PLUS;
    afterSign(parser, &token);
    break;
  case '-':
    token.sign = CONSTRAINT_MINUS;
    afterSign(parser, &token);
    break;
  case 'x':
  case 'X':
    token.type = VARIABLE;
    token.sign = CONSTRAINT_PLUS;
    afterX(parser, &token);
    break;
  case 'v':
    token.type = VARIABLES_KEYWORD;
    variablesKeyword(parser, &token);
    break;
  case '<':
    token.type = LESS_THAN_OR_EQUAL;
    lessThanOrEqualToken(parser, &token);
    break;
  case '\n':
    token.type = NEWLINE;
    break;
  case EOF:
    token.type = END_OF_FILE;
    break;
  default:
    if( isdigit(character) ){
      token.type = INTEGER;
      token.sign = CONSTRAINT_PLUS;
      token.integerComponent = character - '0';
      continuingInteger(parser, &token);
    }
  }
  if( token.type == INTEGER && token.sign == CONSTRAINT_MINUS ){
    token.integerComponent *= -1;
    token.sign = CONSTRAINT_PLUS;
  }
  return token;
}

static void afterSign(Parser * parser, Token * token){
  int character = myfgetc(parser);
  if( isdigit(character) ){
    token->type = INTEGER;
    token->integerComponent = character - '0';
    continuingInteger(parser, token);
  }
  else {
    token->type = SIGN;
    myungetc( character, parser );
  }
}

static void afterX(Parser * parser, Token * token){
  int character = myfgetc(parser);
  if( character == '_' ){
    character = myfgetc(parser);
  }
  if( isdigit(character) ){
    token->integerComponent = character - '0';
    continuingInteger(parser, token);
  }
  else{
    token->type = UNDEFINED;
    myungetc( character, parser );
  }
}

static void variablesKeyword(Parser * parser, Token * token){
  int character;
  char ariables[] = {'a', 'r', 'i', 'a', 'b', 'l', 'e', 's'};
  bool matches = true;
  for(int i = 0; i < 8 && matches; i++){
    character = myfgetc(parser);
    matches = (character == ariables[i]);
  }
  if( !matches ){
    token->type = UNDEFINED;
    myungetc( character, parser );
  }
}

static void lessThanOrEqualToken(Parser * parser, Token * token){
  int character = myfgetc(parser);
  if( character != '=' ){
    token->type = UNDEFINED;
    myungetc( character, parser );
  }
}

static void continuingInteger(Parser * parser, Token * token){
  int character = myfgetc(parser);
  while( isdigit( character ) ){
    token->integerComponent *= 10;
    token->integerComponent += ( character - '0' );
    character = myfgetc(parser);
  }
  if( isalpha( character ) ){
    token->type = UNDEFINED;
  }
  myungetc( character, parser );
}

static int myfgetc(Parser * parser){
  parser->column++;
  if( parser->characterStored ){
    parser->characterStored = false;
    return parser->characterStore;
  }
  else {
    return fgetc(parser->constraintFile);
  }
}

static void myungetc(int character, Parser * parser){
  if( parser->characterStored ){
    fputs("failure\n", stderr);
  }
  parser->column--;
  parser->characterStore = character;
  parser->characterStored = true;
}

void parseError(Parser * parser, const char * message){
  parser->characterStored = false;
  parser->completedSuccessfully = false;
  fputs("Constraint System Error:\n", stderr);
  fprintf(stderr, "  Line %i Column %i\n", parser->line, parser->column);
  fputs("  ", stderr);
  fseek(parser->constraintFile, parser->fileLineOffset, SEEK_SET);
  int character;
  BlockCommentState bcs = NOTHING;
  do {
    character = fgetc(parser->constraintFile);
    if( character == EOF ){
      fputc('\n', stderr);
    }
    else {
      fputc(character, stderr);
      if( parser->inMiddleOfBlockComment ){
        if( character == '*' ){
          bcs = STAR;
        }
        else if( bcs == STAR && character == '/' ){
          parser->inMiddleOfBlockComment = false;
          bcs = NOTHING;
        }
        else {
          bcs = NOTHING;
        }
      }
      else {
        if( character == '/' ){
          bcs = SLASH;
        }
        else if( bcs == SLASH && character == '*' ){
          parser->inMiddleOfBlockComment = true;
          bcs = NOTHING;
        }
        else {
          bcs = NOTHING;
        }
      }
    }
  } while( character != '\n' && character != EOF );
  for(int i = 0; i < parser->column + 1; i++){
    fputc(' ', stderr);
  }
  fputs("^\n", stderr);
  fprintf(stderr, "  %s\n", message);
}
