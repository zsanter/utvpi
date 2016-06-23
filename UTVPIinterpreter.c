#include "UTVPIinterpreter.h"

bool parseFile(FILE * constraintFile,
        void * object,
        void (* initializeSystem)(void * object, int n, Parser * parser), 
        void (* addEdge)(void * object, Constraint * constraint, Parser * parser)){
  
  Parser parser;
  parser.constraintFile = constraintFile;
  parser.object = object;
  parser.initializeSystem = initializeSystem;
  parser.addEdge = addEdge;
  parser.systemMaxIndex = 0;
  parser.fileLineOffset = ftell( constraintFile );
  parser.line = 1;
  parser.column = 0;
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

bool parseVertexCount(Parser * parser){
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

void parseConstraints(Parser * parser){
  Token token;
  do {

    Constraint constraint;
    constraint.sign[0] = NONE;
    constraint.index[0] = 0;
    constraint.sign[1] = NONE;
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
    case SIGNED_VARIABLE:
      constraint.sign[0] = token.sign;
      constraint.index[0] = token.integerComponent;
      if( variableIndexWithinBounds( parser, &token ) ){
        afterFirstVariable(parser, &constraint);
      }
      break;
    case UNSIGNED_VARIABLE:
      constraint.sign[0] = PLUS;
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

void afterInitialSign(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch( token.type ){
  case UNSIGNED_VARIABLE:
    constraint->index[0] = token.integerComponent;
    if( variableIndexWithinBounds( parser, &token ) ){
      afterFirstVariable(parser, constraint);
    }
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

void afterFirstVariable(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch( token.type ){
  case SIGN:
    constraint->sign[1] = token.sign;
    secondVariable(parser, constraint);
    break;
  case SIGNED_VARIABLE:
    constraint->sign[1] = token.sign;
    constraint->index[1] = token.integerComponent;
    if( variableIndexWithinBounds( parser, &token )
        && variableIndecesDiffer( parser, constraint ) ){
      lessThanOrEqual(parser, constraint);
    }
    break;
  case LESS_THAN_OR_EQUAL:
    weight(parser, constraint);
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

void secondVariable(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch(token.type){
  case SIGNED_VARIABLE:
    if( constraint->sign[1] != token.sign ){
      constraint->sign[1] = MINUS;
    }
    else{
      constraint->sign[1] = PLUS;
    }
    constraint->index[1] = token.integerComponent;
    if( variableIndexWithinBounds( parser, &token )
        && variableIndecesDiffer( parser, constraint ) ){
      lessThanOrEqual(parser, constraint);
    }
    break;
  case UNSIGNED_VARIABLE:
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

void lessThanOrEqual(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch(token.type){
  case LESS_THAN_OR_EQUAL:
    weight(parser, constraint);
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

void weight(Parser * parser, Constraint * constraint){
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

void newLine(Parser * parser, Constraint * constraint){
  Token token = getToken(parser);
  switch(token.type){
  case NEWLINE:
  case END_OF_FILE:
    parser->addEdge(parser->object, constraint, parser);
    break;
  default:
    parseError(parser, "Unexpected token.");
  }
}

bool variableIndexWithinBounds(Parser * parser, Token * token){
  if( token->integerComponent < 1 || token->integerComponent > parser->systemMaxIndex ){
    parseError(parser, "Variable index outside of defined range.");
    return false;
  }
  return true;
}

bool variableIndecesDiffer(Parser * parser, Constraint * constraint){
  if( constraint->index[0] == constraint->index[1] ){
    parseError(parser, "Variable indices must differ.");
    return false;
  }
  return true;
}

void newLineAdjustment(Parser * parser){
  parser->fileLineOffset = ftell( parser->constraintFile );
  parser->line++;
  parser->column = 0;
}

Token getToken(Parser * parser){

  Token token;
  token.type = UNDEFINED;
  token.sign = NONE;
  token.integerComponent = 0;

  int character;
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
        ;
        BlockCommentEnd bce = NOTHING;
        do {
          character = myfgetc(parser);
          if(character == '\n'){
            newLineAdjustment(parser);
          }
          if(character == '*'){
            bce = STAR;
          }
          else if(bce == STAR && character == '/'){
            bce = SLASH;
          }
          else{
            bce = NOTHING;
          }
        } while( bce != SLASH && character != EOF );
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
    token.sign = PLUS;
    afterSign(parser, &token);
    break;
  case '-':
    token.sign = MINUS;
    afterSign(parser, &token);
    break;
  case 'x':
  case 'X':
    token.type = UNSIGNED_VARIABLE;
    token.sign = PLUS;
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
      token.sign = PLUS;
      token.integerComponent = character - '0';
      continuingInteger(parser, &token);	
    }
  }
  if( token.type == INTEGER && token.sign == MINUS ){
    token.integerComponent *= -1;
    token.sign = PLUS;
  } 
  return token;
}

void afterSign(Parser * parser, Token * token){
  int character = myfgetc(parser);
  if( character == 'x' || character == 'X'){
    token->type = SIGNED_VARIABLE;
    afterX(parser, token);
  }
  else if( isdigit(character) ){
    token->type = INTEGER;
    token->integerComponent = character - '0';
    continuingInteger(parser, token);
  }
  else if( isspace(character) || character == EOF ){
    token->type = SIGN;
    myungetc( character, parser );
  }
  else {
    token->type = UNDEFINED;
    myungetc( character, parser );
  }
}

void afterX(Parser * parser, Token * token){
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

void variablesKeyword(Parser * parser, Token * token){
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

void lessThanOrEqualToken(Parser * parser, Token * token){
  int character = myfgetc(parser);
  if( character != '=' ){
    token->type = UNDEFINED;
    myungetc( character, parser );
  }
}

void continuingInteger(Parser * parser, Token * token){
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

int myfgetc(Parser * parser){
  parser->column++;
  return fgetc(parser->constraintFile);
}

void myungetc(int character, Parser * parser){
  if(character != EOF){
    parser->column--;
    fseek(parser->constraintFile, -1, SEEK_CUR);
  }
}

void parseError(Parser * parser, const char * message){  
  parser->completedSuccessfully = false;
  fputs("Constraint System Error:\n", stderr);
  fprintf(stderr, "  Line %i Column %i\n", parser->line, parser->column);
  fputs("  ", stderr);
  fseek(parser->constraintFile, parser->fileLineOffset, SEEK_SET);
  int character;
  do {
    character = fgetc(parser->constraintFile);
    if( character != EOF ){
      fputc(character, stderr);
    }
    else{
      fputc('\n', stderr);
    }
  } while( character != '\n' && character != EOF );
  for(int i = 0; i < parser->column + 1; i++){
    fputc(' ', stderr);
  }
  fputs("^\n", stderr);
  fprintf(stderr, "  %s\n", message);
}
