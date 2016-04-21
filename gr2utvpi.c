#include "gr2utvpi.h"

int main(int argc, char *argv[]){
  if( argc != 3 ){
    fprintf(stderr, "Proper use: %s [input gr file] [output utvpi file]\n", argv[0]);
    exit(1);
  }
  
  Parser parser;
  parser.grFile = fopen(argv[1], "r");
  if( parser.grFile == NULL ){
    fprintf(stderr, "\"%s\" not found\n", argv[1]);
    exit(1);
  }
  parser.utvpiFile = fopen(argv[2], "w");
  if( parser.utvpiFile == NULL ){
    fprintf(stderr, "\"%s\" could not be created\n", argv[2]);
    exit(1);
  }
  parser.fileLineOffset = ftell( parser.grFile );
  parser.line = 1;
  parser.column = 0;
  parser.completedSuccessfully = true;
  
  parse(&parser);
  fclose(parser.grFile);
  if( !parser.completedSuccessfully ){
    fputs("Conversion failed.\n", stderr);
    fclose( parser.utvpiFile );
    remove( argv[2] );
  }
}

void parse(Parser * parser){
  Token token;
  do {  
    do {
      token = getToken(parser);
      if( token.type == NEWLINE ){
        fputc('\n', parser->utvpiFile);
	newLineAdjustment(parser);
      }
    } while( token.type == NEWLINE );
    switch(token.type){
    case C_KEYWORD:
      copyComment(parser);
      break;
    case P_KEYWORD:
      sp(parser);
      break;
    case A_KEYWORD:
      tailNodeID(parser);
      break;
    case END_OF_FILE:
      break;
    default:
      copyLine(parser);
    }
    newLineAdjustment(parser);
  } while( token.type != END_OF_FILE );
}

void copyComment(Parser * parser){
  fputs("//", parser->utvpiFile);
  fprintLineFromFilePointer(parser, parser->utvpiFile);
}

void copyLine(Parser * parser){
  fputs("//? ", parser->utvpiFile);
  fseek(parser->grFile, parser->fileLineOffset, SEEK_SET);
  fprintLineFromFilePointer(parser, parser->utvpiFile);
}

void sp(Parser * parser){
  Token token = getToken(parser);
  if( token.type == SP_KEYWORD ){
    token = getToken(parser);
    if( token.type == INTEGER ){
      int nodeCount = token.integerComponent;
      token = getToken(parser);
      if( token.type == INTEGER ){
        int arcCount = token.integerComponent;
	token = getToken(parser);
	if( token.type == NEWLINE || token.type == END_OF_FILE ){
	  fprintf(parser->utvpiFile, "%i variables // %i constraints\n", nodeCount, arcCount);
	  return;
	}
      }
    }
  }
  parseError(parser, "Unexpected token.\n");
}

void tailNodeID(Parser * parser){
  Token token = getToken(parser);
  if( token.type == INTEGER ){
    int tailNode = token.integerComponent;
    token = getToken(parser);
    if( token.type == INTEGER ){
      int headNode = token.integerComponent;
      token = getToken(parser);
      if( token.type == INTEGER ){
        int weight = token.integerComponent;
	token = getToken(parser);
	if( token.type == NEWLINE || token.type == END_OF_FILE ){
	  fprintf(parser->utvpiFile, "-x%i +x%i <= %i\n", tailNode, headNode, weight);
	  return;
	}
      }
    }
  }
  parseError(parser, "Unexpected token.\n");
}

void newLineAdjustment(Parser * parser){
  parser->fileLineOffset = ftell( parser->grFile );
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
  } while( isspace(character) && character != '\n' );
  switch(character){
  case '+':
    token.type = INTEGER;
    token.sign = PLUS;
    continuingInteger(parser, &token);
    break;
  case '-':
    token.type = INTEGER;
    token.sign = MINUS;
    continuingInteger(parser, &token);
    break;
  case 'c':
    if( endOfKeyword(parser) ){
      token.type = C_KEYWORD;
    }
    break;
  case 'p':
    if( endOfKeyword(parser) ){
      token.type = P_KEYWORD;
    }
    break;
  case 's':
    spKeyword(parser, &token);
    break;
  case 'a':
    if( endOfKeyword(parser) ){
      token.type = A_KEYWORD;
    }
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

void spKeyword(Parser * parser, Token * token){
  int character = myfgetc(parser);
  if( character == 'p' && endOfKeyword(parser) ){
    token->type = SP_KEYWORD;
  }
}

bool endOfKeyword(Parser * parser){
  int character = myfgetc(parser);
  myungetc(character, parser);
  return (isspace(character) || character == EOF);
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
  return fgetc(parser->grFile);
}

void myungetc(int character, Parser * parser){
  if(character != EOF){
    parser->column--;
    fseek(parser->grFile, -1, SEEK_CUR);
  }
}

void parseError(Parser * parser, const char * message){  
  parser->completedSuccessfully = false;
  fputs("Constraint System Error:\n", stderr);
  fprintf(stderr, "  Line %i Column %i\n", parser->line, parser->column);
  fputs("  ", stderr);
  fseek(parser->grFile, parser->fileLineOffset, SEEK_SET);
  fprintLineFromFilePointer(parser, stderr);
  for(int i = 0; i < parser->column + 1; i++){
    fputc(' ', stderr);
  }
  fputs("^\n", stderr);
  fprintf(stderr, "  %s\n", message);
}

void fprintLineFromFilePointer(Parser * parser, FILE * output){
  int character;
  do {
    character = fgetc( parser->grFile );
    if( character != EOF ){
      fputc(character, output);
    }
    else{
      fputc('\n', output);
    }
  } while( character != '\n' && character != EOF );
}