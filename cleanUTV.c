#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

typedef enum State {
  NUM_VARS_POUND,
  NUM_VARS_V,
  NUM_VARS_A,
  NUM_VARS_R,
  NUM_VARS_S,
  NUM_VARS_SPACE1,
  NUM_VARS_COLON,
  NUM_VARS_SPACE2,
  NUM_VARS_NUMBER,
  CONS_COLON,
  CONS_C,
  CONS_O,
  CONS_N,
  CONS_S,
  NOTHING,
} State;

void main(int argc, char * argv[]){
  if(argc != 2){
    fprintf(stderr, "Proper use: %s [uncleaned utv file]\n", argv[0]);
    exit(1);
  }
  FILE * input = fopen(argv[1], "r");
  if(input == NULL){
    fprintf(stderr, "\"%s\" not found.\n", argv[1]);
    exit(2);
  }
  char * tempName = "cleanUTVtempFile.txt";
  FILE * temp = fopen(tempName, "w");
  if(temp == NULL){
    fputs("Error creating temporary file.\n", stderr);
    exit(3);
  }
  int numVariables = 0;
  int currChar = fgetc(input);
  State currState = NOTHING;
  fputs("/*\n", temp);
  while( currChar != EOF ){
    fputc(currChar, temp);
    switch(currState){
    case NOTHING:
      switch(currChar){
      case '#':
        currState = NUM_VARS_POUND;
        break;
      case ':':
        currState = CONS_COLON;
      }
      break;
    case NUM_VARS_POUND:
      switch(currChar){
      case 'v':
        currState = NUM_VARS_V;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case NUM_VARS_V:
      switch(currChar){
      case 'a':
        currState = NUM_VARS_A;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case NUM_VARS_A:
      switch(currChar){
      case 'r':
        currState = NUM_VARS_R;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case NUM_VARS_R:
      switch(currChar){
      case 's':
        currState = NUM_VARS_S;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case NUM_VARS_S:
      switch(currChar){
      case ' ':
        currState = NUM_VARS_SPACE1;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case NUM_VARS_SPACE1:
      switch(currChar){
      case ':':
        currState = NUM_VARS_COLON;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case NUM_VARS_COLON:
      switch(currChar){
      case ' ':
        currState = NUM_VARS_SPACE2;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case NUM_VARS_SPACE2:
      if( isdigit( currChar ) ){
        numVariables = currChar - '0';
        currState = NUM_VARS_NUMBER;
      }
      else{
        currState = NOTHING;
      }
      break;
    case NUM_VARS_NUMBER:
      if( isdigit( currChar ) ){
        numVariables = numVariables * 10 + currChar - '0';
      }
      else{
        currState = NOTHING;
      }
      break;
    case CONS_COLON:
      switch(currChar){
      case 'c':
        currState = CONS_C;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case CONS_C:
      switch(currChar){
      case 'o':
        currState = CONS_O;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case CONS_O:
      switch(currChar){
      case 'n':
        currState = CONS_N;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case CONS_N:
      switch(currChar){
      case 's':
        currState = CONS_S;
        break;
      default:
        currState = NOTHING;
      }
      break;
    case CONS_S:
      switch(currChar){
      case '\r':
        break;
      case '\n':
        fprintf(temp, "*/\n%d variables\n", numVariables);
      default:
        currState = NOTHING;
      }
    }
    currChar = fgetc(input);
  }
  fclose(input);
  fclose(temp);
  remove(argv[1]);
  rename(tempName, argv[1]);
}