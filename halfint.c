#include "halfint.h"

half_int intToHalfInt(int input){
  return (half_int)(input << 1);
}

int halfIntToInt(half_int input){
  int output = (int)(input >> 1);
  if( input < 0 && ( input & 0x00000001 ) ){
    output += 1;
  }
  return output;
}

double halfIntToDouble(half_int input){
  double output = (double)(input >> 1);
  if( input & 0x00000001 ){
    output += 0.5;
  }
  return output;
}

half_int intDivBy2ToHalfInt(int input){
  return (half_int)input;
}

bool halfIntIsIntegral(half_int input){
  return !(bool)( input & 0x00000001 );  
}

half_int halfIntCeil(half_int input){
  if( input & 0x00000001 ){
    input += (half_int) 0x00000001;
  }
  return input;
}

half_int halfIntFloor(half_int input){
  if( input & 0x00000001 ){
    input -= (half_int) 0x00000001;
  }
  return input;
}
