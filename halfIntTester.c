/*
 * halfIntTester.c
 *
 * Call with [executable]
 * Was used to ensure that half_ints exhibited the desired behavior.
 */
 
#include "halfint.h"
#include <stdio.h>
#include <stdbool.h>

void main(){

  int inputs[10] = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4};
  for(int i = 0; i < 10; i++){
    printf("for input integer %i\n", inputs[i]);
    printf("\thalf of %i is %i\n", inputs[i], inputs[i]/2);
    half_int divByTwo = intDivBy2ToHalfInt(inputs[i]);
    printf("\tdivByTwo value is %i as a straight integer, and represents %f as a float and %i when converted to an integer\n", divByTwo, halfIntToDouble(divByTwo), halfIntToInt(divByTwo));
    printf("\tdivByTwo value is seen as ");
    if( !halfIntIsIntegral(divByTwo) ){
      printf("not ");
    }
    printf("integral\n");
    half_int ceiling = halfIntCeil( divByTwo );
    printf("\tceil(divByTwo) value is %i as a straight integer, and represents %f as a float and %i when converted to an integer\n", ceiling, halfIntToDouble(ceiling), halfIntToInt(ceiling));
    half_int floor = halfIntFloor( divByTwo );
    printf("\tfloor(divByTwo) value is %i as a straight integer, and represents %f as a float and %i when converted to an integer\n", floor, halfIntToDouble(floor), halfIntToInt(floor));
    half_int cast = intToHalfInt(inputs[i]);
    printf("\tcast value is %i as a straight integer, and represents %f as a float and %i when converted to an integer\n", cast, halfIntToDouble(cast), halfIntToInt(cast));
  }

}
