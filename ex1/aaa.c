#include  <stdio.h>
#include  <math.h>
#include <stdlib.h>
#include <math.h>

int main()
{
  float x;
  
  while(1==1){
    printf("Gi inn et tall 0 avslutter):\n");
    scanf("%f", &x);
    
    if ( x == 0 ) exit(0);
    
    printf("Tallet er ............: %f\n", x);
    printf("Tallet kvarert er ..: : %f\n", x*x);
    printf("Roten av tallet er ...: %f\n", sqrt(x));
  }
}

