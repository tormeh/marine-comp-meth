#include  <stdio.h>
#include  <math.h>
#include <stdlib.h>
#include <math.h>

int main()
{
  int a;
  int b;
  int n;
  
    scanf("%d", &a);
    scanf("%d", &b);
    scanf("%d", &n);
    
    double arr[n];
    double xs[n];
    double stepsize = ((double)(b-a))/(n-1); //assume n bigger than 1
    double first = (double)a;
    double mutbiggest = first*first*(sin(M_PI*first));
    
    for(int i=0; i<n; i++)
    {
      xs[i] = i*stepsize + a;
      arr[i] = xs[i]*xs[i]*(sin(M_PI*xs[i]));
      if (arr[i] > mutbiggest)
      {
        mutbiggest = arr[i];
      }
      else
      {
        mutbiggest = mutbiggest;
      }
    }
    double biggest = mutbiggest;
    
    printf("#biggest = %lf\n", biggest);
    /*for(int i=0; i<n; i++)
    {
      printf("%lf  %lf\n", xs[i], arr[i]);
    }*/
    
    
    
    
    
    FILE * temp = fopen("data.temp", "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;
    
    for (i=0; i < n; i++)
    {
      fprintf(temp, "%lf %lf \n", xs[i], arr[i]); //Write the data to a temporary file
    }

    fprintf(gnuplotPipe, "%s \n", "plot \"data.temp\" with lines"); //Send commands to gnuplot one by one.
}

