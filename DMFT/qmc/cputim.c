#include <sys/types.h>
#include <sys/times.h>
#include <stdio.h>

cputim_(dsec)
double *dsec;
{
      struct tms buffer;

      times(&buffer);
      *dsec = (double)buffer.tms_utime/100.0;
/*      fprintf(stderr,"inside cputim: %lf\n",dsec);*/
}
