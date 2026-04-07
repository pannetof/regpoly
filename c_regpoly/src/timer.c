/*
 * timer.c — CPU time measurement module.
 * Authors: Pierre L'Ecuyer, Benoit Martin.
 *          Adapted to ANSI C by Francis Picard and Jean-Sebastien Senecal, August 1999.
 */

#include <sys/times.h>
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef CLK_TCK
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif
#include "timer.h"

/* Returns the CPU time consumed by the program since its start,
   split into whole seconds (*tsec) and microseconds (*tusec).
   Uses the POSIX times() system call (not ANSI C). */
static void Heure(unsigned long *tsec, unsigned long *tusec)
{
  clock_t y;
  struct tms us;

  y = times(&us);
  *tusec = us.tms_utime + us.tms_stime;  /* CPU time = user time + system time */
  *tsec = *tusec / CLK_TCK;
  *tusec = (*tusec % CLK_TCK) * 1000000 / CLK_TCK;
}

/* Initializes (resets) the stopwatch C to the current CPU time. */
void timer_Init(timer_Chrono *C)
{
  Heure(&C->second, &C->microsec);
}

/* Returns the CPU time elapsed since the last timer_Init call for C,
   expressed in the unit given by Unit. Passing timer_hms is an error. */
double timer_Val(timer_Chrono C, timer_TimeFormat Unit)
{
  double temps;
  timer_Chrono now;

  Heure(&now.second, &now.microsec);
  temps = (((double) now.microsec - (double) C.microsec) / 1.0E+6
           + (double) now.second) - (double) C.second;
  switch (Unit) {
  case timer_sec:
    return temps;
  case timer_min:
    return temps * 1.666666667E-2;
  case timer_hours:
    return temps * 2.777777778E-4;
  case timer_days:
    return temps * 1.157407407E-5;
  case timer_hms:
    printf("timer_Val: timer_hms is not a valid unit for timer_Val\n");
    exit(1);
  }
  return 0.0;
}

/* Prints the CPU time elapsed since the last timer_Init call for C,
   formatted according to Form. Use timer_hms for HH:MM:SS.xx output. */
void timer_Write(timer_Chrono C, timer_TimeFormat Form)
{
  long centieme, minute, heure, seconde;
  double temps;

  if (Form != timer_hms)
    temps = timer_Val(C, Form);
  else
    temps = 0.0;
  switch (Form) {
  case timer_sec:
    printf("%10.2lf seconds", temps);
    break;
  case timer_min:
    printf("%10.2lf minutes", temps);
    break;
  case timer_hours:
    printf("%10.2lf hours", temps);
    break;
  case timer_days:
    printf("%10.2lf days", temps);
    break;
  case timer_hms:
    temps = timer_Val(C, timer_sec);
    heure = (long) (temps * 2.777777778E-4);
    if (heure > 0)
      temps -= (double) heure * 3600.0;
    minute = (long) (temps * 1.666666667E-2);
    if (minute > 0)
      temps -= (double) minute * 60.0;
    seconde = (long) temps;
    centieme = (long) (100.0 * (temps - (double) seconde));
    printf("%02ld:%02ld:%02ld.%02ld", heure, minute, seconde, centieme);
    break;
  }
}
