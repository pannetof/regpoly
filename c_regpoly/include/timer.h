
 
/* timer.h for ANSI C */

#ifndef timer_H_
#define timer_H_
 
typedef struct { unsigned long microsec, second; } timer_Chrono;




typedef enum {timer_sec, timer_min, timer_hours,
              timer_days, timer_hms} timer_TimeFormat;




void timer_Init (timer_Chrono * C);




double timer_Val (timer_Chrono C, timer_TimeFormat Unit);




void timer_Write (timer_Chrono C, timer_TimeFormat Unit);


 
#endif
 

