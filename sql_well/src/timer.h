#ifndef __TIMER_H
#define __TIMER_H

/*!
  \file timer.h
  \brief Time functions declarations. 
*/

typedef struct _timer_
{
  int tic;
  int sec;
  int min;
  int hour;
}
TIMER;

#include <time.h>

void ConvertTime (long clocks, TIMER * t);

// Return string representing elapsed time from the first call
// Argument is buffer for answer, return - answer string (this buffer)
char * get_elapsed_time_as_string (long elapsed_time, char *string);
//! Start timer
void start_elapsed_timer (void);
//! Start thread timer
void start_elapsed_thread_timer ();

//! Get time from the previous call
long get_elapsed_time (void);
//! Get thread time from the previous call
long get_elapsed_thread_time ();

void get_elapsed_time_as_long (long *elapsed_time = 0);

// Returns astronomic time in milliseconds
long int get_full_time (void);

//! Whether current time greater then 31.08.2007 = 1188500400 secs from 01.01.1970

#if 0
#define CHECK_TIME_LIMIT() \
  (get_full_time () > 1188500400)
#endif //0

/*!
  \fn long int get_thread_time ()
  \brief Return time in one hundreds of second.\n
         The is user time + system time
  \return Return time in one hundreds of second.
*/
long int get_thread_time (void);

/*!
  \fn long int get_process_time ()
  \brief Return time in one hundreds of second.\n
         The is only user time, system time does not added
  \return Return time in one hundreds of second.
*/
long int get_process_time (void);

#endif // __TIMER_H
