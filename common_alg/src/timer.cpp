/*!
  \file timer.cpp
  \brief timer - output package
*/

#include <stdio.h>
#include <time.h>

#include "timer.h"

// The way for processor time measurement
#ifdef UNIX
//! All three ususally supported in UNIX systems, but the first is prefferd
#define USE_getrusage
//#define USE_times
//#define USE_clock
//! The way for astronomic time measurement
#define USE_gettimeofday
#else //UNIX
//! Only clock supported in all Windows systems, but measures astronomic time
//! instead of processor one.
//#define USE_clock
// This way supported all NT clones, but not Win9x/ME
#define USE_GetProcessTimes
//! The way for astronomic time measurement
#define USE_time
#endif //UNIX



#ifdef USE_getrusage
#include <sys/time.h>
#include <sys/resource.h>

#ifdef __hpux
#include <sys/syscall.h>
#define getrusage(a,b) syscall(SYS_getrusage,a,b)
#endif //__hpux

/*!
  \fn long int get_thread_time ()
  \brief Return time in one hundreds of second.\n
         The is only user time, system time does not added.
  \return Return time in one hundreds of second.
*/
long int
get_thread_time ()
{
  struct rusage buf;

  getrusage (RUSAGE_SELF, &buf);
  return buf.ru_utime.tv_sec * 100      // convert time in seconsds
    // into one hundreds
    + buf.ru_utime.tv_usec / 10000;     // convert time in microseconds
  // into one hundreds
}
#endif /* USE_getrusage */

#ifdef USE_times
#include <time.h>
#include <sys/times.h>

/*!
  \fn long int get_thread_time ()
  \brief Return time in one hundreds of second.\n
         The is only user time, system time does not added
  \return Return time in one hundreds of second.
*/
long int
get_thread_time ()
{
  struct tms buf;

  times (&buf);

  return buf.tms_utime / (CLK_TCK / 100);       // convert time in CLK_TCK
  // into one hundreds
}
#endif /* USE_times */

#ifdef USE_clock
#include <time.h>

/*!
  \fn long int get_thread_time ()
  \brief Return time in one hundreds of second.\n
         The is user time + system time
  \return Return time in one hundreds of second.
*/
long int
get_thread_time ()
{
  long int t;

  t = (long int) clock ();

  return t / (CLOCKS_PER_SEC / 100);    // convert time in CLOCKS_PER_SEC
  // into one hundreds
}
#endif /* USE_clock */

#ifdef USE_gettimeofday
#include <sys/time.h>
#include <time.h>

/*!
  \fn long int get_full_time ()
  \brief Returns astronomic time in seconds
  \return Returns astronomic time in seconds

*/

long int
get_full_time ()
{
  struct timeval buf;

#if 0                           // Compute 31.08.2005 in secs from 01.01.1970
  static struct tm zero;
  struct tm t;

  t = zero;
  t.tm_mday = 31;
  t.tm_mon = 8 - 1;
  t.tm_year = 105;
  return mktime (&t);
#endif

  gettimeofday (&buf, 0);
#if 1
  return buf.tv_sec;            // return time in seconds
#else
  // convert seconds into milliseconds
  return buf.tv_sec * 100
    // convert microseconds into seconds
    + buf.tv_usec / 10000;
#endif
}

/*!
  \fn long int get_process_time ()
  \brief Return time in one hundreds of second.\n
         The is only user time, system time does not added
  \return Return time in one hundreds of second.
*/
long int
get_process_time ()
{
  struct timeval buf;

  gettimeofday (&buf, 0);
  // convert seconds into milliseconds
  return buf.tv_sec * 100
    // convert microseconds into seconds
    + buf.tv_usec / 10000;
}

#endif /* USE_gettimeofday */

#ifdef USE_time
#include <time.h>

/*!
  \fn long int get_full_time ()
  \brief Returns astronomic time in seconds
*/
long int
get_full_time ()
{
  return (long int)time (0);
}
#endif /* USE_time */


#ifndef UNIX
#ifdef USE_GetProcessTimes
#include "Windows.h"
#ifdef _OPENMP  
#include <omp.h>
#endif //_OPENMP

#ifdef _MPI
#include <mpi.h>
#endif //MPI

long int
get_process_time ()
{
  long int t;
#ifdef _OPENMP
  t = (long)(omp_get_wtime() * 100);
#else //_OPENMP
#ifdef _MPI
  t = (long)(MPI_Wtime () * 100);
#else //_MPI

  __int64 CreationTime, ExitTime, KernelTime, UserTime;
  if (GetProcessTimes
      (GetCurrentProcess (), (FILETIME *) & CreationTime,
       (FILETIME *) & ExitTime, (FILETIME *) & KernelTime,
       (FILETIME *) & UserTime) == 0)
    {
      // error, use clock as last resort
      return ((long int) clock ()) / (CLOCKS_PER_SEC / 100);
    }
  t = (long int) ((KernelTime + UserTime) / 100000);
#endif //_MPI
#endif //_OPENMP
  return t;
}

long int
get_thread_time ()
{
  long int t;
  __int64 CreationTime, ExitTime, KernelTime, UserTime;
  if (GetThreadTimes (GetCurrentThread (), (FILETIME *) & CreationTime,
       (FILETIME *) & ExitTime, (FILETIME *) & KernelTime,
       (FILETIME *) & UserTime) == 0)
    {
      // error, use clock as last resort
      return ((long int) clock ()) / (CLOCKS_PER_SEC / 100);
    }
  t = (long int) ((KernelTime + UserTime) / 100000);
  return t;
}
#endif // USE_GetProcessTimes
#endif // !UNIX

/*!
  \fn static void ConvertTime (long clocks, TIMER *t)
  \brief Convert time in one hundreds into hh.mm.ss:tt
  \param clocks
  \param t
*/
void
ConvertTime (long clocks, TIMER * t)
{
  t->hour = clocks / 360000L;
  clocks %= 360000L;
  t->min = clocks / 6000;
  clocks %= 6000;
  t->sec = clocks / 100;
  t->tic = clocks % 100;
}

/*!
  \fn char * get_elapsed_time_as_string (long elapsed_time, char *string)
  \brief Return string representing elapsed time from the first call\n
         Argument is buffer for answer, return - answer string (this buffer)\n
         To start times pass null pointer as string.
  \param elapsed_time -- returned by get_elapsed_time
  \param string -- buffer for answer
  \return answer string (this buffer)
*/
char *
get_elapsed_time_as_string (long elapsed_time, char *string)
{
  TIMER summ;

  if (!string)
    return 0;

  ConvertTime (elapsed_time, &summ);
  sprintf (string, "%2.2d.%2.2d.%2.2d",
           summ.hour, summ.min, summ.sec);
  return string;
}

static int TimerStarted = 0;    // 1 if timer started
static long int StartTime = 0;  // start time

static int ThreadTimerStarted = 0;      // 1 if thread timer started
static long int ThreadStartTime = 0;    // thread start time
//! Start timer
void
start_elapsed_timer ()
{
  TimerStarted = 1;
  StartTime = get_process_time ();
}

//! Start thread timer
void
start_elapsed_thread_timer ()
{
  ThreadTimerStarted = 1;
  ThreadStartTime = get_thread_time ();
}

//! Get time from the previous call
long
get_elapsed_time ()
{
  long t;
  t = get_process_time ();
  if (!TimerStarted)
    {
      TimerStarted = 1;
      StartTime = t;
    }
  return t - StartTime;
}

//! Get thread time from the previous call
long
get_elapsed_thread_time ()
{
  long t;
  t = get_thread_time ();
  if (!TimerStarted)
    {
      ThreadTimerStarted = 1;
      ThreadStartTime = t;
    }
  return t - ThreadStartTime;
}


