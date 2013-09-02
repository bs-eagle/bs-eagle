#include <stdio.h>
#include <string>
#include <time.h>

void int2str (int a, char *str)
{
if (a < 10)
  sprintf(str, "0%d", a);
else
  sprintf(str, "%d", a);
}


// return date time as string in format dd.MM.yyyy hh:mm:ss
std::string get_date_time_str ()
{
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  char str[1024];
  int2str(timeinfo->tm_mday, &str[0]);
  int2str(1 + timeinfo->tm_mon, &str[3]);
  int2str(1900 + timeinfo->tm_year, &str[6]);
  int2str(timeinfo->tm_hour, &str[11]);
  int2str(timeinfo->tm_min, &str[14]);
  int2str(timeinfo->tm_sec, &str[17]);
  str[2] = str[5] = '.';
  str[10] = '_';
  str[13] = str[16] = ':';

  return std::string (str);
}


