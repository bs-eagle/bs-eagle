/*!
  \file read_class.cpp
  \brief Read functions for Keyword Input Language class #FRead

  Class #FRead methods.
*/

#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#ifdef UNIX
#include <unistd.h>
#else
#include <io.h>
#endif
#include <fcntl.h>

#include "ascii_reader.h"


//#include "date_tools.h"
#include "main_def.h"
#include "timer.h"
#include "toupper.h"
//#include "localization.h"

/*!
  \brief Constructor -- make a new object of class FRead
*/
FRead::FRead ()
{
  fp = 0;
  top = 0;
  F = 0;
  F = new R_FILE[40];
  if (!F)
    fprintf (stderr, "Error: not enough memory!\n");
  incdir = 0;
  inc_num = 0;
  out_inc_list = 0;
}

/*!
  \brief Constructor -- make a new object of class FRead
*                       and open file 'fname'
  \param fname
  \param dir
*/
FRead::FRead (const char *fname, const char *dir)
{
  inc_num = 0;
  fp = 0;
  top = 0;
  F = 0;
  out_inc_list = 0;
  F = new R_FILE[40];
  if (!F)
    fprintf (stderr, "Error: not enough memory!\n");
  
  incdir = 0;
  if (set_incdir (dir) != 0)
    fprintf (stderr, "Error: not enough memory!\n");
  push (fname);
  
}

/*!
  \brief Method store path of main input file for work with include files
  \param dir
  \return
*/
int
FRead::set_incdir (const char *dir)
{
  t_long i, l = 0;
  // check pointer
  if (!dir)
    return -1;
  // allocate memory
  l = (t_long)strlen (dir) + 1;
  incdir = new char[l];
  if (!incdir)
    return -2;
  // copy
  strcpy (incdir, dir);
  // find first left '/' or '\'
  for (i = l - 1; i >= 0; --i)
    {
      if (incdir[i] == '\\' || incdir[i] == '/')
        break;
    }
  incdir[i + 1] = '\0';

  return 0;
}

/*!
  \brief Method return pointer to the string of include path
  \return pointer to the string of include path
*/
char *
FRead::get_incdir () const
{
  return incdir;
}

/*!
  \brief Destructor -- free all allocated memory and close all open files
*/
FRead::~FRead ()
{
  t_long i;

  // Pointer to current file - this is cpy of one of the pointers in F
  fp = 0;
  // Close all files
  if (F)
    {
      for (i = 0; i < top; ++i)
        F[i].close ();
    }
  top = 0;

  // Delete objects: NOTE destructor also close still open files
  if (F)
    delete[]F;
  F = 0;  
  if (incdir)
    delete[]incdir;
  incdir = 0;
  inc_list.clear_list ();
  out_inc_list = 0;
  inc_num = 0;
//  close_bin ();
}

/*!
  \brief close all opened files;
*/
void
FRead::close_all ()
{
  t_long i;
  fp = 0;
  // Close all files
  for (i = 0; i < top; ++i)
    F[i].close ();
  top = 0;
}

/*!
  \brief method try to open file with name 'fname', save this->fp in stack
  \param fname -- pointer to the name of file
  \param main_file_p

  \return if success                                       0
  \return if no name or can not open file                 -1
  \return if file already open                            -3
*/
int
FRead::push (const char *fname, int main_file_p)
{
  t_long i;
  char F_Buf[1024];

  if (fname == 0 || strlen (fname) == 0)
    {
      fprintf (stderr, "Error: input file name is empty\n");
      return -1;
    }

  if (incdir && !main_file_p)
    {
      // If full path spesified then do not prepend directory
#ifdef UNIX
      if (fname[0] != DIR_SYMBOL)
#else
      if (strlen (fname) <= 2 || (fname[1] != ':' && fname[0] != DIR_SYMBOL))
#endif
        {
          strcpy (F_Buf, incdir);
          strcat (F_Buf, fname);
        }
      else
        strcpy (F_Buf, fname);
    }
  else
    strcpy (F_Buf, fname);

  for (i = 0; i < top; ++i)
    {
      if (strcmp (F_Buf, F[i].get_name ()) == 0)
        {
          fprintf (stderr, "Error: in %s: the including file %s is already open...\n",
                   get_prefix (), F_Buf);
          return -3;
        }
    }
  if (!(fp = F[top].open (F_Buf)))
    {
      fprintf (stderr, "Error: in %s: cannot open included file %s\n",
               get_prefix (), F_Buf);
      return -3;
    }
  top++;
  
  // add file to the include files list
  if (inc_list.add_file (F_Buf))
    {
      fprintf (stderr, "Error: cannot add file to the iclude files list.\n");
      return -1;
    }
  // add file in to external include list
  if (out_inc_list && out_inc_list->add_file (F_Buf))
    {
      fprintf (stderr,   
               "cannot add file to the iclude files list.\n");
      return -1;
    }
  
  if (inc_num)
    printf ("Include file name     = %s\n", F_Buf);
  else
    {
      printf ("The following files will be included from %s.\n", F_Buf);
    }
  inc_num ++;

  return 0;
}

/*!
  \brief method close this->fp file and get open file from stack
  \return if success                                      0
  \return if no open files in stack                       -1
*/
int
FRead::pop ()
{
  if (top - 1 > 0)
    {
      F[top - 1].close ();
      F[top - 1].nstr = 0;
      --top;
      fp = F[top - 1].get_fp ();
      return 0;
    }
  else
    return -4;
}

/*!
  \brief Read line from this->fp, pass comments begining from '#' or ';',
         work with include files
  \param buf     Name of read buffer
  \param MaxLen  Maximum number of character to read
  \param flg     flag for case conversion if FREAD_CONVERT_CASE convert string to upper case
                                             FREAD_DONT_CONVERT_CASE do not convert string cases


  \return if success                                      Number of readed characters \n
   if no open files in FRead class                 NO_OPEN_FILE \n
   if 'MaxLen' parameter equals zero               0 \n
   if reading error                                -1 \n
   if end of file and nothing have been read       -2\n
   if can not open include file                    -5
*/
t_long
FRead::read_line (char *buf, const t_long MaxLen, const int flg)
{
  const char *inc =  ("include");
  const char *INC =  ("INCLUDE");
  const t_long inc_len = 7;	// length of word include
  const t_long f_len = 512;	// max length of file name
  char *word_cand;		// temporary buffer for reading            
  t_long word_len;
  char *fln;			// buffer for include file files name
  t_long disp;
  int current_error_num = 0;

  if (!buf)
    return -1;
  if (!this->fp)		// if pointer this->fp is bad return NO_OPEN_FILE
    return -5;
  if (MaxLen == 0)		// Check for max number of character is not a zero
    return 0;

  while (1)
    {
    beg:
      word_cand = buf;
      disp = get_next_not_empty_line (word_cand, MaxLen, this->fp, &word_len);
      SET_LINE_INFO (disp);

      // if exist INCLUDE keyword
      if ((strstr (word_cand, INC) == word_cand)
          || (strstr (word_cand, inc) == word_cand))
        {
          fln = trim_right_s (trim_left_s (word_cand + inc_len, '/'), '/');
          while (*fln == '\0')
            {
              fln = buf;
              disp = get_next_not_empty_line (fln, f_len, this->fp, &word_len);
              SET_LINE_INFO (disp);
              fln = trim_right_s (trim_left_s (fln, '/'), '/');
            }
          if ((current_error_num = this->push (fln)))	// try to open file with name <fln>
            return current_error_num;	// if can not open include file return -5
          continue;
        }
      // pass comment strings
    else if ((word_cand[0] != '#') && (word_cand[0] != ';') &&
	       !(word_cand[0] == '-' && word_cand[1] == '-'))
	{
	  word_cand = trim_right_s (word_cand);
	  word_len = (t_long) strlen(word_cand);
	  if (!((word_len == 1) && (*word_cand == '/')))
	    word_cand = trim_right_s (word_cand, '/');
	  break;
	}
    }
  if (flg == FREAD_CONVERT_CASE)
    locale_ucase (buf);
  return word_len;
}

/*!
  \brief Read text block up to "\", pass comments begining from '#' or ';' or '--'
  \param buf     Name of read buffer
  \param MaxLen  Maximum number of character to read
  \return if success                               Number of readed characters \n
   if no open files in FRead class                 NO_OPEN_FILE \n
   if 'MaxLen' parameter equals zero               0 \n
   if reading error                                -1 \n
   if end of file and nothing have been read       -2\n
*/
t_long
FRead::read_text_block (char *buf, const t_long MaxLen)
{
  char *word_cand;		// temporary buffer for reading            
  t_long word_len, len = 0;
  t_long disp;
  char *end_block;

  if (!buf)
    return -1;
  if (!this->fp)		// if pointer this->fp is bad return NO_OPEN_FILE
    return -5;
  if (MaxLen == 0)		// Check for max number of character is not a zero
    return 0;

  word_cand = buf;
  while (1)
    {
    beg:
      disp = get_next_not_empty_line (word_cand, MaxLen, this->fp, &word_len);
      SET_LINE_INFO (disp);

      // pass comment strings
      if ((word_cand[0] != '#') && (word_cand[0] != ';') &&
             !(word_cand[0] == '-' && word_cand[1] == '-'))
        {
          end_block = strchr (word_cand, '/');
          if (end_block != NULL)
            {
              t_long pos = (t_long)(end_block - word_cand);
              word_cand[pos] = '\0';
              len += pos;
              break;
            }
          
          word_cand[word_len] = ' ';
          word_cand = word_cand + word_len + 1;
          len += word_len + 1;
        }
    }

  return len;
}


/*!
  \brief   Recursive function\n
*          For keyword KEY read array from
*          file stream this->fp to buffer ARRAY.\n
*          string format: 2*15.8 3*{12 2*5.6 2*{1.1 1.2}} is equals\n
*          15.8 15.8 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2

  \param key    Name of calling keyword
  \param array  string buffer
  \param len_array  lenght of buffer
  \param pos number of doubles have been in array
  \param buf

  \return if success                                      number of read doubles\n
*         if bad pointer 'array'                          -1\n
*         if bad pointer 'buf'                            -2\n
*         if cann't allocate memory                       -3\n
*         if string format error                          -4
*/
t_long
FRead::convert_f (t_float *array, const t_long len_array, t_long pos, char *buf,
                  const char *key)
{
  char *sbuf;
  t_long c, i, counter;
  char *start_ptr, *end_ptr = 0;
  t_float t;
  // check section
  if (array == 0)               // check array pointer
    return -1;
  if (buf == 0)                 // check buf pointer
    return -2;
  if (pos >= len_array)         // check for input parameter
    return 0;

  sbuf = new char[strlen (buf) + 1];    // allocate new array
  if (sbuf == 0)                // check allocation error
    {
      fprintf (stderr, "Error: not enough memory!\n");
      return -3;
    }
  start_ptr = buf;              // set start pointer to begin of buf
  counter = 0;                  // set up counter
  
  // main loop
  for (;;)
    {
      // check for garbage
      if (pos >= len_array)
        {
          if (*end_ptr != '\0')
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
            }
          delete[]sbuf;
          return counter;
        }
      if (trim_left (&start_ptr))
        return -50;
      t = (t_float)strtod (start_ptr, &end_ptr); // try to read double from buf
      if (trim_left (&end_ptr))
        return -50;
      if (*start_ptr == '\0')
        {
          delete[]sbuf;
          return counter;
        }
      if (start_ptr == end_ptr) // if have not read return error -4
        {
          delete[]sbuf;
          return -4;
        }
      else if (*end_ptr == '*') // if next character is '*'
        {
          ++end_ptr;
          if (trim_left (&end_ptr))
            return -50;
          start_ptr = end_ptr;
          if (pos < len_array)
            array[pos] = (t_float)strtod (start_ptr, &end_ptr);      // try to read double from buf
          else
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
              delete[]sbuf;
              return counter;
            }
          c = pos;
          ++pos;
          ++counter;
          if (start_ptr == end_ptr) // if have not read return error -4
            {
              delete[]sbuf;
              return -4;
            }
          else
            {
              for (i = 0; i < (t_long) floor (t - 1 + 0.5); ++i)
                {
                  if (pos < len_array)
                    array[pos] = array[c];
                  else
                    {
                      fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                               get_prefix (), end_ptr, key);
                      delete[]sbuf;
                      return counter;
                    }
                  ++pos;
                  ++counter;
                }
            }
        }
      else
        {
          if (pos < len_array)
            array[pos] = t;
          else
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
              delete[]sbuf;
              return counter;
            }
          ++pos;
          ++counter;
        }
      start_ptr = end_ptr;
      if (*start_ptr == '\0')
        {
          delete[]sbuf;
          return counter;
        }
    }
//  delete[] sbuf;
//  return 0;
}
#if 0
/*!
  \brief   Recursive function\n
*          For keyword KEY read array from
*          file stream this->fp to buffer ARRAY.\n
*          string format: 2*15.8 3*{12 2*5.6 2*{1.1 1.2}} is equals\n
*          15.8 15.8 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2

  \param key    Name of calling keyword
  \param array  string buffer
  \param len_array  lenght of buffer
  \param pos number of doubles have been in array
  \param buf

  \return if success                                      number of read doubles\n
*         if bad pointer 'array'                          -1\n
*         if bad pointer 'buf'                            -2\n
*         if cann't allocate memory                       -3\n
*         if string format error                          -4
*/
long
FRead::convert_d (double *array, const long len_array, long pos, char *buf,
                  const char *key)
{
  char *sbuf;
  long c, i, counter;
  char *start_ptr, *end_ptr = 0;
  double t;
  // check section
  if (array == 0)               // check array pointer
    return -1;
  if (buf == 0)                 // check buf pointer
    return -2;
  if (pos >= len_array)         // check for input parameter
    return 0;

  sbuf = new char[strlen (buf) + 1];    // allocate new array
  if (sbuf == 0)                // check allocation error
    {
      fprintf (stderr, "not enough memory!\n");
      return -3;
    }
  start_ptr = buf;              // set start pointer to begin of buf
  counter = 0;                  // set up counter
  
  // main loop
  for (;;)
    {
      // check for garbage
      if (pos >= len_array)
        {
          if (*end_ptr != '\0')
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
            }
          delete[]sbuf;
          return counter;
        }
      if (trim_left (&start_ptr))
        return -50;
      t = strtod (start_ptr, &end_ptr); // try to read double from buf
      if (trim_left (&end_ptr))
        return -50;
      if (*start_ptr == '\0')
        {
          delete[]sbuf;
          return counter;
        }
      if (start_ptr == end_ptr) // if have not read return error -4
        {
          delete[]sbuf;
          return -4;
        }
      else if (*end_ptr == '*') // if next character is '*'
        {
          ++end_ptr;
          if (trim_left (&end_ptr))
            return -50;
          start_ptr = end_ptr;
          if (pos < len_array)
            array[pos] = strtod (start_ptr, &end_ptr);      // try to read double from buf
          else
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
              delete[]sbuf;
              return counter;
            }
          c = pos;
          ++pos;
          ++counter;
          if (start_ptr == end_ptr) // if have not read return error -4
            {
              delete[]sbuf;
              return -4;
            }
          else
            {
              for (i = 0; i < (long) floor (t - 1 + 0.5); ++i)
                {
                  if (pos < len_array)
                    array[pos] = array[c];
                  else
                    {
                      fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                               get_prefix (), end_ptr, key);
                      delete[]sbuf;
                      return counter;
                    }
                  ++pos;
                  ++counter;
                }
            }
        }
      else
        {
          if (pos < len_array)
            array[pos] = t;
          else
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
              delete[]sbuf;
              return counter;
            }
          ++pos;
          ++counter;
        }
      start_ptr = end_ptr;
      if (*start_ptr == '\0')
        {
          delete[]sbuf;
          return counter;
        }
    }
//  delete[] sbuf;
//  return 0;
}
#endif //0
/*!
  \brief   Recursiv function\n
*          For keyword KEY read array from
*          file stream this->fp to buffer ARRAY.\n
*          string format: 2*15.8 3*{12 2*5.6 2*{1.1 1.2}} is equals\n
*          15.8 15.8 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2

  \param key    Name of calling keyword
  \param array  string buffer
  \param len_array  lenght of buffer
  \param pos number of doubles have been in array
  \param buf

  \return if success                                      number of read doubles\n
*         if bad pointer 'array'                          -1\n
*         if bad pointer 'buf'                            -2\n
*         if cann't allocate memory                       -3\n
*         if string format error                          -4
*/

t_long
FRead::convert_u (t_int *array, const t_long len_array, t_long pos, char *buf,
                  const char *key)
{
  char *sbuf;
  t_long c, i, counter;
  char *start_ptr, *end_ptr = 0;
  t_int t;

  // check section
  if (array == 0)               // check array pointer
    return -1;
  if (buf == 0)                 // check buf pointer
    return -2;
  if (pos >= len_array)         // check for input parameter
    return 0;

  sbuf = new char[strlen (buf) + 1];    // allocate new array
  if (sbuf == 0)                // check allocation error
    {
      fprintf (stderr, "not enough memory!\n");
      return -3;
    }
  start_ptr = buf;              // set start pointer to begin of buf
  counter = 0;                  // set up counter
  for (;;)
    {
      if (pos >= len_array)
        {
          if (*end_ptr != '\0')
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
            }
          delete[]sbuf;
          return counter;
        }
      if (trim_left (&start_ptr))
        return -50;
      if (*start_ptr == '-')    // negative values found
        {
          delete[]sbuf;
          return -4;
        }
      t = (t_int)strtol (start_ptr, &end_ptr, 10);     // try to read int from buf
      if (trim_left (&end_ptr))
        return -50;
      if (start_ptr == end_ptr) // if have not read return error -4
        {
          delete[]sbuf;
          return -4;
        }
      else if (*end_ptr == '*') // if next character is '*'
        {
          ++end_ptr;
          if (trim_left (&end_ptr))
            return -50;
          start_ptr = end_ptr;
          if (*start_ptr == '-')    // negative values
            {
              delete[]sbuf;
              return -4;
            }
          if (pos < len_array)
            {
              array[pos] = (t_int)strtol (start_ptr, &end_ptr, 10);        // try to read int from buf
            }
          else
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
              delete[]sbuf;
              return counter;
            }
          c = pos;
          ++pos;
          ++counter;
          if (start_ptr == end_ptr) // if have not read return error -4
            {
              delete[]sbuf;
              return -4;
            }
          else
            {
              for (i = 0; i < t - 1; ++i)
                {
                  if (pos < len_array)
                    array[pos] = array[c];
                  else
                    {
                      fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                               get_prefix (), end_ptr, key);
                      delete[]sbuf;
                      return counter;
                    }
                  ++pos;
                  ++counter;
                }
            }
        }
      else
        {
          if (pos < len_array)
            array[pos] = t;
          else
            {
              fprintf (stderr, "Error: in %s: trailing garbage %s is ignored for keyword %s\n",
                       get_prefix (), end_ptr, key);
              delete[]sbuf;
              return counter;
            }
          ++pos;
          ++counter;
        }
      start_ptr = end_ptr;
      if (*start_ptr == '\0')
        {
          delete[]sbuf;
          return counter;
        }
    }
//  delete[] sbuf;
//  return 0;
}



/*!
  \brief unwrap constuction as '3*' to '* * *'
  \param s                       pointer to string
  \param flag -- if 0 - unwrap whole string, >0 unwrap after FLAG words
  \return if success                              0\n
*         if bad input pointer                    -2
*/
int
FRead::unwrap (char *s, const int flag) const
{
  char *ptr = 0;
  char *ptr_new = 0;
  char *ns_ptr = 0;
  char ns[CHAR_BUF_LEN];
  t_long n, i;
  t_long words_number = 0;

  if (!s)
    {
      return -2;
    }

  ns_ptr = ns;
  ptr = s;
  if (trim_left (&ptr))
    return -4;
  for (; *ptr != '\0' && words_number < flag;)
    {
      if (*ptr == ' ' || *ptr == '\t')
        {
          ++words_number;
          *ns_ptr = ' ';
          ++ns_ptr;
          if (trim_left (&ptr))
            return -4;
        }
      else
        {
          *ns_ptr = *ptr;
          ++ptr;
          ++ns_ptr;
        }

    }
  for (; *ptr != '\0';)         //main loop
    {
      n = strtol (ptr, &ptr_new, 10);   // try to get int
      if (ptr != ptr_new)       // if something read
        {
          if (*ptr_new == '*')  // check for '*'
            {
              *ns_ptr = ' ';
              ++ns_ptr;
              for (i = 0; i < n; ++i)   // copy '*' to new array
                {
                  *ns_ptr = '*';
                  ++ns_ptr;
                  *ns_ptr = ' ';
                  ++ns_ptr;
                }
              ++ptr_new;
              ptr = ptr_new;
            }
          else                  // copy all from ptr to ptr_new
            {
              for (; ptr != ptr_new; ++ptr, ++ns_ptr)
                *ns_ptr = *ptr;
            }
        }
      else
        {
          *ns_ptr = *ptr;
          ++ptr;
          ++ns_ptr;
        }
    }

  *ns_ptr = ' ';
  ++ns_ptr;
  *ns_ptr = '/';
  ++ns_ptr;
  *ns_ptr = '\0';
  strcpy (s, ns);
  return 0;
}

/*!
  \brief try to read string to dest if * read nothing
  \param start_ptr                  pointer to string
  \param end_ptr                    pointer to the next value
  \param dest                       read string
  \param no_default_p               signal error instead of use default values

  \return if success                              0\n
         if cannot allocate memory               -1\n
         if bad input pointer                    -2\n
         if trim_left return error  code         -3
*/
int
FRead::scanf_s (char *start_ptr, char **end_ptr, char *dest, int no_default_p, int *is_default) const
{
  char *dest_ptr = 0;
  char *ptr = 0;
  int current_error_num;
  //check all
  if (!start_ptr || !dest)
    {
      return -2;
    }
  if ((current_error_num = trim_left (&start_ptr)) != 0)
    {
      return current_error_num;
    }
  if (is_default)
    *is_default = 0;
  // If no default values, then read '*' literally
  if (*start_ptr == '*' && !no_default_p)
    {
      *end_ptr = start_ptr;
      ++(*end_ptr);
      if (is_default)
        *is_default = 1;
      return 0;
    }
  if (*start_ptr == '/')
    {
      *end_ptr = start_ptr;
      if (no_default_p)
        return -14;
      if (is_default)
        *is_default = 1;
      return 0;
    }
  *end_ptr = start_ptr;
  while (**end_ptr == '\'' || **end_ptr == '\"')
    (*end_ptr)++;
      
  dest_ptr = dest;
  for (; **end_ptr != ' ' && **end_ptr != '\0' && **end_ptr != '\t';
       ++(*end_ptr), ++dest_ptr)
    *dest_ptr = **end_ptr;
  
  ptr = *end_ptr - 1;
  // move to last not '\0' or ' ' or '\t'
  dest_ptr--;
  while (*ptr == '\'' || *ptr == '\"')
    {
      --ptr;
      --dest_ptr;
    }  
  *(dest_ptr + 1) = '\0';
  return 0;
}

//read word in text with commas
int
FRead::scanf_text (char *start_ptr, char **end_ptr, char *dest) const
{
  char *dest_ptr = 0;
  int current_error_num;
  bool is_commas = false;

  //check all
  if (!start_ptr || !dest)
    {
      return -2;
    }
  if ((current_error_num = trim_left (&start_ptr)) != 0)
    {
      return current_error_num;
    }

  if (strlen (start_ptr) == 0)
    return -14;

  *end_ptr = start_ptr;
  if (**end_ptr == '\'' || **end_ptr == '\"')
    {
      (*end_ptr)++;
      is_commas = true;
    }
      
  dest_ptr = dest;
  for (; **end_ptr != '\0' && **end_ptr != '\t'; ++(*end_ptr), ++dest_ptr)
    {
      if (is_commas)
        {
          if (**end_ptr == '\'' || **end_ptr == '\"')
            {
              ++(*end_ptr);
              break;
            }
        }
      else
        {  
          if (**end_ptr == ' ')
            break;
        }

      *dest_ptr = **end_ptr;
    }
  
  dest_ptr--; 
  *(dest_ptr + 1) = '\0';
  return 0;
}

/*!
  \brief try to read integer to dest if * read nothing
  \param start_ptr                  pointer to string
  \param end_ptr                    pointer to the next value
  \param dest                       int

  \return if success                              0\n
         if cannot allocate memory               -1\n
         if bad input pointer                    -2\n
         if trim_left return error  code         -3
*/
int
FRead::scanf_u (char *start_ptr, char **end_ptr, t_long *dest, int *is_default) const
{
  t_long t;
  int current_error_num;
  //check all
  if (!start_ptr || !dest)
    {
      return -1;
    }
  if ((current_error_num = trim_left (&start_ptr)) != 0)
    {
      return current_error_num;
    }
  if (is_default)
    *is_default = 0;
  if (*start_ptr == '*')
    {
      *end_ptr = start_ptr;
      ++(*end_ptr);
      if (is_default)
        *is_default = 1;
      return 0;

    }
  if (*start_ptr == '/')
    {
      *end_ptr = start_ptr;
      if (is_default)
        *is_default = 1;
      return 0;
    }
  t = (t_long)strtol (start_ptr, end_ptr, 10);
  if (start_ptr == *end_ptr)
    {
      // Cannot read integer: error output should be done by caller
      return -3;
    }
  *dest = t;
  return 0;
}

/*!
  \brief try to read double to dest if * read nothing
  \param start_ptr                  pointer to string
  \param end_ptr                    pointer to the next value
  \param dest                       double

  \return if success                              0\n
         if bad input pointer                    -2\n
         if trim_left return error  code         current_error_num
*/
int
FRead::scanf_d (char *start_ptr, char **end_ptr, t_double *dest, int no_default_p, int *is_default) const
{
  double t;
  int current_error_num;
  //check input pointers
  if (!start_ptr || !dest)
    {
      return -1;
    }
  // Try to trim all left blanks and tab symbols
  if ((current_error_num = trim_left (&start_ptr)) != 0)
    {
      return current_error_num;
    }
  if (is_default)
    *is_default = 0;
  // if '*' have found return 0 and don't change value *dest
  // also move pointer to the next symbol
  if (*start_ptr == '*' && !no_default_p)
    {
      *end_ptr = start_ptr;
      ++(*end_ptr);
      if (is_default)
        *is_default = 1;
      return 0;

    }
  // if '*' have found return 0 and don't change value *dest
  // don't move pointer

  if (*start_ptr == '/')
    {
      *end_ptr = start_ptr;
      if (no_default_p)
        return -14;
      if (is_default)
        *is_default = 1;
      return 0;
    }
  // try to read double value from input string
  t = (t_double)strtod (start_ptr, end_ptr);
  // check reading result
  // if have read nothing return -3
  if (start_ptr == *end_ptr)
    {
      // Cannot read double: error output should be done by caller
      return -3;
    }
  *dest = t;
  return 0;
}

/*!
  \brief build and return prefix
  \return Pointer to prefix string
*/
char *
FRead::get_prefix ()
{
  t_long i;
  char s[512];

  prefix[0] = '\0';
  for (i = 0; i < top; ++i)
    {
      sprintf (s, "%s '%s' %s %ld ",  ("File"), F[i].get_name (),
                ("Line"), (long)F[i].nstr);
      strcat (prefix, s);
    }
  return prefix;
}

/*!
  \brief read file name from source
  \param dest buffer for file name
  \param source buffer to read from
  \return 0 if success < 0 else
*/
int
FRead::scanf_file_name (char *dest, char *source) const
{
  char *d_ptr;
  char *s_ptr = source;

  if (!dest || !source)
    return -1;

  // skip left blanks
  trim_left (&s_ptr);

  if ((*s_ptr)== '\"' || (*s_ptr)== '\'')
    {
      ++s_ptr;
    }
  for (d_ptr = dest; (*s_ptr) != '\"' && (*s_ptr) != '\'' && (*s_ptr) != '\0'; ++s_ptr, ++d_ptr)
    *d_ptr = *s_ptr;
  *d_ptr = '\0';
  return 0;
}
/*!
  \brief skip given number of double arguments
  \param start_ptr -- pointer to the initial position in string
  \param end_ptr -- final position in string
  \param count -- number of skiping arguments
  \return 0 if success < 0 if error occur 
*/
int
FRead::skip_d (char *start_ptr, char **end_ptr, const t_long count) const 
{
  int current_error_num;
  t_long i;
  
  //check input pointers
  if (!start_ptr || !end_ptr)
    {
      return -1;
    }
  // main loop
  for (i = count; i > 0; --i)
    {
      // Try to trim all left blanks and tab symbols
      if ((current_error_num = trim_left (&start_ptr)) != 0)
        {
          return current_error_num;
        }
      // check for default argument if found continue the loop  
      if (*start_ptr == '*')
        {
          ++(start_ptr);
          continue;
        }
      // if found final default argument set final pointer and return success
      if (*start_ptr == '/')
        {
          *end_ptr = start_ptr;
          return 0;
        }
      // check for the end of string
      if (*start_ptr == '\0') 
        {
          *end_ptr = start_ptr;
          return -1;
        }
      // try to read double value from input string
      double t = strtod (start_ptr, end_ptr);

      if (start_ptr == *end_ptr)
        {
          return -1 * (int)t * 0 - 1;  // just make compiller happy (unused t)
        }
      start_ptr = *end_ptr;
    }
  // set final position in string
  *end_ptr = start_ptr;
  return 0;
}
/*!
  \brief skip given number of int arguments
  \param start_ptr -- pointer to the initial position in string
  \param end_ptr -- final position in string
  \param count -- number of skiping arguments
  \return 0 if success < 0 if error occur 
*/
int
FRead::skip_u (char *start_ptr, char **end_ptr, const t_long count) const
{
  t_long i;
  t_long buf;

  for (i = count; i > 0; --i)
    {
      if (scanf_u (start_ptr, end_ptr, &buf))
        return -1;
      start_ptr = *end_ptr;
    }

  return 0;
}

/*!
  \brief skip given number of string arguments
  \param start_ptr -- pointer to the initial position in string
  \param end_ptr -- final position in string
  \param count -- number of skiping arguments
  \return 0 if success < 0 if error occur 
*/
int
FRead::skip_s (char *start_ptr, char **end_ptr, const t_long count) const
{
  int i;
  char buf[CHAR_BUF_LEN];

  for (i = count; i > 0; --i)
    {
      if (scanf_s (start_ptr, end_ptr, buf))
        return -1;
      start_ptr = *end_ptr;
    }

  return 0;
}

