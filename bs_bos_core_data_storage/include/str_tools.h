/** 
 * \file str_tools.h
 * \brief tools for string processing for FRead
 * \author Sergey Miryanov
 * \date 25.05.2009
 * */
#ifndef BS_BOS_CORE_DATA_STORAGE_STR_TOOLS_H_
#define BS_BOS_CORE_DATA_STORAGE_STR_TOOLS_H_

#include <boost/format.hpp>

namespace blue_sky {

// trim all right blanks
  inline char *
  trim_left_s (char *s, const char ch = ' ')
  {
    size_t str_pos, length;

    if (s == NULL)
      return NULL;
    if ((length = strlen (s)) == 0)
      return s;

    str_pos = 0;
    while (s[str_pos] == ' ' ||
           s[str_pos] == '\t' || s[str_pos] == '\'' || s[str_pos] == ch)
      {
        ++str_pos;
        if (length - 1 < str_pos)
          break;
      }
    return (s + str_pos);
  }

// trim all left blanks
  inline char *
  trim_right_s (char *s, const char ch = ' ')
  {
    int str_pos, length;

    if (s == NULL)
      return NULL;
    if ((length = (int)strlen (s)) == 0)
      return s;

    str_pos = length - 1;
    while (s[str_pos] == ' ' ||
           s[str_pos] == '\t' ||
           s[str_pos] == '"' || s[str_pos] == '\'' || s[str_pos] == ch)
      {
        s[str_pos--] = '\0';
        if (str_pos < 0)
          break;
      }
    return s;
  }

  inline int
  get_next_not_empty_line (char *string, int MaxLen, FILE * fp, int *len)
  {
    char *ptr_s = string;
    int ch;
    int disp = 0;

    if ((ptr_s == NULL) || (fp < 0))
      return -1;

    while ((ch = fgetc (fp)) != EOF)
      {
        if (ch == '\n' || ch == '\r' || ch == char (10) || ch == char (13))
          ++disp;
        else if (!(ch == ' ' || ch == '\t'))
          break;
      }

    if (ch == EOF)
      return -1;

    *len = 0;
    while (--MaxLen > 0)
      {
        *ptr_s++ = ch;
        ++(*len);
        if ((ch = fgetc (fp)) == EOF || (ch == '\n') || (ch == char (10)) || (ch == char (13)))
          break;
      }
    *ptr_s = '\0';
    return (ptr_s == string) ? -1 : ++disp;
  }

  /*!
    \brief change pointer to the next not blanks
    \param start_ptr -- pointer to pointer to string
    \return if success  number of read doubles\n
  *         if bad pointer '*start_ptr' -1
  */
  inline void
  trim_left (char **start_ptr)
  {
    if (!(*start_ptr))
      {
        bs_throw_exception ("Invalid argument: start_ptr points to null");
      }

    while ((**start_ptr == '\t' 
            || **start_ptr == ' ' 
            || **start_ptr == ','
            || **start_ptr == ':') && (**start_ptr != '\0'))
      {
        ++(*start_ptr);
      }
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
  inline void
  scanf_s (char *start_ptr, char **end_ptr, char *dest, int no_default_p = 0, int *is_default = 0)
  {
    char *dest_ptr = 0;
    char *ptr = 0;
    //check all
    if (!start_ptr || !dest)
      {
        bs_throw_exception ((boost::format ("Invalid argument: start_ptr or dest is null [%p, %p]") % start_ptr % dest).str ());
      }
    trim_left (&start_ptr);

    if (is_default)
      *is_default = 0;
    // If no default values, then read '*' literally
    if (*start_ptr == '*' && !no_default_p)
      {
        *end_ptr = start_ptr;
        ++(*end_ptr);
        if (is_default)
          {
            *is_default = 1;
          }
      }
    if (*start_ptr == '/')
      {
        *end_ptr = start_ptr;
        if (no_default_p)
          bs_throw_exception ("no_default_p is true");

        if (is_default)
          {
            *is_default = 1;
          }
      }
    *end_ptr = start_ptr;
    while (**end_ptr == '\'' || **end_ptr == '\"')
      {
        (*end_ptr)++;
      }

    dest_ptr = dest;
    for (; **end_ptr != ' ' && **end_ptr != '\0' && **end_ptr != '\t';
         ++(*end_ptr), ++dest_ptr)
      {
        *dest_ptr = **end_ptr;
      }

    ptr = *end_ptr - 1;
    // move to last not '\0' or ' ' or '\t'
    dest_ptr--;
    while (*ptr == '\'' || *ptr == '\"')
      {
        --ptr;
        --dest_ptr;
      }
    *(dest_ptr + 1) = '\0';
  }


  /*!
    \brief skip given number of double arguments
    \param start_ptr -- pointer to the initial position in string
    \param end_ptr -- final position in string
    \param count -- number of skiping arguments
    \return 0 if success < 0 if error occur
  */
  inline void
  skip_d (char *start_ptr, char **end_ptr, int count)
  {
    double t;
    int i;

    //check input pointers
    if (!start_ptr || !end_ptr)
      {
        bs_throw_exception ((boost::format ("Invalid argument start_ptr or end_ptr is null [%p, %p]") % start_ptr % end_ptr).str ());
      }

    // main loop
    for (i = count; i > 0; --i)
      {
        // Try to trim all left blanks and tab symbols
        trim_left (&start_ptr);

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
            return ;
          }
        // check for the end of string
        if (*start_ptr == '\0')
          {
            *end_ptr = start_ptr;
            bs_throw_exception ("Reach end of string");
          }
        // try to read double value from input string
        t = (float)strtod (start_ptr, end_ptr);
        if (start_ptr == *end_ptr)
          {
            bs_throw_exception ("Empty string");
          }
      }
    // set final position in string
    *end_ptr = start_ptr;
  }

  //read word in text with commas
  inline void 
  scanf_text (char *start_ptr, char **end_ptr, char *dest)
  {
    char *dest_ptr = 0;
    bool is_commas = false;

    //check all
    if (!start_ptr || !dest)
      {
        bs_throw_exception ((boost::format ("Invalid argument: start_ptr or dest is null [%p, %p]") % start_ptr % dest).str ());
      }

    trim_left (&start_ptr);
    if (strlen (start_ptr) == 0)
      bs_throw_exception ("start_ptr is empty");

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
  template <typename T>
  inline void
  scanf_u (char *start_ptr, char **end_ptr, T *dest, int *is_default = 0)
  {
    //check all
    if (!start_ptr || !dest)
      {
        bs_throw_exception ((boost::format ("Invalid argument: start_ptr or dest is null [%p, %p]") % start_ptr % dest).str ());
      }

    trim_left (&start_ptr);
    if (is_default)
      *is_default = 0;
    if (*start_ptr == '*')
      {
        *end_ptr = start_ptr;
        ++(*end_ptr);
        if (is_default)
          *is_default = 1;

        return ;
      }
    if (*start_ptr == '/')
      {
        *end_ptr = start_ptr;
        if (is_default)
          *is_default = 1;

        return ;
      }
    int t = strtol (start_ptr, end_ptr, 10);
    if (start_ptr == *end_ptr)
      {
        // Cannot read integer: error output should be done by caller
        bs_throw_exception ("Can't read int");
      }
    *dest = t;
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
  inline void 
  scanf_d (char *start_ptr, char **end_ptr, double *dest, int no_default_p = 0, int *is_default = 0, const char *exception_context = "")
  {
    //check input pointers
    if (!start_ptr || !dest)
      {
        bs_throw_exception (boost::format ("%s: Invalid argument: start_ptr or dest is null [%p, %p]") % exception_context % start_ptr % dest);
      }

    // Try to trim all left blanks and tab symbols
    trim_left (&start_ptr);
    if (is_default)
      *is_default = 0;

    // if '*' have found return YS_SUCCESS and don't change value *dest
    // also move pointer to the next symbol
    if (*start_ptr == '*' && !no_default_p)
      {
        *end_ptr = start_ptr;
        ++(*end_ptr);
        if (is_default)
          *is_default = 1;

        return ;
      }

    // if '*' have found return YS_SUCCESS and don't change value *dest
    // don't move pointer
    if (*start_ptr == '/')
      {
        *end_ptr = start_ptr;
        if (no_default_p)
          bs_throw_exception (boost::format ("%s: End of string reach and no default value") % exception_context);

        if (is_default)
          *is_default = 1;

        return ;
      }

    // try to read double value from input string
    double t = (float)strtod (start_ptr, end_ptr);
    // check reading result
    // if have read nothing return -3
    if (start_ptr == *end_ptr)
      {
        bs_throw_exception (boost::format ("%s: Can't read double") % exception_context);
      }
    *dest = t;
  }

  /*!
    \brief read file name from source
    \param dest buffer for file name
    \param source buffer to read from
    \return 0 if success < 0 else
  */
  inline void 
  scanf_file_name (char *dest, char *source)
  {
    char *d_ptr;
    char *s_ptr = source;
    int flag = 0;

    if (!dest || !source)
      {
        bs_throw_exception ((boost::format ("Invalid argument: dest or source is null [%p, %p]") % dest % source).str ());
      }

    // skip left blanks
    trim_left (&s_ptr);
    if ((*s_ptr)== '\"' || (*s_ptr)== '\'')
      {
        flag = 1;
        ++s_ptr;
      }

    for (d_ptr = dest; (*s_ptr) != '\"' && (*s_ptr) != '\'' && (*s_ptr) != '\0'; ++s_ptr, ++d_ptr)
      {
        *d_ptr = *s_ptr;
      }

    *d_ptr = '\0';
  }

  /*!
    \brief skip given number of int arguments
    \param start_ptr -- pointer to the initial position in string
    \param end_ptr -- final position in string
    \param count -- number of skiping arguments
    \return 0 if success < 0 if error occur
  */
  inline void 
  skip_u (char *start_ptr, char **end_ptr, int count)
  {
    int buf;

    if (!start_ptr || !end_ptr)
      {
        bs_throw_exception ((boost::format ("Invalid argument: dest or source is null [%p, %p]") % start_ptr % end_ptr).str ());
      }

    for (int i = count; i > 0; --i)
      {
        scanf_u (start_ptr, end_ptr, &buf);
        start_ptr = *end_ptr;
      }
  }

  /*!
    \brief unwrap constuction as '3*' to '* * *'
    \param s                       pointer to string
    \param flag -- if 0 - unwrap whole string, >0 unwrap after FLAG words
    \return if success                              0\n
  *         if bad input pointer                    -2
  */
  inline void 
  unwrap (char *s, const int flag = 0)
  {
    char *ptr = 0;
    char *ptr_new = 0;
    char *ns_ptr = 0;
    char ns[CHAR_BUF_LEN];
    int n, i;
    int words_number = 0;

    if (!s)
      {
        bs_throw_exception ((boost::format ("Invalid argument: s is null [%p, %p]") % s).str ());
      }

    ns_ptr = ns;
    ptr = s;
    trim_left (&ptr);
      
    for (; *ptr != '\0' && words_number < flag;)
      {
        if (*ptr == ' ' || *ptr == '\t')
          {
            ++words_number;
            *ns_ptr = ' ';
            ++ns_ptr;
            trim_left (&ptr);
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
  }

} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_DATA_STORAGE_STR_TOOLS_H_

