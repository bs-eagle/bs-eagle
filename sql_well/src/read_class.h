#ifndef __READ_CLASS
#define __READ_CLASS

/*!
  \file read_class.h
  \brief Read functions for Keyword Input Language -- class #FRead, list of file included -- class #include_files
*/
#include <stdio.h>
#include <string.h>

#include "date_tools.h"
#include "localization.h"

#define FREAD_CONVERT_CASE 1
#define FREAD_DONT_CONVERT_CASE 0

// trim all right blanks
static inline char *
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
static inline char *
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

static inline int
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

#define SET_LINE_INFO(disp)                                                                     \
  if (disp < 0)                                                                                 \
    {                                                                                           \
      if (feof (this->fp))                                                                      \
        {                                                                                       \
          if (this->pop () != 0)                                                                \
            return -32;                                                          \
          goto beg;                                                                             \
        }                                                                                       \
      else                                                                                      \
        return -42;                                                        \
    }                                                                                           \
  IncrLineNumInFile (disp);


/*!
  \brief change pointer to the next not blanks
  \param start_ptr -- pointer to pointer to string
  \return if success  number of read doubles\n
*         if bad pointer '*start_ptr' -1
*/

static inline int
trim_left (char **start_ptr)
{
  if (!(*start_ptr))
    return -1;
  while ((**start_ptr == '\t' || **start_ptr == ' ' || **start_ptr == ','
          || **start_ptr == ':') && (**start_ptr != '\0'))
    ++(*start_ptr);
  return 0;
}

static inline int
get_phrase (char **next, char delim = ';')
{
  if (!(*next))
    return -1;
  char *start = *next;
  while ((**next) != delim && (**next) != '\0')
    {
      ++(*next);
    }
  if (**next != '\0')
    {
      **next = '\0';
      ++(*next);
    }
  trim_right_s (start);
  trim_left (next);
  return 0;
}

static inline int
get_phrase_str (char **next, char *buf, char delim = ';')
{
  trim_left (next);
  char *start = *next;
  if (get_phrase (next, delim))
    return -2;
  //printf ("NAMEMMMM: %s\n", start);
  if (*start == '*' || *start == '\0')
    return 0;
  strcpy (buf, start);
  return 0;
}
static inline int
get_phrase_int (char **next, int *i, char delim = ';')
{
  trim_left (next);
  char *start = *next;
  if (get_phrase (next, delim))
    return -2;
  if (*start == '*' || *start == '\0')
    return 0;
  if (sscanf (start, "%d", i) < 1)
    {
      fprintf (stderr, "Error: can not read int from %s\n", start);
      return -1;
    }
  return 0;
}
static inline int
get_phrase_double (char **next, double *d, char delim = ';')
{
  trim_left (next);
  char *start = *next;
  if (get_phrase (next, delim))
    {
      //printf ("get_phrase kkkkkk\n");
      return -2;
    }
  if (*start == '*' || *start == '\0')
    return 0;
  if (sscanf (start, "%lf", d) < 1)
    {
      fprintf (stderr, "Error: can not read double from %s\n", start);
      return -1;
    }
  return 0;
}

/*!
  \class R_FILE
  \brief class stores information about input keyword files
*/
class R_FILE
{
  public:
    R_FILE ()                     //! Constructor
      {
        name = 0;
        fp = 0;
        nstr = 0;
      }
    ~R_FILE ()                   //! Destructor
      {
        if (name)
          delete[]name;
        name = 0;
        if (fp)
          fclose (fp);
        fp = 0;
        nstr = 0;
      }
    FILE *get_fp () const         //! return file stream
      {
        return fp;
      }
    FILE *open (const char *Name) //! open --
      {
        if (!Name)
          return 0;
        if (fp)
          fclose (fp);
        fp = fopen (Name, "rt");
        if (!fp)
          return 0;
        if (set_name (Name))
          return 0;
        return fp;
      }
    void close ()
      {
        if (fp)
          fclose (fp);
        fp = 0;
        if (name)
          delete[]name;
        name = 0;
      }
    int set_name (const char *Name)       //! set new name
      {
        int l;
        if (!Name)
          return -1;
        if (name)
          delete[]name;
        l = (int) strlen (Name) + 1;
        name = new char[l];
        if (!name)
          return -2;
        strcpy (name, Name);
        return 0;                   // success
      }
    char *get_name () const       //! return pointer to the name
      {
        return name;
      }
    int nstr;                     //! position in the file
  private:
    char *name;                   //! file name
    FILE *fp;                     //! file descriptor
};


#ifdef UNIX
#define DIR_SYMBOL '/'
#else // !UNIX
#define DIR_SYMBOL '\\'
#endif //  UNIX

/*!
  \class include_file_node
  \ingroup KeywordLanguage
  \brief include_file_node -- list node for storing file name
*/
class include_file_node
{
public:
  //! default constructor
  include_file_node (const char *new_file = 0)
  {
    full_file_name = 0;
    file_name = 0;
    next = 0;
    if (new_file)
      set_file (new_file);
  }
  //! default destructor
  ~include_file_node ()
  {
    if (full_file_name)
      delete[]full_file_name;
    full_file_name = 0;
    file_name = 0;
    next = 0;
  }
  //! set new file name, return 0 if success, <0 if error occur
  int set_file (const char *new_file)
  {
    int len;

    if (full_file_name)
      delete[]full_file_name;
    full_file_name = 0;
    file_name = 0;
    if (!new_file)
      return -1;
    len = (int) strlen (new_file) + 1;

    // memory allocation
    full_file_name = new char[len];
    if (!full_file_name)
      return -2;
    // copy
    strcpy (full_file_name, new_file);
    // set up file name
    file_name = strrchr (full_file_name, DIR_SYMBOL);
    if (!file_name)
      file_name = full_file_name;
    else
      ++file_name;
    return 0;
  }
  //! return pointer to the file name with full path
  char *get_full_file_name () const
  {
    return full_file_name;
  }
  //! return pointer to the file name only
  char *get_file_name () const
  {
    return file_name;
  }
  //! return pointer to the next node, 0 if no next node
  include_file_node *get_next () const
  {
    return next;
  }
  //! set pointer to the next node
  void set_next (include_file_node * new_next)
  {
    next = new_next;
  }
private:
  char *full_file_name;         //!< file name with full puth
  char *file_name;              //!< pointer to place in full_file_name variable. Consist only from file name
  include_file_node *next;      //!< pointer to the next node, can be 0 if no next node
};

/*!
  \class include_files
  \ingroup KeywordLanguage
  \brief include_files -- list of included files
*/
class include_files
{
public:
  //! default constructor
  include_files ()
  {
    files_head = NULL;
    files_tail = NULL;
  }
  //! default destructor
   ~include_files ()
  {
    clear_list ();
  }
  //! get pointer to the first node
  include_file_node *get_first () const
  {
    return files_head;
  }
  //! clear files list
  void clear_list ()
  {
    include_file_node *node, *next;
    for (node = get_first (); node; node = next)
      {
        next = node->get_next ();
        delete node;
      }
    files_head = 0;
    files_tail = 0;
  }
  //! add new file to the list
  int add_file (const char *new_file)
  {
    include_file_node *node = 0;

    if (!new_file)
      return -1;
    node = new include_file_node (new_file);
    if (!node)
      return -2;
    if (!files_head)
      {
        files_head = node;
        //files_tail = node;
      }
    else
      files_tail->set_next (node);
    files_tail = node;
    return 0;
  }
  //! get full file path using given fule name
  char *get_full_path (const char *fname)
  {
    include_file_node *node = 0;
    for (node = get_first (); node; node = node->get_next ())
      {
        if (!strcmp (node->get_file_name (), fname))
          return node->get_full_file_name ();
      }
    return 0;
  }
private:
  //! pointer to the first node in list
  include_file_node * files_head;
  //! pointer to the last node in list
  include_file_node *files_tail;
};


/*!
  \class FRead
  \ingroup KeywordLanguage
  \brief FRead is a set of read function needed to read information from
*               input file, support include files
*/

class FRead
{
public:

  FRead ();
  FRead (const char *fname, const char *dir);
  ~FRead ();

  int set_incdir (const char *dir);
  char *get_incdir () const;
  void close_all ();
  int push (const char *fname, int main_file_p = 0);
  int pop ();
  int skip_d (char *start_ptr, char **end_ptr, int count);
  int skip_u (char *start_ptr, char **end_ptr, int count);
  int skip_s (char *start_ptr, char **end_ptr, int count);


  // read line from file
  int read_line (char *line, int max_len, int flg = FREAD_CONVERT_CASE);
  // read text block from file
  int read_text_block (char *line, int max_len);

  // For keyword KEY read array from file stream FP to buffer ARRAY.
  int read_int_array (const char *key,  //!< Name of calling keyword
                      int *array,       //!< string buffer
                      int len_array     //!< lenght of buffer
    );

  //! For keyword KEY read array from file stream FP to buffer ARRAY.
  int read_float_array (const char *key,       //!< Name of calling keyword
                         float *array, //!< string buffer
                         int len_array, //!< lenght of buffer
                         int first_flag = 0);
  //! For keyword KEY read array from file stream FP to buffer ARRAY.
  int read_double_array (const char *key,       //!< Name of calling keyword
                         double *array, //!< string buffer
                         int len_array, //!< lenght of buffer
                         int first_flag = 0);

  //! For keyword KEY read table from file stream FP to buffer DBUF.
  //! Each string contained NUM_COL double values.
  int read_double_table (const char *key,       //!< Input string of table
                         double *dbuf,  //!< Array for output
                         int max_len,   //!< Dimension for array DBUF
                         int num_col    //!< Number of column in string of table
    );

  // Turn on conversion and set conversion multiplier
  static int unwrap (char *s, const int flag = 0);
  static int scanf_s (char *start_ptr, char **end_ptr, char *dest, int no_default_p = 0, int *is_default = 0);
  static int scanf_text (char *start_ptr, char **end_ptr, char *dest);
  static int scanf_u (char *start_ptr, char **end_ptr, int *dest, int *is_default = 0);
  static int scanf_d (char *start_ptr, char **end_ptr, double *dest, int no_default_p = 0, int *is_default = 0);
  static int scanf_file_name (char *dest, char *source);
  char *get_prefix ();

  int eof ()
  {
    return feof (fp);
  }

  //! set pointer to the external include list
  void set_out_inc_list (include_files * ext_list)
  {
    out_inc_list = ext_list;
  }
private:

  void IncrLineNumInFile (const int val)        //!
  {
    if (top > 0)
      F[top - 1].nstr += val;
  }
  int convert_f (float *array, int len_array, int pos, char *buf,
                 const char *key);
  int convert_d (double *array, int len_array, int pos, char *buf,
                 const char *key);
  int convert_u (int *array, int len_array, int pos, char *buf,
                 const char *key);

  //---------------------------------------
  // VARIABLES
  // --------------------------------------
public:

  int inc_num;                  //!< number of included files
  include_files inc_list;
private:
  FILE *fp;
  char *incdir;
  include_files *out_inc_list;
  int top;
  int total;
  char prefix[2048];
  R_FILE *F;
};

#endif
