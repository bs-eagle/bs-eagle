#ifndef ASCII_READER_XCUW9UDS

#define ASCII_READER_XCUW9UDS

/*!
  \file ascii_reader.h
  \brief Read functions for Keyword Input Language -- class #FRead, list of file included -- class #include_files
*/
#include <stdio.h>
#include <string.h>
#include "conf.h"

//#include "date_tools.h"

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
  t_long str_pos, length;

  if (s == NULL)
    return NULL;
  if ((length = (t_long)strlen (s)) == 0)
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
get_next_not_empty_line (char *string, t_long MaxLen, FILE * fp, t_long *len)
{
  char *ptr_s = string;
  int ch;
  t_long disp = 0;

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
        t_long l;
        if (!Name)
          return -1;
        if (name)
          delete[]name;
        l = (t_long) strlen (Name) + 1;
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
    t_long nstr;                     //! position in the file
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
    t_long len;

    if (full_file_name)
      delete[]full_file_name;
    full_file_name = 0;
    file_name = 0;
    if (!new_file)
      return -1;
    len = (t_long) strlen (new_file) + 1;

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
  int skip_d (char *start_ptr, char **end_ptr, const t_long count) const;
  int skip_u (char *start_ptr, char **end_ptr, const t_long count) const;
  int skip_s (char *start_ptr, char **end_ptr, const t_long count) const;


  // read line from file
  t_long read_line (char *line, const t_long max_len, const int flg = FREAD_CONVERT_CASE);
  // read text block from file
  t_long read_text_block (char *line, t_long const max_len);

  // Turn on conversion and set conversion multiplier
  int unwrap (char *s, const int flag = 0) const;
  int scanf_s (char *start_ptr, char **end_ptr, char *dest, int no_default_p = 0, int *is_default = 0) const;
  int scanf_text (char *start_ptr, char **end_ptr, char *dest) const;
  int scanf_u (char *start_ptr, char **end_ptr, t_long *dest, int *is_default = 0) const;
  int scanf_d (char *start_ptr, char **end_ptr, t_double *dest, int no_default_p = 0, int *is_default = 0) const;
  int scanf_file_name (char *dest, char *source) const;
  
  t_long convert_f (t_float *array, const t_long len_array, t_long pos, char *buf,
                 const char *key);
  //t_long convert_d (t_double *array, const long len_array, long pos, char *buf,
  //               const char *key);
  t_long convert_u (t_int *array, const t_long len_array, t_long pos, char *buf,
                 const char *key);

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

  void IncrLineNumInFile (const t_long val)        //!
  {
    if (top > 0)
      F[top - 1].nstr += val;
  }

  //---------------------------------------
  // VARIABLES
  // --------------------------------------
public:

  t_long inc_num;                  //!< number of included files
  include_files inc_list;
private:
  FILE *fp;
  char *incdir;
  include_files *out_inc_list;
  t_long top;
  t_long total;
  char prefix[2048];
  R_FILE *F;
};

#endif /* end of include guard: ASCII_READER_XCUW9UDS */
