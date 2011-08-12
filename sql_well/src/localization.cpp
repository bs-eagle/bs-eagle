/*!
  \file localization.cpp
  \brief localization functions: set language tables and get translate
*/
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "localization.h"

#include "table_upper.h"

static int locale_lang = LC_ENG;        //!< Current locale lang by default 
//static const char **current_locale;     //!< pointer to current locale table
static unsigned char *upper_table = 0;  //!< pointer to upper case encoding table 

#ifdef _UFASOLVER_FACE_
static int face_locale_lang = LC_ENG;        //!< Current locale face lang by default 
static const char **face_current_locale;     //!< pointer to current face locale table
static unsigned char *face_upper_table = 0;  //!< pointer to face upper case encoding table 
#endif // _UFASOLVER_FACE_

#if 0
/*!
  \brief set output language\n
         set pointer to language table\n
         set pointer to UPPER TABLE
  \param id  - language id
  \return if success                     0\n
          if not supported language     -1 
*/
int
locale_set_lang (int id)
{
  switch (id)
    {
    case LC_ENG:
      locale_lang = id;
      current_locale = lc_messages_eng;
      upper_table = 0;
      break;
    case LC_RUS_WIN:
      locale_lang = id;
      current_locale = lc_messages_rus_win;
      upper_table = table_toupper_win;
      break;
    case LC_RUS_ALT:
      locale_lang = id;
      current_locale = lc_messages_rus_alt;
      upper_table = table_toupper_alt;
      break;
    case LC_RUS_KOI:
      locale_lang = id;
      current_locale = lc_messages_rus_koi;
      upper_table = table_toupper_koi;
      break;
    default:
      return -1;
    }
  return 0;
}

/*!
  \brief binary search of text
  \param word -- pointer to searching text
  \param n    -- lenght of table
  \return if success   -- index of word in table\n
          else --  -1
*/
int
binary_search_loc (const char *word, const int n)
{
  int mid;
  int low = 0;
  int high = n - 1;
  while (low != high)
    {
      mid = (low + high) / 2;
      if (strcmp (word, lc_messages_eng[mid]) > 0)
        low = mid + 1;
      else
        high = mid;
    }
  if (strcmp (word, lc_messages_eng[low]) != 0)
    return -1;
  else
    return low;
}

/*!
  \brief find translate in table
  \param s text to translate
  \return translated text 
*/
const char *
locale_get_text (const char *s)
{
  if (locale_lang != LC_ENG)    // current_locale != lc_messages_eng)
    {
      int index;
      index = binary_search_loc (s, lc_messages_len);
      if (index == -1)
        {
          return s;
        }
      return current_locale[index];
    }
  else
    return s;
}
#endif //0
/*!
  \brief compare key with one of keywords 
  \param key -- word
  \param keyword  -- keyword
  \return result of compare 
*/
int
compare_keyword (const char *key, const char *keyword)
{
  if (locale_lang == LC_ENG)
    return (strcmp (key, keyword));
  if (strcmp (key, keyword) == 0)
    return 0;
  return (strcmp (key, GET_TEXT (keyword)));

}

/*!
  \brief toupper for russian letters
  \param c letter
  \return toupper(c) 
*/
inline char
locale_toupper (unsigned char c)
{
  if (c < 128 || !upper_table)
    return toupper (c);
  return upper_table[c - 128];
}

/*!
  \brief toupper for russian letters
  \param s string
  \return string pointer
*/
char *
locale_ucase (char *s)
{
  if (!s)
    return s;
  char *ptr = s;
  for (; *ptr != '\0'; ++ptr)
    *ptr = locale_toupper (*ptr);

  return s;
}

#ifdef _UFASOLVER_FACE_
/*!
  \brief set output language for face\n
         set pointer to language table\n
         set pointer to UPPER TABLE
  \param id  - language id
  \return if success                     0\n
          if not supported language     -1 
*/
int
face_locale_set_lang (int id)
{
  switch (id)
    {
    case LC_ENG:
      face_locale_lang = id;
      face_current_locale = lc_messages_eng;
      face_upper_table = 0;
      break;
    case LC_RUS_WIN:
      face_locale_lang = id;
      face_current_locale = lc_messages_rus_win;
      face_upper_table = table_toupper_win;
      break;
    case LC_RUS_ALT:
      face_locale_lang = id;
      face_current_locale = lc_messages_rus_alt;
      face_upper_table = table_toupper_alt;
      break;
    case LC_RUS_KOI:
      face_locale_lang = id;
      face_current_locale = lc_messages_rus_koi;
      face_upper_table = table_toupper_koi;
      break;
    default:
      return -1;
    }
  return 0;
}

/*!
  \brief find translate in table
  \param s text to translate
  \return translated text 
*/
const char *
face_locale_get_text (const char *s)
{
  if (face_locale_lang != LC_ENG)    // current_locale != lc_messages_eng)
    {
      int index;
      index = binary_search_loc (s, lc_messages_len);
      if (index == -1)
        {
          return s;
        }
      return face_current_locale[index];
    }
  else
    return s;
}
#endif // _UFASOLVER_FACE_

