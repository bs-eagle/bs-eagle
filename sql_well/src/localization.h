/*!
  \file localization.h
  \brief defines for localization
*/
#ifndef __LOCALIZATION_H
#define __LOCALIZATION_H

//! English language
#define LC_ENG 0
//! Russian language
#define LC_RUS 5

//! Russian WIN coding  
#define LC_RUS_WIN  1
//! Russian ALT coding
#define LC_RUS_ALT  2
//! Russian KOI8 coding                
#define LC_RUS_KOI  3

int locale_set_lang (int);      // Set output language
const char *locale_get_text (const char *);     // Get text translate and encoding
int binary_search_loc (const char *, const int);
int compare_keyword (const char *, const char *);
char locale_toupper (unsigned char);    // toupper 

char *locale_ucase (char *s);   // convert string <s> to upper case

//#define GET_TEXT(X)     (X)
//! macros of locale_get_text (s)
//#define GET_TEXT(X)     (locale_get_text (X))
#define GET_TEXT(X)     (X)

#ifdef _UFASOLVER_FACE_
int face_locale_set_lang (int);      // Set output language for face
const char *face_locale_get_text (const char *);     // Get text translate and encoding
#define FACE_GET_TEXT(X)   (face_locale_get_text (X))
#define SYSINFO_GET_TEXT(FL,X) ((FL) ? FACE_GET_TEXT(X) : GET_TEXT(X)) 
#else
#define SYSINFO_GET_TEXT(FL,X) (GET_TEXT(X))
#endif // _UFASOLVER_FACE_

//! macros of compare with keyword (compare_keyword (X, Y))
#define COMPARE_KEYWORD(X,Y)  (compare_keyword (X, Y))
//#define COMPARE_KEYWORD(X,Y)  (strcmp (X,Y))

//! macros of to upper function
#define TOUPPER(X)      (locale_toupper (X))

#endif // __LOCALIZATION_H
