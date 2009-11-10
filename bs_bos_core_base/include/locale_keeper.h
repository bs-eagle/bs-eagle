/**
 * \file locale_keeper.h
 * \brief Helper class for set locale and reset locale on dtor
 * \author Miryanov Sergey
 * \date 14.04.2008
 */

#ifndef BS_LOCALE_KEEPER_H_
#define BS_LOCALE_KEEPER_H_

#include <locale.h>

namespace blue_sky
  {

  struct locale_keeper
    {
      char *locale_;
      int category_;

      locale_keeper (const char *new_name, int category_=LC_ALL)
          : locale_ (0)
          , category_ (category_)
      {
        locale_ = new char [strlen (setlocale (category_, 0)) + 1];
        memset (locale_, 0, strlen (setlocale (category_, 0)) + 1);
        memcpy (locale_, setlocale (category_, 0), strlen (setlocale (category_, 0)));
        setlocale (category_, new_name);
      }

      ~locale_keeper()
      {
        setlocale (category_, locale_);
        delete [] locale_;
      }
    };

}

#endif // #ifndef BS_LOCALE_KEEPER_H_
