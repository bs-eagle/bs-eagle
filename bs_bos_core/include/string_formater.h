/**
 * \file string_formater.h
 * \brief helper to format string (via ellipses) (for debug purpose only)
 * \author Sergey Miryanov
 * \date 31.10.2008
 * */
#ifndef BS_TOOLS_STRING_FORMATER_H_
#define BS_TOOLS_STRING_FORMATER_H_

namespace blue_sky
  {
  namespace tools
    {

    struct string_formater
      {
        string_formater (const std::string &format, size_t i)
        {
          memset (str, 0, sizeof (str));
          sprintf (str, format.c_str (), i);
        }

        string_formater (const char *format, double dt, double ct, int i)
        {
          memset (str, 0, sizeof (str));
          sprintf (str, format, dt, ct, i);
        }

        string_formater (const char *format, const std::string &s)
        {
          memset (str, 0, sizeof (str));
          sprintf (str, format, s.c_str ());
        }

        char str[255];
      };



  } // namespace tools
} // namespace blue_sky



#endif  // #ifndef BS_TOOLS_STRING_FORMATER_H_

