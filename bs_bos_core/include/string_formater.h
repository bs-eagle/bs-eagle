/**
 *       \file  string_formater.h
 *      \brief  Helper to format string (via ellipses),
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  31.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_TOOLS_STRING_FORMATER_H_
#define BS_TOOLS_STRING_FORMATER_H_

namespace blue_sky
  {
  namespace tools
    {

      /**
       * \class string_formater
       * \brief Formats strings...
       * */
    struct string_formater
      {
        /**
         * \brief  format string should contains only one %d specifier
         * \param  format
         * \param  i
         * */
        string_formater (const std::string &format, size_t i)
        {
          memset (str, 0, sizeof (str));
          sprintf (str, format.c_str (), i);
        }

        /**
         * \brief  format string should contains only two %f and 
         *         one %d specifiers, in such order
         * \param  format
         * \param  dt
         * \param  ct
         * \param  i
         * */
        string_formater (const char *format, double dt, double ct, int i)
        {
          memset (str, 0, sizeof (str));
          sprintf (str, format, dt, ct, i);
        }

        /**
         * \brief  format string should contains only one %s specifier
         * \param  format
         * \param  s
         * */
        string_formater (const char *format, const std::string &s)
        {
          memset (str, 0, sizeof (str));
          sprintf (str, format, s.c_str ());
        }

        //! Formated string
        char str[255];
      };



  } // namespace tools
} // namespace blue_sky



#endif  // #ifndef BS_TOOLS_STRING_FORMATER_H_

