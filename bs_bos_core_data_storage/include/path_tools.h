/**
 * \file path_tools.h
 * \brief
 * \author Sergey Miryanov
 * \date 27.08.2009
 * */
#ifndef BS_BOS_CORE_DATA_STORAGE_PATH_TOOLS_H_
#define BS_BOS_CORE_DATA_STORAGE_PATH_TOOLS_H_

#define DIR_SYMBOL_UNIX '/'
#define DIR_SYMBOL_WIN  '\\'

#ifdef UNIX
  #define DIR_SYMBOL DIR_SYMBOL_UNIX
#else // !UNIX
  #define DIR_SYMBOL DIR_SYMBOL_WIN
#endif //  UNIX


namespace blue_sky 
  {
    namespace detail {

      inline const char *
      strrchr (const char *str)
      {
        const char *x = ::strrchr (str, DIR_SYMBOL);
        if (x)
          return x;

        x = ::strrchr (str, DIR_SYMBOL_UNIX);
        if (x)
          return x;

        return 0;
      }

    } // namespace detail
    namespace path 
      {
        inline std::string
        dirname (const std::string &path)
        {
          const char *x = detail::strrchr (path.c_str ());
          if (x)
            return std::string (path.c_str (), x + 1);

          return "";
        }

        inline bool
        is_absolute_path (const std::string &path)
        {
          BS_ASSERT (path.length ());

  #ifdef UNIX
          return path.c_str ()[0] == DIR_SYMBOL;
  #else
          return path.length () > 2 && (path.c_str ()[1] == ':' || path.c_str ()[0] == DIR_SYMBOL || path.c_str ()[0] == DIR_SYMBOL_UNIX);
  #endif
        }

        inline std::string
        join (std::string path, const std::string &suffix)
        {
          if (path.length ())
            return path + DIR_SYMBOL + suffix;
          else
            return suffix;
        }

      } // namespace path
  } // namespace blue_sky

#endif //#ifndef BS_BOS_CORE_DATA_STORAGE_PATH_TOOLS_H_

