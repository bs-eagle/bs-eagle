/**
 * \file path_tools.h
 * \brief
 * \author Sergey Miryanov
 * \date 27.08.2009
 * */
#ifndef BS_BOS_CORE_DATA_STORAGE_PATH_TOOLS_H_
#define BS_BOS_CORE_DATA_STORAGE_PATH_TOOLS_H_

namespace blue_sky 
  {
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

