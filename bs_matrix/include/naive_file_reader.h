/**
 * \file naive_file_reader.h
 * \brief read plain file breaking on sections
 * \author Sergey Miryanov
 * \date 11.06.2008
 * */
#ifndef BS_NAIVE_FILE_READER_H_
#define BS_NAIVE_FILE_READER_H_

#include <cstdio>

#include "shared_vector.h"
#include "locale_keeper.h"
#include "bs_assert.h"

namespace blue_sky
  {

  struct naive_file_reader
    {
      naive_file_reader (const char *filename)
          : lkeeper ("C", LC_ALL)
      {
        file = fopen (filename, "r");
        BS_ASSERT (file) (filename);

        memset (line, 0, sizeof (line));
      }

      ~naive_file_reader ()
      {
        if (file)
          fclose (file);
      }

      naive_file_reader &
      locate_section (const char *section)
      {
        if (!file)
          return *this;

        while (true)
          {
            if (strstr (line, section))
              break;

            char *read = fgets (line, sizeof (line), file);
            BS_ASSERT (read) (ftell (file));

            if (feof (file))
              break;
          }

        return *this;
      }

      template <class array_t> naive_file_reader &
      read_list (array_t &array)
      {
        if (!file)
          return *this;

        while (!feof (file))
          {
            typename array_t::value_type value = 0;
            if (!fgets (line, sizeof (line), file))
              break;

            int n = 0;
            sscanf (line, get_format (value), &value, &n);
            if (n == 0)
              break;

            array.push_back (value);
          }

        return *this;
      }

      template <class item_t> naive_file_reader &
      read_item (item_t &item)
      {
        if (!file || feof (file))
          return *this;

        int n = 0;
        if (fscanf (file, get_format (item), &item, &n) != 2)
          {
            //TODO: raise error
          }
        BS_ASSERT (n);
        return *this;
      }

      naive_file_reader &
      skip_line ()
      {
        if (!file || feof (file))
          return *this;

        static char tmp[1000] = {0};
        if (!fgets (tmp, sizeof (tmp), file))
          {
            //TODO: raise error
          }

        return *this;
      }

private:

      static const char *
      get_format (int)
      {
        return "%d%n";
      }
      static const char *
      get_format (float)
      {
        return "%f%n";
      }
      static const char *
      get_format (double)
      {
        return "%lf%n";
      }

private:

      locale_keeper   lkeeper;
      FILE            *file;
      char            line[4096];
    };

} // namespace blue_sky


#endif  // #ifndef BS_NAIVE_FILE_READER_H_
