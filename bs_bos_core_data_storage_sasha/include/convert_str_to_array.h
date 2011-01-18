/** 
 * \file convert_str_to_array.h
 * \brief convert string to array of items
 * \author Sergey Miryanov
 * \date 08.06.2009
 * */
#ifndef BS_BOS_CORE_DATA_STORAGE_CONVERT_STR_TO_ARRAY_H_
#define BS_BOS_CORE_DATA_STORAGE_CONVERT_STR_TO_ARRAY_H_

namespace blue_sky {
namespace data_reader {
namespace detail {

  template <typename T>
  struct convert_helper 
  {

  };

  template <>
  struct convert_helper <float>
  {
    static inline double 
    read (const char *start, char **end)
    {
      return strtod (start, end);
    }
    static inline size_t
    loop_condition (double t)
    {
      return (size_t) floor (t - 1 + 0.5);
    }
  };
  template <>
  struct convert_helper <double>
  {
    static inline double 
    read (const char *start, char **end)
    {
      return strtod (start, end);
    }
    static inline size_t
    loop_condition (double t)
    {
      return (size_t) floor (t - 1 + 0.5);
    }
  };

  template <>
  struct convert_helper <int>
  {
    static inline size_t 
    read (const char *start, char **end)
    {
      return strtol (start, end, 10);
    }
    static inline size_t 
    loop_condition (double t)
    {
    return t - 1;
    }
  };

  template <>
  struct convert_helper <unsigned char>
  {
    static inline size_t
    read (const char *start, char **end)
    {
      return strtol (start, end, 10);
    }
    static inline size_t 
    loop_condition (double t)
    {
      return t - 1;
    }
  };
  }

  /*!
    \brief   Recursive function\n
  *          For keyword KEY read array from
  *          file stream this->fp to buffer ARRAY.\n
  *          string format: 2*15.8 3*{12 2*5.6 2*{1.1 1.2}} is equals\n
  *          15.8 15.8 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2 12 5.6 5.6 1.1 1.2 1.1 1.2

    \param key    Name of calling keyword
    \param array  string buffer
    \param len_array  lenght of buffer
    \param pos number of doubles have been in array
    \param buf

    \return if success                                      number of read doubles\n
  *         if bad pointer 'array'                          -1\n
  *         if bad pointer 'buf'                            -2\n
  *         if cann't allocate memory                       -3\n
  *         if string format error                          -4
  */
  template <typename array_t, typename buffer_t>
  size_t 
  convert (array_t &array, 
    size_t offset,
    size_t pos, 
    buffer_t &buf,
    const std::string & key,
    const std::string &prefix)
  {
    typedef typename array_t::value_type item_t;
    buffer_t sbuf;

    size_t len_array = array.size();
    size_t cb, c, i, j, counter;
    size_t k;
    char *start_ptr, *end_ptr = 0;
    double t;
    

    if (array.empty ())               // check array pointer
      bs_throw_exception ("Array is empty");

    if (buf.empty ())
      bs_throw_exception ("Buf is null");

    if (pos >= len_array)         // check for input parameter
      {
        BOSWARN (section::read_data, level::low) << "Pos is out of array length" << bs_end;
        return YS_SUCCESS;
      }

    start_ptr = &buf[0];          // set start pointer to begin of buf
    counter = 0;                  // set up counter

    // main loop
    for (;;)
      {
        // check for garbage
        if (pos >= len_array)
          {
            if (*end_ptr != '\0')
              {
                BOSWARN (section::read_data, level::warning)
                << "Warning in " << prefix
                << ": trailing garbage " << end_ptr
                << " is ignored for keyword "
                << key << bs_end;
              }
            return counter;
          }
        trim_left (&start_ptr);

        t = detail::convert_helper <item_t>::read (start_ptr, &end_ptr); 
        trim_left (&end_ptr);

        if (*start_ptr == '\0')
          {
            return counter;
          }
        if (start_ptr == end_ptr) // if have not read return error -4
          {
            bs_throw_exception ("Nothing to read");
          }
        else if (*end_ptr == '*') // if next character is '*'
          {
            ++end_ptr;
            trim_left (&end_ptr);

            if (*end_ptr == '{')
              {
                ++end_ptr;
                cb = 1;
                k = 0;
                while (cb != 0 && *end_ptr != '\0')
                  {
                    if (*end_ptr == '{')
                      ++cb;
                    else if (*end_ptr == '}')
                      --cb;
                    sbuf[k] = *end_ptr;
                    ++k;
                    ++end_ptr;
                  }
                if (sbuf[k - 1] == '}' && k > 0)
                  {
                    --k;
                    sbuf[k] = '\0';
                  }
                else
                  {
                    bs_throw_exception ("return -40;");
                  }
                k = convert (array, offset, pos, sbuf, key, prefix);
                if (k <= 0)
                  {
                    return k;
                  }
                c = pos;
                pos += k;
                counter += k;
                for (i = 0; i < detail::convert_helper <item_t>::loop_condition (t); ++i)
                  {
                    for (j = 0; j < k; ++j, ++pos, ++counter)
                      {
                        if (pos < len_array)
                          array[pos + offset] = array[c + j + offset];
                        else
                          {
                            BOSWARN (section::read_data, level::warning)
                            << "Warning in " << prefix
                            << ": trailing garbage "
                            << end_ptr << "is ignored for keyword "
                            << key << bs_end;
                            return counter;
                          }
                      }
                  }
              }
            else
              {
                start_ptr = end_ptr;
                if (pos < len_array)
                  array[pos+offset] = static_cast <typename array_t::value_type> (detail::convert_helper <item_t>::read (start_ptr, &end_ptr));
                else
                  {
                    BOSWARN (section::read_data, level::warning)
                    << "Warning in " << prefix
                    << ": trailing garbage " << end_ptr
                    << " is ignored for keyword "
                    << key << bs_end;
                    return counter;
                  }
                c = pos;
                ++pos;
                ++counter;
                if (start_ptr == end_ptr) // if have not read return error -4
                  {
                    bs_throw_exception ("Nothing to read");
                  }
                else
                  {
                    for (i = 0; i < detail::convert_helper <item_t>::loop_condition (t); ++i, ++pos, ++counter)
                      {
                        if (pos < len_array)
                          array[pos+offset] = array[c+offset];
                        else
                          {
                            BOSWARN (section::read_data, level::warning)
                            << "Warning in " << prefix
                            << ": trailing garbage " << end_ptr
                            << " is ignored for keyword "
                            << key << bs_end;
                            return counter;
                          }
                      }
                  }
              }
          }
        else
          {
            if (pos < len_array)
              array[pos+offset] = t;
            else
              {
                BOSWARN (section::read_data, level::warning)
                << "Warning in " << prefix
                << ": trailing garbage " << end_ptr
                << " is ignored for keyword "
                << key << bs_end;
                return counter;
              }
            ++pos;
            ++counter;
          }
        start_ptr = end_ptr;
        if (*start_ptr == '\0')
          {
            return counter;
          }
      }
  }


} // namespace data_reader
} // namespace blue_sky



#endif  // #ifndef BS_BOS_CORE_DATA_STORAGE_CONVERT_STR_TO_ARRAY_H_
