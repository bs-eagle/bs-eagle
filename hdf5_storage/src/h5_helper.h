/** 
 * @file h5_helper.h
 * @brief helper functions for #h5_pool
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-12
 */
#ifndef H5_HELPER_WLU6BMTS

#define H5_HELPER_WLU6BMTS

#include <sys/stat.h>
#include <string>

bool file_exists (const std::string &fname) 
{
  struct stat file_info;
  bool ret;
  int int_stat;

  // Attempt to get the file attributes
  int_stat = stat (fname.c_str (), &file_info);
  if (int_stat == 0) 
    {
      // We were able to get the file attributes
      // so the file obviously exists.
      ret = true;
    } 
  else 
    {
      // We were not able to get the file attributes.
      // This may mean that we don't have permission to
      // access the folder which contains this file. If you
      // need to do that level of checking, lookup the
      // return values of stat which will give you
      // more details on why stat failed.
      ret = false;
    }
  
  return (ret);
}

#endif /* end of include guard: H5_HELPER_WLU6BMTS */
