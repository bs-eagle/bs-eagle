/** 
 * @file misc.h
 * @brief misc functions for common types
 * @author w63f
 * @version 0.1
 * @date 2012-12-19
 */
#ifndef COMMON_TYPES_MISC_H

#define COMMON_TYPES_MISC_H

#include <string>

//std::wstring stow(const std::string &s);
//std::string wtos(const std::wstring &s);

#define wtos(ret, s_)\
  {\
    std::wstring s = std::wstring (s_);\
    std::string d (s.length(), L' ');\
    std::copy (s.begin(), s.end(), d.begin());\
    ret = d;\
  }\

#define stow(ret, s_)\
  {\
    std::string s = std::string (s_);\
    std::wstring d (s.length(), L' ');\
    std::copy (s.begin(), s.end(), d.begin());\
    ret = d;\
  }\

#endif 

