/**
 * @file property_base.cpp
 * @brief
 * @author Borschuk Oleg
 * @date 2008-03-29
 */
#include "property_base.h"

namespace blue_sky
  {
  //constructors
  property_base::property_base(bs_type_ctor_param /*param*/)
      : bs_refcounter(),objbase()
  {}

  property_base::property_base(const property_base& src)
      : bs_refcounter(),objbase(src)
  {
    *this = src;
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(property_base)
  BLUE_SKY_TYPE_STD_COPY(property_base)
  BLUE_SKY_TYPE_IMPL_SHORT(property_base, objbase, "BOS_Core property_base class")

  namespace allowed_types
    {
    BS_API_PLUGIN void cast_array_type(const int&) {}
    BS_API_PLUGIN void cast_array_type(const double&) {}
    BS_API_PLUGIN void cast_array_type(const bool&) {}
    BS_API_PLUGIN void cast_array_type(const std::string &) {}
  }
}//ns bs
