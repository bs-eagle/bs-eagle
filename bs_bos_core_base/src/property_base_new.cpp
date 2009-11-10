/**
 * @file property_base.cpp
 * @brief
 * @author Borschuk Oleg
 * @date 2008-03-29
 */
#include "stdafx.h"

#include "property_base.h"
#include "bs_kernel.h"
#include "bs_link.h"

namespace blue_sky
  {
  //constructors
  property_base::property_base(bs_type_ctor_param param)
      : bs_refcounter(),objbase()
  {}

  property_base::property_base(const property_base& src)
      : bs_refcounter(),objbase(src)
  {
    *this = src;
  }

  void property_base::resize_int (int size)
  {
    resize(size,ivals,is_inited_i);
  }

  void property_base::resize_float (int size)
  {
    resize(size,fvals,is_inited_f);
  }

  void property_base::resize_bool (int size)
  {
    resize(size,bvals,is_inited_b);
  }

  void property_base::resize_string (int size)
  {
    resize(size,svals,is_inited_s);
  }

  void property_base::clear ()
  {
    ivals.clear();
    fvals.clear();
    bvals.clear();
    svals.clear();

    is_inited_i.clear ();
    is_inited_f.clear ();
    is_inited_b.clear ();
    is_inited_i.clear ();
  }

  property_base& property_base::operator= (const property_base &rhs)
  {
    ivals = rhs.ivals;
    fvals = rhs.fvals;
    bvals = rhs.bvals;
    svals = rhs.svals;

    is_inited_i = rhs.is_inited_i;
    is_inited_f = rhs.is_inited_f;
    is_inited_b = rhs.is_inited_b;
    is_inited_s = rhs.is_inited_s;

    return *this;
  }

  property_base& property_base::operator+= (const property_base &rhs)
  {
    BS_ASSERT ("In operator += " && (ivals.size() == rhs.ivals.size())
               && (fvals.size() == rhs.fvals.size())
               && (bvals.size() == rhs.bvals.size())
               && (svals.size() == rhs.svals.size())
               && "can be used for arrays of equil size only.");
    for (int i = 0; i < ivals.size(); ++i)
      {
        if (rhs.check_value ((property_base_idxs_int)i))
          {
            ivals[i] = rhs.ivals[i];
            is_inited_i.set (i);
          }
      }

    for (int i = 0; i < fvals.size(); ++i)
      {
        if (rhs.check_value ((property_base_idxs_float)i))
          {
            fvals[i] = rhs.fvals[i];
            is_inited_f.set (i);
          }
      }

    for (int i = 0; i < bvals.size(); ++i)
      {
        if (rhs.check_value ((property_base_idxs_bool)i))
          {
            bvals[i] = rhs.bvals[i];
            is_inited_b.set (i);
          }
      }

    for (int i = 0; i < svals.size(); ++i)
      {
        if (rhs.check_value ((property_base_idxs_string)i))
          {
            svals[i] = rhs.svals[i];
            is_inited_s.set (i);
          }
      }

    return *this;
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(property_base)
  BLUE_SKY_TYPE_STD_COPY(property_base)
  BLUE_SKY_TYPE_IMPL_SHORT(property_base, objbase, "BOS_Core property_base class")

  namespace allowed_types
    {
    BS_API_PLUGIN void cast_array_type(const int&) {}
    BS_API_PLUGIN void cast_array_type(const float&) {}
    BS_API_PLUGIN void cast_array_type(const bool&) {}
    BS_API_PLUGIN void cast_array_type(const std::string&) {}
  }
}//ns bs
