/*! \file named_pbase_access.cpp
		\brief named_pbase class methods implementations
		\author Nikonov Max
*/
#include "bs_bos_core_base_stdafx.h"

#include "named_pbase_access.h"
#include "bs_kernel.h"

using namespace std;

namespace blue_sky
  {
//! constructor
  named_pbase::named_pbase (bs_type_ctor_param /*param*/)
      : bs_refcounter(), property_base()
  {
    set_default_values ();
    setup_assoc_names();
  }

  named_pbase::named_pbase(const named_pbase& prop)
      : bs_refcounter(), property_base(prop)
  {
    if (&prop != this)
      *this = prop;

    set_default_values ();
  }

  //! destructor
  named_pbase::~named_pbase ()
  {
    set_default_values ();
  }

  void named_pbase::setup_assoc_names()
  {
    for (int i = 0; i < (int)params_names.size(); ++i)
      {
        const std::string &tmp = params_names[i].name; 
        if (tmp.length ())
          assoc_names[tmp] = (named_pbase_idxs)i;
      }
  }

  bool
  named_pbase::set_value (int index, const std::string &str)
  {
    if (index >= (int)params_names.size())
      {
        BOSERR (section::read_data, level::error) << "Index out of range (" << index << ", " << params_names.size () << ")" << bs_end;
        return false;
      }

    if (params_names[index].type == PT_INT)
      {
        try
          {
            this->set_int ((idx_type)index, boost::lexical_cast<int> (str));
          }
        catch (const boost::bad_lexical_cast &)
          {
            bs_throw_exception (boost::format ("Wrong value type (index: %d, conversion to int failed)") % index);
          }
      }
    else if (params_names[index].type == PT_FLOAT)
      {
        try
          {
            this->set_float ((idx_type)index, boost::lexical_cast<double> (str));
          }
        catch (const boost::bad_lexical_cast &)
          {
            bs_throw_exception (boost::format ("Wrong value type (index: %d, conversion to double failed)") % index);
          }
      }
    else if (params_names[index].type == PT_STR)
      {
        this->set_param ((idx_type)index, str);
      }
    else
      {
        bs_throw_exception ("I have no idea how it happens");
      }

    return true;
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(named_pbase)
  BLUE_SKY_TYPE_STD_COPY(named_pbase)
  BLUE_SKY_TYPE_IMPL_SHORT(named_pbase, objbase, "BOS_Core named_pbase class")
}
