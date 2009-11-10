/**
 * \file well_serializer.h
 * \brief class for serialize well data
 * \author Sergey Miryanov
 * \date 21.07.2008
 * */
#ifndef BS_WELL_SERIALIZER_H_
#define BS_WELL_SERIALIZER_H_

#include "data_storage_interface.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class well;

  namespace wells
    {

    template <typename strategy_t>
    class BS_API_PLUGIN well_serializer : public data_serializer
      {
      public:

        typedef well <strategy_t>         well_t;
        typedef smart_ptr <well_t, true>  sp_well_t;

        virtual ~well_serializer () {}

        virtual void save (const sp_storage_t &storage, const sp_obj &obj) const;

        BLUE_SKY_TYPE_DECL_T (well_serializer);
      };

    //////////////////////////////////////////////////////////////////////////
    bool
    well_serializer_register_type (const plugin_descriptor &pd);


  }	// namespace wells
}	// namespace blue_sky


#endif	// #ifndef BS_WELL_SERIALIZER_H_
