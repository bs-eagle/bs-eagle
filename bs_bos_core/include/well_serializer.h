/**
 *       \file  well_serializer.h
 *      \brief  Class for well data serialization
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete
 * */
#ifndef BS_WELL_SERIALIZER_H_
#define BS_WELL_SERIALIZER_H_

#include "data_storage_interface.h"

namespace blue_sky
  {

  namespace wells
    {

    class BS_API_PLUGIN well_serializer : public data_serializer
      {
      public:
        virtual ~well_serializer () {}

        virtual void save (const sp_storage_t &storage, const sp_obj &obj) const;

        BLUE_SKY_TYPE_DECL (well_serializer);
      };

    //////////////////////////////////////////////////////////////////////////
    bool
    well_serializer_register_type (const plugin_descriptor &pd);


  } // namespace wells
} // namespace blue_sky


#endif  // #ifndef BS_WELL_SERIALIZER_H_
