/**
 *       \file  well_limit_operation.h
 *      \brief  Class well_limit_operation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  14.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete
 * */
#ifndef BS_WELL_LIMIT_OPERATION_H_
#define BS_WELL_LIMIT_OPERATION_H_

namespace blue_sky
  {
  namespace wells
    {

    enum limit_operation_type
    {
      operation_none,
      operation_con,
      operation_well,
      operation_total,
    };

    limit_operation_type
    limit_operation_cast (const std::string &str);

    class BS_API_PLUGIN well_limit_operation : public objbase
      {
      public:

        typedef double item_t;

      public:

        void set_min_oil_rate (item_t min_oil_rate);
        void set_max_water_cut (item_t max_water_cut);
        void set_min_water_cut (item_t min_water_cut);

        void set_value (int value_type, item_t value);

        bool fi_check_limits () const;

        BLUE_SKY_TYPE_DECL (well_limit_operation);

      private:

        //friend class well_limit_operation;

        item_t      min_oil_rate_;
        item_t      max_water_cut_;
        item_t      min_water_cut_;

      };

    class BS_API_PLUGIN well_limit_operation_factory : public objbase
      {
      public:

        typedef smart_ptr <well_limit_operation>    sp_well_limit_operation_t;

      public:

        virtual ~well_limit_operation_factory () {}

        virtual sp_well_limit_operation_t create_limit (limit_operation_type type) const;

        BLUE_SKY_TYPE_DECL (well_limit_operation_factory);
      };

    bool
    well_limit_operation_register_type (const blue_sky::plugin_descriptor &pd);

    bool
    well_limit_operation_factory_register_type (const blue_sky::plugin_descriptor &pd);


  } // namespace wells
} // namespace blue_sky


#endif  // #ifndef BS_WELL_LIMIT_OPERATION_H_
