/**
 * \file scale_array_holder.h
 * \brief
 * \author Sergey Miryanov
 * \date 02.12.2008
 * */
#ifndef BS_SCAL_SCALE_ARRAY_HOLDER_H_
#define BS_SCAL_SCALE_ARRAY_HOLDER_H_

#include "scale_arrays_placement_strategies.h"
#include "shared_vector.h"
#include "vartype_table_iface.h"
#include "bs_serialize_decl.h"

namespace blue_sky
  {

  struct value_accessor
    {
      value_accessor (const t_float *array_, int count_, double value_)
          : count_ (count_)
          , array_ (array_)
          , value_ (value_)
      {

      }

      double operator [] (int index) const
      {
        return count_ ? array_[index] : value_;
      }

private:

      int             count_;
      const t_float * array_;
      double          value_;

    };

  class BS_API_PLUGIN scale_array_holder : public scale_array_holder_iface
    {
    public:
      typedef vartype_table_iface <t_float>     table_t;
      typedef table_t::vector_t                 vector_t;
      typedef smart_ptr<table_t, true>          sp_vartype_table; 

    public:

      value_accessor 
      get (scale_array_name array, t_float default_value) const
      {
        return value_accessor (data_pool->get_col_ptr (array), data_pool->get_n_rows (array), default_value);
      }

      void
      set (scale_array_name array, std::wstring const &name, spv_float const &data)
      {
        if (data->size ())
          data_pool->add_col_vector (array, name, data);
      }

      // FIXME:
      t_double
      get_su (t_long cell, t_double value) const
      {
        return get (su, value)[cell];
      }

      // FIXME:
      t_double
      get_sl (t_long cell, t_double value) const
      {
        return get (sl, value)[cell];
      }

      void remove (scale_array_name array)
      {
        data_pool->remove_col_vector (array);
      }
      
      int 
      is_prop_valid (const t_long col) const
        {
          return (data_pool->get_n_rows (col) == 0) ? 0 : 1;
        }
      
    private:

      sp_vartype_table data_pool;
    public:

      BLUE_SKY_TYPE_DECL (scale_array_holder);

      friend class blue_sky::bs_serialize;
    };


} // namespace blue_sky



#endif  // #ifndef BS_SCAL_SCALE_ARRAY_HOLDER_H_

