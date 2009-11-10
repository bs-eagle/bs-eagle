/**
 * \file well_reporter.cpp
 * \brief prints well data
 * \date 29.09.2009
 * \author Sergey Miryanov
 * */
#include "stdafx.h"
#include "well_reporter.h"
#include "calc_well.h"
#include "facility_manager.h"
#include "well_connection.h"

#include <iostream>

#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/preprocessor/tuple.hpp>

#define LIQUID_RATE_MULT            data->output_units_converter.liquid_rate_mult ()
#define GAS_RATE_MULT               data->output_units_converter.gas_rate_mult ()
#define GAS_OIL_RATIO_MULT          data->output_units_converter.gas_liquid_mult ()
#define WATER_GAS_RATIO_MULT        1.0 / data->output_units_converter.gas_liquid_mult ()
#define PRESSURE_MULT               data->output_units_converter.pressure_mult ()
#define GAS_LIQUID_RATE_MULT        data->output_units_converter.gas_liquid_rate_mult ()
#define GRAVITY                     data->input_units_converter.output_constants.gravity_constant

#define RATE_LIQUID_SURFACE_NOTE    convert_units::rate_liquid_surface_note (data->output_units_converter.get_output_units ()) 
#define RATE_GAS_SURFACE_NOTE       convert_units::rate_gas_surface_note (data->output_units_converter.get_output_units ()) 
#define RATE_LIQUID_RESERVOIR_NOTE  convert_units::rate_liquid_reservoir_note (data->output_units_converter.get_output_units ()) 
#define GAS_LIQUID_NOTE             convert_units::gas_liquid_note (data->output_units_converter.get_output_units ())
#define LIQUID_GAS_NOTE             convert_units::liquid_gas_note (data->output_units_converter.get_output_units ())
#define PRESSURE_NOTE               convert_units::pressure_note (data->output_units_converter.get_output_units ())
#define VOL_LIQUID_SURFACE_NOTE     convert_units::volume_liquid_surface_note (data->output_units_converter.get_output_units ())
#define VOL_GAS_SURFACE_NOTE        convert_units::volume_gas_surface_note (data->output_units_converter.get_output_units ())
#define VOL_LIQUID_RES_NOTE         convert_units::volume_liquid_reservoir_note (data->output_units_converter.get_output_units ())

#define DECL_COLUMN_HOLDER_I(r, _, i, elem)                                           \
  column_t elem;

#define DECL_COLUMN_HOLDER_PTR_I(r, _, i, elem)                                       \
  columns_[i] = &elem;

#define DATA_LEN_I(r, _, idx, elem)                                                   \
  + 1 + columns.BOOST_PP_TUPLE_ELEM (3, 0, elem).width

#define DATA_I(r, func, i, elem)                                                      \
  sprintf (buffer + idx,                                                              \
    BOOST_PP_TUPLE_ELEM (3, 1, elem).format_string (                                  \
      columns.BOOST_PP_TUPLE_ELEM (3, 0, elem)),                                      \
    BOOST_PP_TUPLE_ELEM (3, 1, elem).data (                                           \
      columns.BOOST_PP_TUPLE_ELEM (3, 0, elem),                                       \
      ds,                                                                             \
      columns.BOOST_PP_TUPLE_ELEM (3, 0, elem).mult,                                  \
      BOOST_PP_TUPLE_ELEM (3, 2, elem)));                                             \
  idx += 1 + columns.BOOST_PP_TUPLE_ELEM (3, 0, elem).width;

#define DECL_COLUMN_HOLDER(holder_name, len, tuple)                                   \
  struct holder_name                                                                  \
  {                                                                                   \
    BOOST_PP_SEQ_FOR_EACH_I (DECL_COLUMN_HOLDER_I, _,                                 \
      BOOST_PP_TUPLE_TO_SEQ (len, tuple))                                             \
                                                                                      \
    column_t        *columns_[len];                                                   \
    size_t          columns_count_;                                                   \
    mutable char    *line_buffer;                                                     \
    mutable char    *unit_buffer;                                                     \
    mutable char    *header_buffer[column_t::line_count];                             \
    mutable size_t  header_buffer_count;                                              \
    holder_name ()                                                                    \
    : columns_count_ (len)                                                            \
    , line_buffer (0)                                                                 \
    , unit_buffer (0)                                                                 \
    , header_buffer_count (0)                                                         \
    {                                                                                 \
      BOOST_PP_SEQ_FOR_EACH_I (DECL_COLUMN_HOLDER_PTR_I, _,                           \
        BOOST_PP_TUPLE_TO_SEQ (len, tuple))                                           \
      detail::init_header_buffer (header_buffer);                                     \
    }                                                                                 \
    ~holder_name ()                                                                   \
    {                                                                                 \
      delete [] line_buffer;                                                          \
      delete [] unit_buffer;                                                          \
      detail::free_header_buffer (header_buffer);                                     \
    }                                                                                 \
  };

#define DATA_LINE(columns_type, name, data_type, type, seq)                           \
  struct BOOST_PP_CAT (name, _class)                                                  \
  {                                                                                   \
    mutable char *buffer;                                                             \
    BOOST_PP_CAT (name, _class) ()                                                    \
    : buffer (0)                                                                      \
    {                                                                                 \
    }                                                                                 \
    ~BOOST_PP_CAT (name, _class) ()                                                   \
    {                                                                                 \
      delete[] buffer;                                                                \
    }                                                                                 \
    void                                                                              \
    print (const columns_type &columns, const data_type &data, const type &ds,        \
      bool do_print_line = true) const                                                \
    {                                                                                 \
      size_t len = 1                                                                  \
        BOOST_PP_SEQ_FOR_EACH_I (DATA_LEN_I, _, seq);                                 \
      if (!buffer)                                                                    \
        {                                                                             \
          buffer = new char [len + 1];                                                \
        }                                                                             \
      memset (buffer, 0, len + 1);                                                    \
      size_t idx = 0;                                                                 \
      BOOST_PP_SEQ_FOR_EACH_I (DATA_I, _, seq)                                        \
      buffer[idx] = '|';                                                              \
      BOSOUT (section::print_reports, level::highest) << buffer << bs_end;            \
      if (do_print_line)                                                              \
        print_line (columns);                                                         \
    }                                                                                 \
  };                                                                                  \
  BOOST_PP_CAT (name, _class) name;

namespace blue_sky {

  namespace detail {
    template <size_t line_count>
    void
    init_header_buffer (char *(&header_buffer)[line_count])
    {
      for (size_t i = 0; i < line_count; ++i)
        {
          header_buffer[i] = 0;
        }
    }

    template <size_t line_count>
    void
    free_header_buffer (char *(&header_buffer)[line_count])
    {
      for (size_t i = 0; i < line_count; ++i)
        {
          delete[] header_buffer[i];
        }
    }
  }


  struct column_info 
  {
    size_t      w, w2;
    size_t      header_len;
    const char  *header;
    double      mult;

    column_info (size_t w, size_t w2, size_t len, const char *header, double mult = 1.0)
    : w (w), w2 (w2), header_len (len), header (header), mult (mult)
    {
    }
  };

  template <size_t N>
  column_info
  column (size_t width, size_t width2, const char (&header) [N])
  {
    return column_info (width, width2, N, header);
  }

  template <size_t N>
  column_info
  column (size_t w, size_t w2, const char (&header) [N], double mult)
  {
    return column_info (w, w2, N, header, mult);
  }

  const char *
  unit (const char *str)
  {
    return str;
  }

  template <typename ds_t, typename function_t>
  struct unit_t
  {
    const ds_t        &ds;
    const function_t  &function;
    const char        *note;

    unit_t (const ds_t &ds, const function_t &function, const char *note)
    : ds (ds)
    , function (function)
    , note (note)
    {
    }
  };

  template <typename ds_t, typename function_t>
  unit_t <ds_t, function_t>
  unit (const ds_t &ds, const function_t &function, const char *note)
  {
    return unit_t <ds_t, function_t> (ds, function, note);
  }

  struct column_t
  {
    enum {
      line_count = 4, 
    };

    const char  *line[line_count];
    size_t      line_idx;
    size_t      width;
    double      mult;

  public:
    char        string_format[15 + 1];
    char        float_format[15 + 1];
    char        *unit_format;

    char        *header;
    char        *unit;

    column_t ()
    {
      line_idx  = 0;
      width     = 0;

      memset (line,           0, sizeof (line));
      memset (string_format,  0, sizeof (string_format));
      memset (float_format,   0, sizeof (float_format));

      unit_format = 0;

      header    = 0;
      unit      = 0;
    }

    ~column_t ()
    {
      delete[] header;
      delete[] unit;
    }

    template <typename column_info_t>
    void
    make_header (const column_info_t &c)
    {
      width = c.w;
      header = new char [c.header_len + 1];
      memcpy (header, c.header, c.header_len);
      make_header (header, 0, c.header_len, c.w - 1);
    }

    void
    make_header (char *str, size_t space_idx, size_t len, size_t width)
    {
      static const char *empty_str = "";
      size_t idx = 0;

      line[line_idx] = str + space_idx;
      if (line_idx < line_count)
        line[line_idx + 1] = empty_str;

      for (size_t i = space_idx; i < len; ++i)
        {
          char c = str[i];
          if (c == ' ')
            {
              space_idx = i;
            }

          if (idx >= width)
            {
              str[space_idx] = 0;
              ++line_idx;
              make_header (str, space_idx + 1, len, width);
              break;
            }

          idx++;
        }
    }

    void
    operator () (const column_info &c, const char *unit_str)
    {
      mult = c.mult;
      sprintf (string_format, "| %%%ds ", (int)c.w - 2);
      sprintf (float_format, "|%%%d.%df", (int)c.w, (int)c.w2);
      make_header (c);

      unit_format = string_format;
      unit        = new char [strlen (unit_str) + 1];

      memset (unit, 0, strlen (unit_str) + 1);
      memcpy (unit, unit_str, strlen (unit_str) + 1);
    }

    template <typename ds_t, typename function_t>
    void
    operator () (const column_info &c, const unit_t <ds_t, function_t> &u)
    {
      mult = c.mult;
      sprintf (string_format, "| %%%ds ", (int)c.w - 2);
      sprintf (float_format, "|%%%d.%df", (int)c.w, (int)c.w2);
      make_header (c);

      unit_format   = float_format;
      double value  = u.function.data (*this, u.ds, mult, 1.0);

      if (value > 1.e9)
        {
          mult *= 1.e-6;
          unit = new char[strlen (u.note) + 1 + 1];
          unit[0] = 'M';
          memset (unit + 1, 0, strlen (u.note) + 1);
          memcpy (unit + 1, u.note, strlen (u.note));
        }
      else if (value > 1.e6)
        {
          mult *= 1.e-3;
          unit = new char[strlen (u.note) + 1 + 1];
          unit[0] = 'K';
          memset (unit + 1, 0, strlen (u.note) + 1);
          memcpy (unit + 1, u.note, strlen (u.note));
        }
      else
        {
          unit = new char[strlen (u.note) + 1];
          memset (unit, 0, strlen (u.note) + 1);
          memcpy (unit, u.note, strlen (u.note));
        }
    }
  };
  namespace detail {
    struct rate_rc
    {
      template <typename ds_t>
      static typename ds_t::rate_data_t
      get (const ds_t &ds)
      {
        return ds.rate_rc_;
      }
    };
    struct rate_rc_wefac
    {
      template <typename ds_t>
      static typename ds_t::rate_data_t
      get (const ds_t &ds)
      {
        return ds.rate_rc_wefac_;
      }
    };
    struct rate_total
    {
      template <typename ds_t>
      static typename ds_t::rate_data_t
      get (const ds_t &ds)
      {
        return ds.rate_total_;
      }
    };
    struct rate
    {
      template <typename ds_t>
      static typename ds_t::rate_data_t
      get (const ds_t &ds)
      {
        return ds.rate_;
      }
    };
    struct rate_wefac 
    {
      template <typename ds_t>
      static typename ds_t::rate_data_t
      get (const ds_t &ds)
      {
        return ds.rate_wefac_;
      }
    };

    struct prod_rate_data
    {
      template <typename data_t>
      static typename data_t::rate_data_inner
      get (const data_t &data)
      {
        return data.prod;
      }
    };

    struct inj_rate_data
    {
      template <typename data_t>
      static typename data_t::rate_data_inner
      get (const data_t &data)
      {
        return data.inj;
      }
    };
  }

  struct string_data
  {
    const char *
    format_string (const column_t &c) const
    {
      return c.string_format;
    }
  };

  struct float_data
  {
    const char *
    format_string (const column_t &c) const
    {
      return c.float_format;
    }
  };

  struct reservoir_name : string_data
  {
    template <typename data_source_t>
    const char *
    data (const column_t &c, const data_source_t &ds, float, float) const
    {
      return "RESERVOIR";
    }
  };

  struct empty_cell : string_data
  {
    template <typename data_source_t>
    const char *
    data (const column_t &c, const data_source_t &ds, float, float) const
    {
      return "";
    }
  };

  template <typename get_rate_data_t, typename get_inner_rate_data_t>
  struct liquid_rate : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float a_mult) const
    {
      const typename data_source_t::rate_data_t::rate_data_inner &rate = get_inner_rate_data_t::get (get_rate_data_t::get (ds));
      return rate.liquid * mult * a_mult;
    }
  };
  template <typename get_rate_data_t, typename get_inner_rate_data_t>
  struct oil_rate : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float a_mult) const
    {
      const typename data_source_t::rate_data_t::rate_data_inner &rate = get_inner_rate_data_t::get (get_rate_data_t::get (ds));
      return rate.oil * mult * a_mult;
    }
  };
  template <typename get_rate_data_t, typename get_inner_rate_data_t>
  struct water_rate : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float a_mult) const
    {
      const typename data_source_t::rate_data_t::rate_data_inner &rate = get_inner_rate_data_t::get (get_rate_data_t::get (ds));
      return rate.water * mult * a_mult;
    }
  };
  template <typename get_rate_data_t, typename get_inner_rate_data_t>
  struct gas_rate : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float a_mult) const
    {
      const typename data_source_t::rate_data_t::rate_data_inner &rate = get_inner_rate_data_t::get (get_rate_data_t::get (ds));
      return rate.gas * mult * a_mult;
    }
  };

  template <typename rate_t>
  struct prod_fluid_rate : float_data
  {
    float gas_liquid_rate_mult;

    prod_fluid_rate (float gas_liquid_rate_mult)
    : gas_liquid_rate_mult (gas_liquid_rate_mult)
    {
    }

    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float a_mult) const
    {
      return (rate_t::get (ds).free_gas * gas_liquid_rate_mult 
        + rate_t::get (ds).prod.water 
        + rate_t::get (ds).prod.oil) * mult * a_mult;
    }
  };
  template <typename rate_t>
  struct inj_fluid_rate : float_data
  {
    float gas_liquid_rate_mult;

    inj_fluid_rate (float gas_liquid_rate_mult)
    : gas_liquid_rate_mult (gas_liquid_rate_mult)
    {
    }

    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float a_mult) const
    {
      return (rate_t::get (ds).inj.gas * gas_liquid_rate_mult 
        + rate_t::get (ds).inj.water 
        + rate_t::get (ds).inj.oil) * mult * a_mult;
    }
  };

  template <typename get_rate_data_t, typename get_inner_rate_data_t>
  struct water_cut : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float, float) const
    {
      const typename data_source_t::rate_data_t::rate_data_inner &rate = get_inner_rate_data_t::get (get_rate_data_t::get (ds));
      return (rate.water * rate.oil) != 0.0 ? (rate.water / (rate.water + rate.oil)) : 0;
    }
  };

  template <typename get_rate_data_t, typename get_inner_rate_data_t>
  struct gas_oil_ratio : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float) const
    {
      const typename data_source_t::rate_data_t::rate_data_inner &rate = get_inner_rate_data_t::get (get_rate_data_t::get (ds));
      return fabs (rate.oil) > 0 ? (rate.gas * mult / rate.oil) : 0;
    }
  };

  template <typename get_rate_data_t, typename get_inner_rate_data_t>
  struct water_gas_ratio : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float) const
    {
      const typename data_source_t::rate_data_t::rate_data_inner &rate = get_inner_rate_data_t::get (get_rate_data_t::get (ds));
      return fabs (rate.gas) > 0 ? (rate.water * mult / rate.gas) : 0;
    }
  };

  struct wefac : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float, float) const
    {
      return ds.exploitation_factor_;
    }
  };

  struct well_name : string_data
  {
    template <typename data_source_t>
    const char *
    data (const column_t &c, const data_source_t &ds, float, float) const
    {
      return ds.name ().c_str ();
    }
  };

  struct ctrl_mode : string_data
  {
    template <typename data_source_t>
    const char *
    data (const column_t &c, const data_source_t &ds, float, float) const
    {
      if (ds.is_shut ())
        return "SHUT";

      if (ds.is_bhp ())
        {
          return "BHP";
        }
      else
        {
          if (ds.get_well_controller ()->is_production ())
            {
              using namespace wells;
              rate_control_type type = ds.get_well_controller ()->get_control_type ();
              switch (type)
              {
                case liquid_rate_control:
                  return "LRAT";
                case oil_rate_control:
                  return "ORAT";
                case water_rate_control:
                  return "WRAT";
                case gas_rate_control:
                  return "GRAT";
                default:
                  return "UNKNOWN";
              }
            }
          else
            {
              return "RATE";
            }
        }
    }
  };
  struct total_ctrl_mode : string_data
  {
    template <typename data_source_t>
    const char *
    data (const column_t &c, const data_source_t &ds, float, float) const
    {
      if (ds.is_shut ())
        return "SHUTED WELL";

      if (ds.get_well_controller ()->is_production ())
        {
          if (!ds.is_bhp ())
            {
              using namespace wells;
              rate_control_type type = ds.get_well_controller ()->get_control_type ();
              switch (type)
              {
                case liquid_rate_control:
                  return "LIQUID PRODUCER";
                case oil_rate_control:
                  return "OIL PRODUCER";
                case water_rate_control:
                  return "WATER PRODUCER";
                case gas_rate_control:
                  return "GAS PRODUCER";
                default:
                  return "UNKNOWN PRODUCER";
              }
            }
          else
            {
              return "BHP PRODUCER";
            }
        }
      else
        {
          if (ds.is_bhp ())
            {
              return "BHP INJECTOR";
            }
          else
            {
              using namespace wells;
              injection_type type = ds.get_well_controller ()->injection ();
              switch (type)
                {
                  case injection_oil:
                    return "OIL INJECTOR";
                  case injection_water:
                    return "WATER INJECTOR";
                  case injection_gas:
                    return "GAS INJECTOR";
                  default:
                    return "UNKNOWN INJECTOR";
                }
            }
        }
    }
  };

  struct well_gas_oil_ratio : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float a_mult) const
    {
      return fabs (ds.rate ().prod.oil) > 0 ? (ds.gor_ * mult) : 0;
    }
  };

  struct first_con_bhp : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float) const
    {
      return ds.get_connections_count () ? (ds.get_connection (0)->get_cur_bhp () * mult) : 0;
    }
  };

  struct first_con_bulkp : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float) const
    {
      return ds.get_connections_count () ? (ds.get_connection (0)->get_bulkp () * mult) : 0;
    }
  };

  struct grid_block 
  {
    mutable char buffer [128 * 3 + 1];
    mutable char format [128 * 3 + 1];

    const char *
    format_string (const column_t &c)
    {
      return "%s";
    }

    template <typename data_source_t>
    const char *
    data (const column_t &c, const data_source_t &ds, float mult, float) const
    {
      memset (buffer, 0, sizeof (buffer));
      memset (format, 0, sizeof (format));
      int width = (int)c.width - 2;
      int size  = width / 3;

      sprintf (format, "|%%%dd,%%%dd,%%%dd", size, size, width - size - size);
      sprintf (buffer, format, ds.i_coord () + 1, ds.j_coord () + 1, ds.k_coord () + 1);

      return buffer;
    }
  };

  struct con_bhp : float_data
  {
    float gravity;

    con_bhp (float gravity)
    : gravity (gravity)
    {
    }

    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float) const
    {
      return (ds.get_cur_bhp () + ds.get_density () * gravity * ds.get_connection_depth () * 0.5) * mult;
    }
  };

  struct con_bulkp : float_data
  {
    template <typename data_source_t>
    float
    data (const column_t &c, const data_source_t &ds, float mult, float) const
    {
      return ds.get_bulkp () * mult;
    }
  };

  namespace detail {
    template <typename holder_t>
    size_t 
    buffer_len (const holder_t &holder)
    {
      size_t len = 1;
      for (size_t i = 0; i < holder.columns_count_; ++i)
        {
          len += holder.columns_[i]->width + 1;
        }

      return len;
    }

    template <typename T>
    struct aptr
    {
      aptr (T *t = 0) : t (t) {}
      ~aptr () {delete[] t;}

      T *get () { return t; }
      T &operator[] (size_t i) { return t[i];}

      T *t;
    };

  }

  template <typename holder_t>
  void
  print_line (const holder_t &holder)
  {
    if (!holder.line_buffer)
      {
        size_t len = detail::buffer_len (holder);
        holder.line_buffer = new char [len + 1];
        memset (holder.line_buffer, '-', len);

        size_t idx = 0;
        for (size_t i = 0; i < holder.columns_count_; ++i)
          {
            holder.line_buffer[idx] = '|';
            idx += holder.columns_[i]->width + 1;
          }

        holder.line_buffer[idx + 0] = '|';
        holder.line_buffer[idx + 1] = 0;
      }

    BS_ASSERT (holder.line_buffer);
    BOSOUT (section::print_reports, level::highest) << holder.line_buffer << bs_end;
  }

  template <typename holder_t>
  void
  print_units (const holder_t &holder)
  {
    if (!holder.unit_buffer)
      {
        size_t len = detail::buffer_len (holder);
        holder.unit_buffer = new char [len + 1];
        memset (holder.unit_buffer, 0, len + 1);

        size_t idx = 0;
        for (size_t i = 0; i < holder.columns_count_; ++i)
          {
            sprintf (holder.unit_buffer + idx, holder.columns_[i]->string_format, holder.columns_[i]->unit);
            idx += holder.columns_[i]->width + 1;
          }
        holder.unit_buffer[idx + 0] = '|';
      }

    BOSOUT (section::print_reports, level::highest) << holder.unit_buffer << bs_end;
  }

  template <typename holder_t>
  void
  print_header (const holder_t &holder)
  {
    print_line (holder);

    if (!holder.header_buffer_count)
      {
        size_t line_idx = 0;
        for (size_t i = 0; i < holder.columns_count_; ++i)
          {
            if (holder.columns_[i]->line_idx > line_idx)
              line_idx = holder.columns_[i]->line_idx;
          }

        holder.header_buffer_count = line_idx;
        size_t len = detail::buffer_len (holder);
        for (size_t j = 0; j <= line_idx; ++j)
          {
            holder.header_buffer[j] = new char[len + 1];
            memset (holder.header_buffer[j], 0, len + 1);
            size_t idx = 0;
            for (size_t i = 0; i < holder.columns_count_; ++i)
              {
                sprintf (holder.header_buffer[j] + idx, holder.columns_[i]->string_format, holder.columns_[i]->line[j]);
                idx += holder.columns_[i]->width + 1;
              }
            holder.header_buffer[j][idx + 0] = '|';
          }
      }

    for (size_t j = 0; j <= holder.header_buffer_count; ++j)
      {
        BOSOUT (section::print_reports, level::highest) << holder.header_buffer[j] << bs_end;
      }

    print_units (holder);
    print_line (holder);
  }

  namespace detail {
    struct is_prod
    {
      template <typename T>
      inline bool
      filter (const T &well) const
      {
        return well->get_well_controller ()->is_production ();
      }
    };

    struct is_inj
    {
      template <typename T>
      inline bool
      filter (const T &well) const
      {
        return !well->get_well_controller ()->is_production ();
      }
    };

    template <typename printer_t, typename data_t, typename rs_t, typename well_filter_t>
    void
    print_well_data (const printer_t &printer, const data_t &data, const rs_t &rs, const well_filter_t &well_filter)
    {
      typedef typename rs_t::well_t                                     well_t;
      typedef typename rs_t::connection_t                               connection_t;
      typedef typename rs_t::facility_manager_t::well_const_iterator_t  well_iterator_t;

      typedef data_t                                                    sp_data_t;
      typedef typename rs_t::sp_well_t                                  sp_well_t;
      typedef typename rs_t::sp_connection_t                            sp_connection_t;

      print_header (printer.columns);
      printer.reservoir_data.print (printer.columns, data, rs);

      well_iterator_t wb = rs.get_facility_list ()->wells_begin ();
      well_iterator_t we = rs.get_facility_list ()->wells_end ();

      for (; wb != we; ++wb)
        {
          sp_well_t well (wb->second, bs_dynamic_cast ());

          if (well_filter.filter (well))
            {
              printer.well_data.print (printer.columns, data, *well);

              for (size_t j = 0, jcnt = well->get_connections_count (); j < jcnt; ++j)
                {
                  printer.connection_data.print (printer.columns, data, *well->get_connection (j), false);
                }
              if (well->get_connections_count ())
                print_line (printer.columns);
            }
        }
    }
  }

  template <typename strategy_t>
  struct prod_printer
  {
    typedef reservoir <strategy_t>                reservoir_t;
    typedef well <strategy_t>                     well_t;
    typedef wells::connection <strategy_t>        connection_t;
    typedef smart_ptr <idata, true>               sp_data_t;

    DECL_COLUMN_HOLDER (prod_columns_holder_t, 13, 
      (WELL_NAME,
        GRID_BLOCK,
        CTRL_MODE,
        WEFAC,
        OIL_RATE,
        WATER_RATE,
        GAS_RATE,
        FLUID_RES_VOL,
        WATER_CUT,
        GAS_OIL_RATIO,
        WATER_GAS_RATIO,
        BHP,
        BULK_PRESSURE));

    prod_columns_holder_t columns;

    typedef liquid_rate <detail::rate,            detail::prod_rate_data> liquid_rate_t;
    typedef oil_rate <detail::rate,               detail::prod_rate_data> oil_rate_t;
    typedef water_rate <detail::rate,             detail::prod_rate_data> water_rate_t;
    typedef gas_rate <detail::rate,               detail::prod_rate_data> gas_rate_t;
    typedef water_cut <detail::rate,              detail::prod_rate_data> water_cut_t;
    typedef gas_oil_ratio <detail::rate,          detail::prod_rate_data> gas_oil_ratio_t;
    typedef water_gas_ratio <detail::rate,        detail::prod_rate_data> water_gas_ratio_t;

    typedef oil_rate <detail::rate_wefac,         detail::prod_rate_data> oil_rate_wefac_t;
    typedef water_rate <detail::rate_wefac,       detail::prod_rate_data> water_rate_wefac_t;
    typedef gas_rate <detail::rate_wefac,         detail::prod_rate_data> gas_rate_wefac_t;
    typedef water_cut <detail::rate_wefac,        detail::prod_rate_data> water_cut_wefac_t;
    typedef gas_oil_ratio <detail::rate_wefac,    detail::prod_rate_data> gas_oil_ratio_wefac_t;
    typedef water_gas_ratio <detail::rate_wefac,  detail::prod_rate_data> water_gas_ratio_wefac_t;

    typedef prod_fluid_rate <detail::rate_rc>                             fluid_rate_t;
    typedef prod_fluid_rate <detail::rate_rc_wefac>                       fluid_rate_rc_wefac_t;

    prod_printer (const sp_data_t &data, const smart_ptr <reservoir_t, true> &rs)
    {
      columns.WELL_NAME           (column (12, 0, "WELL NAME"),                               unit (""));
      columns.GRID_BLOCK          (column (14, 0, "LOCATION (I,J,K)"),                        unit (""));
      columns.CTRL_MODE           (column (12, 0, "CTRL MODE"),                               unit (""));
      columns.WEFAC               (column (12, 5, "WEFAC"),                                   unit (""));
      columns.OIL_RATE            (column (13, 3, "OIL RATE", LIQUID_RATE_MULT),              unit (*rs, oil_rate_wefac_t (), RATE_LIQUID_SURFACE_NOTE));
      columns.WATER_RATE          (column (13, 3, "WATER RATE", LIQUID_RATE_MULT),            unit (*rs, water_rate_wefac_t (), RATE_LIQUID_SURFACE_NOTE));
      columns.GAS_RATE            (column (13, 3, "GAS RATE", GAS_RATE_MULT),                 unit (*rs, gas_rate_wefac_t (), RATE_GAS_SURFACE_NOTE));
      columns.FLUID_RES_VOL       (column (16, 3, "FLUID RES.VOL.", LIQUID_RATE_MULT),        unit (*rs, fluid_rate_rc_wefac_t (GAS_LIQUID_RATE_MULT), RATE_LIQUID_RESERVOIR_NOTE));
      columns.WATER_CUT           (column (15, 3, "WATER CUT"),                               unit (""));
      columns.GAS_OIL_RATIO       (column (15, 3, "GAS/OIL RATIO", GAS_OIL_RATIO_MULT),       unit (*rs, gas_oil_ratio_wefac_t (), GAS_LIQUID_NOTE));
      columns.WATER_GAS_RATIO     (column (15, 4, "WATER/GAS RATIO", WATER_GAS_RATIO_MULT),   unit (*rs, water_gas_ratio_wefac_t (), LIQUID_GAS_NOTE));
      columns.BHP                 (column (15, 2, "BHP OR CON.PR.", PRESSURE_MULT),           unit (PRESSURE_NOTE));
      columns.BULK_PRESSURE       (column (12, 2, "Bulk Pressure", PRESSURE_MULT),            unit (PRESSURE_NOTE));
    }

    DATA_LINE (prod_columns_holder_t, reservoir_data, sp_data_t, reservoir_t, 
      ((WELL_NAME,       reservoir_name (),                              1.0))
      ((GRID_BLOCK,      empty_cell (),                                  1.0))
      ((CTRL_MODE,       empty_cell (),                                  1.0))
      ((WEFAC,           empty_cell (),                                  1.0))
      ((OIL_RATE,        oil_rate_wefac_t (),                           -1.0))
      ((WATER_RATE,      water_rate_wefac_t (),                         -1.0))
      ((GAS_RATE,        gas_rate_wefac_t (),                           -1.0))
      ((FLUID_RES_VOL,   fluid_rate_rc_wefac_t (GAS_LIQUID_RATE_MULT),  -1.0))
      ((WATER_CUT,       water_cut_wefac_t (),                           1.0))
      ((GAS_OIL_RATIO,   gas_oil_ratio_wefac_t (),                       1.0))
      ((WATER_GAS_RATIO, water_gas_ratio_wefac_t (),                     1.0))
      ((BHP,             empty_cell (),                                  1.0))
      ((BULK_PRESSURE,   empty_cell (),                                  1.0))
    );

    DATA_LINE (prod_columns_holder_t, well_data, sp_data_t, well_t, 
      ((WELL_NAME,       well_name (),                                   1.0))
      ((GRID_BLOCK,      empty_cell (),                                  1.0))
      ((CTRL_MODE,       ctrl_mode (),                                   1.0))
      ((WEFAC,           wefac (),                                       1.0))
      ((OIL_RATE,        oil_rate_t (),                                 -1.0))
      ((WATER_RATE,      water_rate_t (),                               -1.0))
      ((GAS_RATE,        gas_rate_t (),                                 -1.0))
      ((FLUID_RES_VOL,   fluid_rate_t (GAS_LIQUID_RATE_MULT),           -1.0))
      ((WATER_CUT,       water_cut_t (),                                 1.0))
      ((GAS_OIL_RATIO,   well_gas_oil_ratio (),                          1.0))
      ((WATER_GAS_RATIO, water_gas_ratio_t (),                           1.0))
      ((BHP,             first_con_bhp (),                               1.0))
      ((BULK_PRESSURE,   first_con_bulkp (),                             1.0))
    );

    DATA_LINE (prod_columns_holder_t, connection_data, sp_data_t, connection_t, 
      ((WELL_NAME,       empty_cell (),                                  1.0))
      ((GRID_BLOCK,      grid_block (),                                  1.0))
      ((CTRL_MODE,       empty_cell (),                                  1.0))
      ((WEFAC,           empty_cell (),                                  1.0))
      ((OIL_RATE,        oil_rate_t (),                                 -1.0))
      ((WATER_RATE,      water_rate_t (),                               -1.0))
      ((GAS_RATE,        gas_rate_t (),                                 -1.0))
      ((FLUID_RES_VOL,   fluid_rate_t (GAS_LIQUID_RATE_MULT),           -1.0))
      ((WATER_CUT,       water_cut_t (),                                 1.0))
      ((GAS_OIL_RATIO,   gas_oil_ratio_t (),                             1.0))
      ((WATER_GAS_RATIO, water_gas_ratio_t (),                           1.0))
      ((BHP,             con_bhp (GRAVITY),                              1.0))
      ((BULK_PRESSURE,   con_bulkp (),                                   1.0))
    );
  };

  template <typename strategy_t>
  struct inj_printer
  {
    typedef reservoir <strategy_t>                reservoir_t;
    typedef well <strategy_t>                     well_t;
    typedef wells::connection <strategy_t>        connection_t;
    typedef smart_ptr <idata, true>               sp_data_t;

    DECL_COLUMN_HOLDER (inj_columns_holder_t, 10,
      (
        WELL_NAME,
        GRID_BLOCK,
        CTRL_MODE,
        WEFAC,
        WATER_RATE,
        OIL_RATE,
        GAS_RATE,
        FLUID_RES_VOL,
        BHP,
        BULK_PRESSURE
      ));

    inj_columns_holder_t columns;

    typedef oil_rate <detail::rate,       detail::inj_rate_data>  oil_rate_t;
    typedef water_rate <detail::rate,     detail::inj_rate_data>  water_rate_t;
    typedef gas_rate <detail::rate,       detail::inj_rate_data>  gas_rate_t;

    typedef oil_rate <detail::rate_wefac, detail::inj_rate_data>  oil_rate_wefac_t;
    typedef oil_rate <detail::rate_wefac, detail::inj_rate_data>  water_rate_wefac_t;
    typedef oil_rate <detail::rate_wefac, detail::inj_rate_data>  gas_rate_wefac_t;

    typedef inj_fluid_rate <detail::rate_rc>                      fluid_rate_t;
    typedef inj_fluid_rate <detail::rate_rc_wefac>                fluid_rate_rc_wefac_t;

    inj_printer (const sp_data_t &data, const smart_ptr <reservoir_t, true> &rs)
    {
      columns.WELL_NAME           (column (12, 0, "WELL NAME"),                               unit (""));
      columns.GRID_BLOCK          (column (15, 0, "LOCATION (I,J,K)"),                        unit (""));
      columns.CTRL_MODE           (column (20, 0, "CTRL MODE"),                               unit (""));
      columns.WEFAC               (column (13, 5, "WEFAC"),                                   unit (""));
      columns.WATER_RATE          (column (13, 3, "WATER RATE", LIQUID_RATE_MULT),            unit (*rs, water_rate_wefac_t (), RATE_LIQUID_SURFACE_NOTE));
      columns.OIL_RATE            (column (14, 3, "OIL RATE", LIQUID_RATE_MULT),              unit (*rs, oil_rate_wefac_t (), RATE_LIQUID_SURFACE_NOTE));
      columns.GAS_RATE            (column (14, 3, "GAS RATE", GAS_RATE_MULT),                 unit (*rs, gas_rate_wefac_t (), RATE_GAS_SURFACE_NOTE));
      columns.FLUID_RES_VOL       (column (16, 3, "FLUID RES.VOL.", LIQUID_RATE_MULT),        unit (*rs, fluid_rate_rc_wefac_t (GAS_LIQUID_RATE_MULT), RATE_LIQUID_RESERVOIR_NOTE));
      columns.BHP                 (column (15, 2, "BHP OR CON.PR.", PRESSURE_MULT),           unit (PRESSURE_NOTE));
      columns.BULK_PRESSURE       (column (13, 2, "Bulk Pressure", PRESSURE_MULT),            unit (PRESSURE_NOTE));
    }

    DATA_LINE (inj_columns_holder_t, reservoir_data, sp_data_t, reservoir_t,
      ((WELL_NAME,       reservoir_name (),                             1.0))
      ((GRID_BLOCK,      empty_cell (),                                 1.0))
      ((CTRL_MODE,       empty_cell (),                                 1.0))
      ((WEFAC,           empty_cell (),                                 1.0))
      ((WATER_RATE,      water_rate_wefac_t (),                         1.0))
      ((OIL_RATE,        oil_rate_wefac_t (),                           1.0))
      ((GAS_RATE,        gas_rate_wefac_t (),                           1.0))
      ((FLUID_RES_VOL,   fluid_rate_rc_wefac_t (GAS_LIQUID_RATE_MULT),  1.0))
      ((BHP,             empty_cell (),                                 1.0))
      ((BULK_PRESSURE,   empty_cell (),                                 1.0))
    );

    DATA_LINE (inj_columns_holder_t, well_data, sp_data_t, well_t,
      ((WELL_NAME,       well_name (),                                  1.0))
      ((GRID_BLOCK,      empty_cell (),                                 1.0))
      ((CTRL_MODE,       ctrl_mode (),                                  1.0))
      ((WEFAC,           wefac (),                                      1.0))
      ((WATER_RATE,      water_rate_t (),                               1.0))
      ((OIL_RATE,        oil_rate_t (),                                 1.0))
      ((GAS_RATE,        gas_rate_t (),                                 1.0))
      ((FLUID_RES_VOL,   fluid_rate_t (GAS_LIQUID_RATE_MULT),           1.0))
      ((BHP,             first_con_bhp (),                              1.0))
      ((BULK_PRESSURE,   first_con_bulkp (),                            1.0))
    );

    DATA_LINE (inj_columns_holder_t, connection_data, sp_data_t, connection_t,
      ((WELL_NAME,       empty_cell (),                                 1.0))
      ((GRID_BLOCK,      grid_block (),                                 1.0))
      ((CTRL_MODE,       empty_cell (),                                 1.0))
      ((WEFAC,           empty_cell (),                                 1.0))
      ((WATER_RATE,      water_rate_t (),                               1.0))
      ((OIL_RATE,        oil_rate_t (),                                 1.0))
      ((GAS_RATE,        gas_rate_t (),                                 1.0))
      ((FLUID_RES_VOL,   fluid_rate_t (GAS_LIQUID_RATE_MULT),           1.0))
      ((BHP,             con_bhp (GRAVITY),                             1.0))
      ((BULK_PRESSURE,   con_bulkp (),                                  1.0))
    );
  };

  template <typename strategy_t>
  struct prod_total_printer
  {
    typedef reservoir <strategy_t>                reservoir_t;
    typedef well <strategy_t>                     well_t;
    typedef wells::connection <strategy_t>        connection_t;
    typedef smart_ptr <idata, true>               sp_data_t;
    typedef smart_ptr <well_t, true>              sp_well_t;

    DECL_COLUMN_HOLDER (columns_holder, 9,
      (
        WELL_NAME,
        CTRL_MODE,
        OIL_PROD,
        WATER_PROD,
        GAS_PROD,
        LIQUID_PROD,
        WATER_INJ,
        OIL_INJ,
        GAS_INJ
      ));

    columns_holder columns;

    typedef oil_rate <detail::rate_total,     detail::prod_rate_data>  prod_oil_rate_t;
    typedef water_rate <detail::rate_total,   detail::prod_rate_data>  prod_water_rate_t;
    typedef gas_rate <detail::rate_total,     detail::prod_rate_data>  prod_gas_rate_t;
    typedef liquid_rate <detail::rate_total,  detail::prod_rate_data>  prod_liquid_rate_t;

    typedef oil_rate <detail::rate_total,     detail::inj_rate_data>   inj_oil_rate_t;
    typedef water_rate <detail::rate_total,   detail::inj_rate_data>   inj_water_rate_t;
    typedef gas_rate <detail::rate_total,     detail::inj_rate_data>   inj_gas_rate_t;

    prod_total_printer (const sp_data_t &data, const smart_ptr <reservoir_t, true> &rs)
    {
      columns.WELL_NAME   (column (12, 0, "WELL NAME"),                                     unit (""));
      columns.CTRL_MODE   (column (22, 0, "Current CTRL MODE and status:"),                 unit (""));
      columns.OIL_PROD    (column (17, 0, "Total Oil production", LIQUID_RATE_MULT),        unit (*rs, prod_oil_rate_t    (), VOL_LIQUID_SURFACE_NOTE));
      columns.WATER_PROD  (column (15, 0, "Total Water production", LIQUID_RATE_MULT),      unit (*rs, prod_water_rate_t  (), VOL_LIQUID_SURFACE_NOTE));
      columns.GAS_PROD    (column (20, 0, "Total Gas production", GAS_RATE_MULT),           unit (*rs, prod_gas_rate_t    (), VOL_GAS_SURFACE_NOTE));
      columns.LIQUID_PROD (column (24, 0, "Total Reser Liq. production", LIQUID_RATE_MULT), unit (*rs, prod_liquid_rate_t (), VOL_LIQUID_RES_NOTE));
      columns.WATER_INJ   (column (20, 0, "Total Water injection", LIQUID_RATE_MULT),       unit (*rs, inj_water_rate_t   (), VOL_LIQUID_SURFACE_NOTE));
      columns.OIL_INJ     (column (20, 0, "Total Oil injection", LIQUID_RATE_MULT),         unit (*rs, inj_oil_rate_t     (), VOL_LIQUID_SURFACE_NOTE));
      columns.GAS_INJ     (column (20, 0, "Total Gas injection", GAS_RATE_MULT),            unit (*rs, inj_gas_rate_t     (), VOL_GAS_SURFACE_NOTE));
    }

    DATA_LINE (columns_holder, reservoir_data, sp_data_t, reservoir_t, 
      ((WELL_NAME,      reservoir_name (),                          1.0))
      ((CTRL_MODE,      empty_cell (),                              1.0))
      ((OIL_PROD,       prod_oil_rate_t (),                        -1.0))
      ((WATER_PROD,     prod_water_rate_t (),                      -1.0))
      ((GAS_PROD,       prod_gas_rate_t (),                        -1.0))
      ((LIQUID_PROD,    prod_liquid_rate_t (),                     -1.0))
      ((WATER_INJ,      inj_water_rate_t (),                        1.0))
      ((OIL_INJ,        inj_oil_rate_t (),                          1.0))
      ((GAS_INJ,        inj_gas_rate_t (),                          1.0))
    );

    DATA_LINE (columns_holder, well_data, sp_data_t, well_t,
      ((WELL_NAME,      well_name (),                               1.0))
      ((CTRL_MODE,      total_ctrl_mode (),                         1.0))
      ((OIL_PROD,       prod_oil_rate_t (),                        -1.0))
      ((WATER_PROD,     prod_water_rate_t (),                      -1.0))
      ((GAS_PROD,       prod_gas_rate_t (),                        -1.0))
      ((LIQUID_PROD,    prod_liquid_rate_t (),                     -1.0))
      ((WATER_INJ,      inj_water_rate_t (),                        1.0))
      ((OIL_INJ,        inj_oil_rate_t (),                          1.0))
      ((GAS_INJ,        inj_gas_rate_t (),                          1.0))
    );

    void
    print (const sp_data_t &data, const reservoir_t &rs)
    {
      print_header (columns);
      reservoir_data.print (columns, data, rs);

      typedef typename reservoir_t::facility_manager_t::well_const_iterator_t  well_iterator_t;

      well_iterator_t wb = rs.get_facility_list ()->wells_begin ();
      well_iterator_t we = rs.get_facility_list ()->wells_end ();

      for (; wb != we; ++wb)
        {
          sp_well_t well (wb->second, bs_dynamic_cast ());
          well_data.print (columns, data, *well);
        }
    }
  };

  template <typename strategy_t>
  void
  well_data_printer <strategy_t>::print_prod (const smart_ptr <idata, true> &data, const smart_ptr <reservoir <strategy_t>, true> &rs)
  {
    static prod_printer <strategy_t> prod_printer_ (data, rs);

    BOSOUT (section::print_reports, level::highest) <<
      "    PRODUCTION REPORT\n"
      "    -----------------" 
      << bs_end;
    detail::print_well_data (prod_printer_, data, *rs, detail::is_prod ());
  }

  template <typename strategy_t>
  void
  well_data_printer <strategy_t>::print_inj (const smart_ptr <idata, true> &data, const smart_ptr <reservoir <strategy_t>, true> &rs)
  {
    static inj_printer <strategy_t> inj_printer_ (data, rs);

    BOSOUT (section::print_reports, level::highest) <<
      "    INJECTOR REPORT\n"
      "    ---------------" 
      << bs_end;
    detail::print_well_data (inj_printer_, data, *rs, detail::is_inj ());
  }

  template <typename strategy_t>
  void
  well_data_printer <strategy_t>::print_total_prod (const smart_ptr <idata, true> &data, const smart_ptr <reservoir <strategy_t>, true> &rs)
  {
    static prod_total_printer <strategy_t> prod_total_printer_ (data, rs);

    BOSOUT (section::print_reports, level::highest) <<
      "    CUMULATIVE INJECTION/PRODUCTION TOTAL\n"
      "    -------------------------------------" 
      << bs_end;
    prod_total_printer_.print (data, *rs);
  }

  template struct well_data_printer <base_strategy_fi>;
  template struct well_data_printer <base_strategy_di>;
  template struct well_data_printer <base_strategy_mixi>;

} // namespace blue_sky

