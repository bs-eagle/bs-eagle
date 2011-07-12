/**
 * \file pvt_base.h
 * \brief
 * \author Miryanov Sergey
 * \date 29.04.2008
 */
#ifndef BS_PVT_BASE_H_
#define BS_PVT_BASE_H_

namespace blue_sky
  {

  namespace pvt
    {

    template <typename lhs_t, typename rhs_t>
    inline void
    set_pvt_pointer (lhs_t *dst, rhs_t value)
    {
      if (dst)
        {
          *dst = value;
        }
    }

  }	// namespace pvt

  /**
   * \brief pvt_base
   */
  class BS_API_PLUGIN pvt_base : public objbase
    {
      friend class table;
    public:

      typedef table_iface::vector_t                   vector_t;
      typedef smart_ptr<table_iface,true>             sp_table; 
      /**
       * \brief constructor
       */
      pvt_base ();

      /**
       * \brief destructor
       */
      virtual ~pvt_base () {}

      /**
       * \brief parse string and store values into data
       *
       * \param char_line line of chars
       */
      //virtual void parse_char_line (const char * char_line) = 0;

      /**
       * \brief store values into data
       *
       * \param seq_vector
       */
      virtual void insert_vector (const v_double &vec) = 0;

      /**
       * \brief set density and molar density
       *
       * \param density
       * \param molar density
       * */
      virtual void set_density (t_double density, t_double molar_density);

      /**
       * \brief build dependent data
       *
       * \param atm_p
       * \param min_p
       * \param max_p
       * \param n_intervals
       * \return
       */
      virtual void build (t_double atm_p, t_double min_p, t_double max_p, t_long n_intervals) = 0;

      /**
       * \brief
       * */
      virtual t_double interpolate_and_fix (t_double cell_pbub) const;

      /**
       * \brief
       * */
      virtual t_double get_gor_for_pressure (t_double pressure_data) const;

      //! get p_step value
      t_double get_p_step () const;

      //! get surface_density value
      t_double get_surface_density () const;

      //! set surface density
      void set_surface_density (t_double surface_density);

      //! print pvt table
      virtual void print () const = 0;
      
      //! return pvt_input 
      sp_table get_pvt_input_table () const
        {  return pvt_input_props; }
 
    protected:

      /**
       * \brief check pressure interval and throw exception if interval too small
       *
       * \param min_p minimal value of pressure
       * \param max_p maximal value of pressure
       * \return
       */
      void check_pressure_interval (t_double min_p, t_double max_p);

      /**
       * \brief check number of intervals and if n_intervals is invalid correct its value
       *
       * \param n_intervals number of intervals
       * \return
       */
      void check_interval_numbers (t_long &n_intervals);

      void check_common ();

      void check_gas_common (const vector_t &pressure, const vector_t &fvf, const vector_t &visc);

      void check_oil_common (const vector_t &pressure, const vector_t &fvf, const vector_t &visc);

    protected:

      //! interpolation step (usually step of PRESSURE)
      t_double p_step;

      //! density of phase at surface condition
      t_double surface_density;

      t_double molar_density;

      //!
      bool init_dependent;
    public:
      //! input main pvt properties 
      sp_table pvt_input_props; 
      //! pvt properties in table format
      sp_table pvt_props_table;
      
    };

  //! register all pvt_* types
  bool pvt_register_types (const plugin_descriptor &pd);

} // namespace blue_sky

#endif  // #ifndef BS_PVT_BASE_H_
