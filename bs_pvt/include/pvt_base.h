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
    public:

      typedef t_double                              	item_t;
      typedef t_long                        			  	index_t;
      typedef v_double                                array_item_t;
      typedef v_long                                  array_index_t;
      typedef smart_ptr <v_double, true>       			  sp_array_item_t;
      typedef smart_ptr <v_long, true>                sp_array_index_t;
      typedef v_double                   		          input_vector_t;

      typedef std::vector <t_double>                  vector_t;

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
      virtual void insert_vector (const input_vector_t &vec) = 0;

      /**
       * \brief set density and molar density
       *
       * \param density
       * \param molar density
       * */
      virtual void set_density (item_t density, item_t molar_density);

      /**
       * \brief build dependent data
       *
       * \param atm_p
       * \param min_p
       * \param max_p
       * \param n_intervals
       * \return
       */
      virtual void build (item_t atm_p, item_t min_p, item_t max_p, t_int n_intervals) = 0;

      /**
       * \brief
       * */
      virtual item_t interpolate_and_fix (item_t cell_pbub) const;

      /**
       * \brief
       * */
      virtual item_t get_gor_for_pressure (item_t pressure_data) const;

      //! get p_step value
      item_t get_p_step () const;

      //! get surface_density value
      item_t get_surface_density () const;

      //! set surface density
      void set_surface_density (item_t surface_density);

      //! print pvt table
      virtual void print () const = 0;

    protected:

      /**
       * \brief check pressure interval and throw exception if interval too small
       *
       * \param min_p minimal value of pressure
       * \param max_p maximal value of pressure
       * \return
       */
      void check_pressure_interval (item_t min_p, item_t max_p);

      /**
       * \brief check number of intervals and if n_intervals is invalid correct its value
       *
       * \param n_intervals number of intervals
       * \return
       */
      void check_interval_numbers (t_int &n_intervals);

      void check_common ();

      void check_gas_common (const vector_t &pressure, const vector_t &fvf, const vector_t &visc);

      void check_oil_common (const vector_t &pressure, const vector_t &fvf, const vector_t &visc);

    protected:

      //! interpolation step (usually step of PRESSURE)
      item_t p_step;

      //! density of phase at surface condition
      item_t surface_density;

      item_t molar_density;

      //!
      bool init_dependent;
    };

  //! register all pvt_* types
  bool pvt_register_types (const plugin_descriptor &pd);

} // namespace blue_sky

#endif  // #ifndef BS_PVT_BASE_H_
