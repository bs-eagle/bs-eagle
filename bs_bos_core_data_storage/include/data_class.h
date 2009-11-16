#ifndef BS_DATA_CLASS_H
#define BS_DATA_CLASS_H
/*!
  \file data_class.h
  \brief initial data holder
	\author Nikonov Maxim
*/

#include "pool.h"
#include "rocktab_table.h"
#include "prvd_table.h"
#include "read_class.h"
#include "arrays.h"
#include "arrays_tables.h"
#include "convert_units.h"

#include "data_dimens.h"

namespace blue_sky {

//!< values of geometry definition type
#define GEOM_FLAG_NOT_INITIALIZED  0
#define GEOM_FLAG_DX_DY_DZ         1
#define GEOM_FLAG_ZCORN_COORD      2

#define BOUNDARY_CONDITION_TYPE_I  1
#define BOUNDARY_CONDITION_TYPE_II 2

#define DEFAULT_BOUNDARY_CONDITION 2

//! Maximal number of time steps in one TSTEP expression
#define MAX_TIME_STEPS_DEF      1000

  /*!
  \brief convert  -- convert 'carray' to UfaSolver format and write it to 'array'
  \param msh    -- mesh
  \param array  -- pointer to destination array
  \param carray -- pointer to source array
  */
  template <class index_array_t, class array_t, class carray_t>
  void
  convert_arrays (typename index_array_t::value_type cells_count, const index_array_t &index_map, array_t &array, const carray_t &carray)
  {
    for (typename index_array_t::value_type i = 0; i < cells_count; ++i)
      array[i] = carray[index_map[i]];
  }

  /*!
  \class data
  \ingroup KeywordLanguage
  \brief Main class for read, save and store input information
  */
  class BS_API_PLUGIN idata : public bs_node
    {
    private:
      struct idata_traits;

    public:
      //! typedefs
      typedef idata                                       this_t;
      typedef smart_ptr<this_t , true>                    sp_this_t;

      typedef bs_array_map <int, uint8_t>                 amap_i;
      typedef bs_array_map <int, float16_t>               amap_d;

      typedef smart_ptr<amap_i>                           sp_amap_i;
      typedef smart_ptr<amap_d>                           sp_amap_d;

      typedef seq_vector <val_vs_depth>                   vval_vs_depth;
      typedef seq_vector <float16_t>                      vec_d;
      typedef seq_vector <uint8_t>                        vec_i;

      typedef smart_ptr <FRead, true>										  sp_reader_t;

      struct pvt_info
        {
          typedef seq_vector <double> main_pvt_data_t;

          main_pvt_data_t							main_data_;
          auto_value <bool, false>		has_density_;
          auto_value <double>         density_;
          auto_value <double>         molar_density_;
        };

      typedef seq_vector <pvt_info>											pvt_vector;

			struct BS_API_PLUGIN arrays_helper {
				typedef std::map < std::string, int > names_map_t;

				names_map_t names_map_i;
				names_map_t names_map_d;

				uint8_t     *dummy_array_i;
				float16_t   *dummy_array_d;

				arrays_helper ();
				~arrays_helper ();

				void init_names_maps ();

				void add_correspondence_i (const std::string &name, int index);
				void add_correspondence_d (const std::string &name, int index);

				int get_idx_i (const std::string &name) const;
				int get_idx_d (const std::string &name) const;
			};

		public:

      ~idata ();

      //! \brief init idata members method
      void init();

      //! \brief Assignment operator
      this_t &operator=(const this_t&);

      void set_defaults_in_pool();
      void set_region (int r_pvt, int r_sat, int r_eql, int r_fip);

      void set_density (const seq_vector <double> &density);
      void set_density_internal (const double *density);

      //! return prvd vector
      vval_vs_depth &get_prvd();
      //! return rsvd vector
      vval_vs_depth &get_rsvd();
      //! return pbvd vector
      vval_vs_depth &get_pbvd();

      vec_d &
      get_rock()
      {
        return rock;
      }

      vec_d &
      get_p_ref()
      {
        return p_ref;
      }

      //int no_blanks (char *buf) const;
      //void build_argument_list (const sp_reader_t &reader);
			//void build_operator_list ();
			//void output_argument_list ();
			//void clear_argument_list ();
      //void read_left (const sp_reader_t &reader,
      //                char *s, char **right, std::string &name, int *i1,
      //                int *i2, int *j1, int *j2, int *k1, int *k2, const char *keyword);
      //void algorithm_read_and_done (const sp_reader_t &reader, char *buf, const char *keyword);
      //void calculate (const sp_reader_t &reader, ar_args_t &l, char *right, const char *keyword);
      //int read_numbers (ar_stack <ar_args_t> &st, char **buf);
      //int read_operator (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op,
      //                    char **start_ptr, int *flag_brek, const char *keyword,
      //                    const char *buf, int *Arguments_count);
      //int do_calculate (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op, const char *keyword);
      //int prior_calculate (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op, int priority, const char *keyword);
      //int read_arg_func (const sp_reader_t &reader, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op, char **start_ptr, const char *keyword, int *arguments_count, int flag_brek);
      //int test_token (int prev, int cur);

      array_uint8_t   get_int_non_empty_array (int array_index);
      array_float16_t get_float_non_empty_array (int array_index);
      
      array_uint8_t   get_int_array (const std::string &array_name);
      array_float16_t get_float_array (const std::string &array_name);
      
      array_uint8_t & get_int_non_empty_array (const std::string &array_name);
      array_float16_t get_float_non_empty_array (const std::string &array_name);
      
      bool 
      contain (const std::string &array_name) const;
      
    public:
			arrays_helper ahelper;
      int rpo_model;                //!< 3-ph oil relative permeability model: flag 0, 1 or 2 (stone model)

      double minimal_pore_volume;     //!< Minimal pore volume allowed for active cells
      double minimal_splice_volume;   //!< Minimal pore volume allowed for active cells to splice with other cells
      double maximum_splice_thickness;//!< Default maximum thickness allowed between active cells to be coupled

      data_dimens dimens;                                  //!< dimension description

      auto_value <long long> pvt_region;               //!< Number of PVT regions in simulation
      auto_value <long long> sat_region;               //!< Number of saturation regions in simulation
      auto_value <long long> eql_region;               //!< Number of equilibrium regions in simulation
      auto_value <long long> fip_region;               //!< Number of FIP regions in simulation

      auto_value <int> fi_n_phase;               //!< number of phases if full implicit simulation
      auto_value <int> fi_phases;                //!< sizeof (int) bit fields (1 -- phase present, 0 -- do not present)

      auto_value <long long> rock_region;              //!< Number of ROCK regions

      vec_i equil_regions;

      //! \brief class's public data area
      sp_amap_i   i_map;
      sp_amap_d   d_map;


      auto_value <int> init_section;             //!< flag indicating whether we have init section
      //!< if init_section == 0, cdata::check_sat will not go.

      convert_units input_units_converter;          //!< Input units converter
      convert_units output_units_converter;         //!< Output units converter

      std::vector<rocktab_table <base_strategy_fi> > rocktab;    //!< rocktab tables for all rock regions

      //TITLE
      std::string title;

      pvt_vector pvto, pvtdo, pvtg, pvtw;

      vec_d rock;            //!< Array (pvt_region) - compressibility of rock for each region
      vec_d equil;
      vec_d p_ref;           //!< Array (pvt_region) - reference pressure for compressibility of rock for each region

      //! pressure points at reference depth used for PRVD keyword content,
      //! and initial pressure initialization
      //! for all pvt regions
      //! array (eql_region)
      vval_vs_depth prvd;

      //! RS points at reference depth used for RSVD keyword content,
      //! and initial RS initialization
      //! for all pvt regions
      //! array (eql_region)
      vval_vs_depth rsvd;

      //! PBUB points at reference depth used for PBVD keyword content,
      //! and initial RS initialization
      //! for all pvt regions
      //! array (eql_region)
      vval_vs_depth pbvd;

      BLUE_SKY_TYPE_DECL_T(idata)
    };

  bool register_idata(const plugin_descriptor &pd);
}

#endif // BS_DATA_CLASS_H
