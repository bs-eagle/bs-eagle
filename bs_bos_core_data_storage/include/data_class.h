#ifndef BS_DATA_CLASS_H
#define BS_DATA_CLASS_H
/*!
  \file data_class.h
  \brief initial data holder
	\author Nikonov Maxim
*/
#include BS_FORCE_PLUGIN_IMPORT ()

#include "pool_iface.h"
#include "prop_iface.h"

#include BS_STOP_PLUGIN_IMPORT ()

#include "rocktab_table.h"
#include "prvd_table.h"
#include "read_class.h"
#include "arrays.h"
#include "arrays_tables.h"
#include "convert_units.h"
#include "auto_value.h"

#include "data_dimens.h"

#include "strategies.h"
#include "bs_tree.h"


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

enum   //! indexes for dimension parameters
  {
    ARRAY_POOL_NX_A,
    ARRAY_POOL_NX_B,
    ARRAY_POOL_NY_A,
    ARRAY_POOL_NY_B,
    ARRAY_POOL_NZ_A,
    ARRAY_POOL_NZ_B,
    ARRAY_POOL_TOTAL
  };

  /*!
  \brief convert  -- convert 'carray' to UfaSolver format and write it to 'array'
  \param msh    -- mesh
  \param array  -- pointer to destination array
  \param carray -- pointer to source array
  */
  template <class index_t, class sp_index_array_t, class array_t, class carray_t>
  void
  convert_arrays (index_t cells_count, const sp_index_array_t index_map, array_t &array, const carray_t &carray)
  {
    index_t *index_map_data = &(*index_map)[0];
    //typename carray_t::pointed_t::value_type *carray_data = &(*carray)[0];
    
    for (index_t i = 0; i < cells_count; ++i)
      array[i] = carray[index_map_data[i]];
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
      typedef idata                           this_t;
      typedef smart_ptr<this_t , true>                    sp_this_t;
      
      typedef std::vector <val_vs_depth>                  vval_vs_depth;
      
      typedef smart_ptr <FRead, true>										  sp_reader_t;
      typedef smart_ptr <h5_pool_iface, true>							sp_h5_pool_t;
      typedef smart_ptr <prop_iface, true>							  sp_prop_t;

      struct pvt_info
        {
          spv_float							      main_data_;
          auto_value <bool, false>		has_density_;
          auto_value <t_float>      density_;
          auto_value <t_float>      molar_density_;
        };

      typedef std::vector <pvt_info>									  pvt_vector;

		public:

      ~idata ();

      //! \brief init idata members method
      void init();

      //! \brief flush pool data
      void flush_pool();

      //! \brief Assignment operator
      this_t &operator=(const this_t&);

      void set_region (int r_pvt, int r_sat, int r_eql, int r_fip);

      void set_density (spv_float density);
      void set_density_internal (const t_float *density);

      //! return prvd vector
      vval_vs_depth &get_prvd();
      //! return rsvd vector
      vval_vs_depth &get_rsvd();
      //! return pbvd vector
      vval_vs_depth &get_pbvd();

      spv_float
      get_rock()
      {
        return rock;
      }

      spv_float
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

      
      spv_int    get_i_array (const std::string &array_name);
      spv_float   get_fp_array (const std::string &array_name);

      bool contains_i_array (const std::string &array_name);
      bool contains_fp_array (const std::string &array_name);

      spv_int create_i_array (const std::string &name, t_int *dimens, t_int def_value = 0);
      spv_float create_fp_array (const std::string &name, t_int *dimens, t_float def_value = 0);

      int set_i_array (const std::string & array_name,  spv_int array);
      int set_fp_array (const std::string & array_name,  spv_float array);

    public:
      sp_prop_t props;

      
      

      /*

      int rpo_model;                   //!< 3-ph oil relative permeability model: flag 0, 1 or 2 (stone model)
      t_float minimal_pore_volume;      //!< Minimal pore volume allowed for active cells
      t_float minimal_splice_volume;    //!< Minimal pore volume allowed for active cells to splice with other cells
      t_float maximum_splice_thickness; //!< Default maximum thickness allowed between active cells to be coupled

      data_dimens dimens;              //!< dimension description

      t_int pvt_region;               //!< Number of PVT regions in simulation
      t_int sat_region;               //!< Number of saturation regions in simulation
      t_int eql_region;               //!< Number of equilibrium regions in simulation
      t_int fip_region;               //!< Number of FIP regions in simulation

      t_int fi_n_phase;               //!< number of phases if full implicit simulation
      t_int fi_phases;                //!< sizeof (int) bit fields (1 -- phase present, 0 -- do not present)

      t_int rock_region;              //!< Number of ROCK regions


      auto_value <int> init_section;             //!< flag indicating whether we have init section
      //!< if init_section == 0, cdata::check_sat will not go.

      //TITLE
      std::string title;


      */

      
      

      //! \brief class's public data area
      sp_h5_pool_t h5_pool;

      
      convert_units input_units_converter;          //!< Input units converter
      convert_units output_units_converter;         //!< Output units converter

      std::vector<rocktab_table> rocktab;    //!< rocktab tables for all rock regions

      
      pvt_vector pvto, pvtdo, pvtg, pvtw;

      spv_int equil_regions;

      spv_float rock;            //!< Array (pvt_region) - compressibility of rock for each region
      spv_float equil;
      spv_float p_ref;           //!< Array (pvt_region) - reference pressure for compressibility of rock for each region

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
}

#endif // BS_DATA_CLASS_H
