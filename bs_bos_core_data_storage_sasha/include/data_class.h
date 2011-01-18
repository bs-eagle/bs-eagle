#ifndef BS_DATA_CLASS_PN_H
#define BS_DATA_CLASS_PN_H
/*!
  \file data_class.h
  \brief initial data holder
	\author Nikonov Maxim
*/

#include "pool_treeish.h"
#include "ar_operat.h"
#include "rocktab_table.h"
#include "prvd_table.h"
#include "read_class.h"
#include "arrays.h"
//#include "arrays_tables.h"

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
template <class strategy_t>
class BS_API_PLUGIN idata : public bs_pool_node
{
private:
	struct idata_traits;

public:
	//! typedefs
	typedef idata<strategy_t>                           this_t;
	typedef smart_ptr<this_t , true>                    sp_this_t;

	typedef typename strategy_t::i_type_t                index_t;
	typedef typename strategy_t::fp_type_t                 item_t;
	//typedef typename strategy_t::index_array_t          index_array_t;
	//typedef typename strategy_t::item_array_t           item_array_t;

	typedef ar_args <strategy_t>                        ar_args_t;
	typedef ar_operat <strategy_t>                      ar_operat_t;

	typedef std::vector<val_vs_depth>                   vval_vs_depth;
	typedef std::vector< double >                       vec_d;
	typedef std::vector< int >                          vec_i;

      typedef smart_ptr <FRead, true>										  sp_reader_t;

	struct pvt_info
	{
		typedef std::vector<item_t> main_pvt_data_t;

		main_pvt_data_t							main_data_;
		auto_value <bool, false>		has_density_;
		auto_value <item_t>         density_;
		auto_value <item_t>         molar_density_;
	};

	typedef std::vector < pvt_info >											pvt_vector;

public:

	~idata ();

	//! \brief init idata members method
	void init();

	//! \brief Assignment operator
	this_t &operator=(const this_t&);

	void set_defaults_in_pool();
	void set_region (int r_pvt,int r_sat, int r_eql, int r_fip);
	void set_density (const std::vector <item_t> &density);
	void set_density_internal (const item_t *density);

	//! return prvd vector
	vval_vs_depth &get_prvd();
	//! return rsvd vector
	vval_vs_depth &get_rsvd();
	//! return pbvd vector
	vval_vs_depth &get_pbvd();

      std::vector<float> &get_rock()
	{
		return rock;
	}

      std::vector<float> &get_p_ref()
	{
		return p_ref;
	}

	int no_blanks (char *buf) const;
      void build_argument_list (const sp_reader_t &reader);
	void build_operator_list ();
	void output_argument_list ();
	void clear_argument_list ();
      void read_left (const sp_reader_t &reader,
			char *s, char **right, std::string &name, int *i1,
			int *i2, int *j1, int *j2, int *k1, int *k2, const char *keyword);
      void algorithm_read_and_done (const sp_reader_t &reader, char *buf, const char *keyword);
      void calculate (const sp_reader_t &reader, ar_args_t &l, char *right, const char *keyword);
	int read_numbers (ar_stack <ar_args_t> &st, char **buf);
      int read_operator (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op,
			char **start_ptr, int *flag_brek, const char *keyword,
			const char *buf, int *Arguments_count);
      int do_calculate (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op, const char *keyword);
      int prior_calculate (const sp_reader_t &reader, ar_args_t &res, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op, int priority, const char *keyword);
      int read_arg_func (const sp_reader_t &reader, ar_stack <ar_args_t> &ar, ar_stack <ar_operat_t> &op, char **start_ptr, const char *keyword, int *arguments_count, int flag_brek);
	int test_token (int prev, int cur);

	bs_array_i& get_int_array (const std::string& array_name);
	bs_array_fp& get_float_array (const std::string& array_name);

public:
	int rpo_model;                //!< 3-ph oil relative permeability model: flag 0, 1 or 2 (stone model)

	double minimal_pore_volume;     //!< Minimal pore volume allowed for active cells
	double minimal_splice_volume;   //!< Minimal pore volume allowed for active cells to splice with other cells
	double maximum_splice_thickness;//!< Default maximum thickness allowed between active cells to be coupled

	auto_value <int, false> restart;                  //!< flag show from which step start simulation
	auto_value <index_t> nx;                       //!< Number of nodes in X dimension (NX)
	auto_value <index_t> ny;                       //!< Number of nodes in Y dimension (NY)
	auto_value <index_t> nz;                       //!< Number of nodes in Z dimension (NZ)
	auto_value <index_t> pvt_region;               //!< Number of PVT regions in simulation
	auto_value <index_t> sat_region;               //!< Number of saturation regions in simulation
	auto_value <index_t> eql_region;               //!< Number of equilibrium regions in simulation
	auto_value <index_t> fip_region;               //!< Number of FIP regions in simulation

	auto_value <index_t> fi_n_phase;               //!< number of phases if full implicit simulation
	auto_value <int> fi_phases;                //!< sizeof (int) bit fields (1 -- phase present, 0 -- do not present)

	auto_value <index_t> rock_region;              //!< Number of ROCK regions

	vec_i equil_regions;

	//! \brief class's public data area
	//sp_pool_subnode d_pool() const;
	//sp_pool_subnode i_pool() const;

	auto_value <int> init_section;             //!< flag indicating whether we have init section
	//!< if init_section == 0, cdata::check_sat will not go.

	auto_value <int> units_in;                 //!< Number of units used for input data
	auto_value <int> units_out;                //!< Number of units used for output data

	//! \brief arguments
	std::map <std::string, ar_args_t> args;
	std::set <ar_operat_t> ops;

	int ar_tokens [2];

	std::vector<rocktab_table <strategy_t> > rocktab;    //!< rocktab tables for all rock regions

	//TITLE
	std::string title;

	pvt_vector pvto, pvtdo, pvtg, pvtw;

      std::vector < float > rock;            //!< Array (pvt_region) - compressibility of rock for each region
      std::vector < float > equil;
      std::vector < float > p_ref;           //!< Array (pvt_region) - reference pressure for compressibility of rock for each region

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

#endif // BS_DATA_CLASS_PN_H

