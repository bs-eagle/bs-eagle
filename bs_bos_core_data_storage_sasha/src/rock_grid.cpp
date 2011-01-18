#include "bs_bos_core_data_storage_stdafx.h"

#include "rock_grid.h"
#include "arrays.h"
#include "data_class.h"

#include "plane_orientation.h"


#ifndef PI
#define PI      3.141592653589793238463
#endif

namespace blue_sky {
// Special blue-sky objbase class implementations
template <class strategy_t/*, class mtrx_fp_type*/>
rock_grid<strategy_t /*,mtrx_fp_type*/>::rock_grid (bs_type_ctor_param /*param*/)
: bs_refcounter(), objbase()
{
	n_elements = 0;
	n_pvt_regions = 0;
	//n_planes = 0;
	//n_boundary_planes = 0;
}

template <class strategy_t/*, class mtrx_fp_type*/>
rock_grid<strategy_t /*,mtrx_fp_type*/>::rock_grid (const rock_grid& prop)
: bs_refcounter(), objbase(prop)
{
	if (&prop != this)
		*this = prop;
}

template <class strategy_t/*, class mtrx_fp_type*/>
rock_grid<strategy_t /*,mtrx_fp_type*/>::~rock_grid ()
{
}

// init rock_grid
	template <class strategy_t>
void rock_grid<strategy_t>::init (const sp_idata &input_data, index_t n_els, index_t n_pvt_regs)
{
	n_elements = n_els;
	n_pvt_regions = n_pvt_regs;

	net_to_gros.resize(n_elements);
	permeability.resize(n_elements * PLANE_ORIENTATION_TOTAL);
	porosity_p_ref.resize(n_elements);
	comp_const.resize(n_pvt_regions);
	comp_ref_pressure.resize(n_pvt_regions);

	if (input_data->d_pool()->contain("MULTPV"))
	{
		multpv.resize(n_elements);
	}
}

/*!
  \brief convert  -- convert 'carray' to UfaSolver format and write it to 'array'

  \param msh    -- mesh
  \param array  -- pointer to destination array
  \param carray -- pointer to source array
  */
	template <class index_array_t, class array_t, class src_t>
void convert_permeability (typename index_array_t::value_type cells_count, const index_array_t &index_map, array_t &perm, const src_t &permx, const src_t &permy, const src_t &permz)
{
	typedef typename index_array_t::value_type index_t;
	for (index_t i = 0; i < cells_count; ++i)
	{
		index_t ind = index_map [i];

		if (permx[ind] < 0)
			throw bs_exception ("convert_permeability", "the value of absolute permeability tensor element for node number [d] in X direction couldn't be less than 0");

		if (permy[ind] < 0)
			throw bs_exception ("convert_permeability", "the value of absolute permeability tensor element for node number [d] in Y direction couldn't be less than 0");

		if (permz[ind] < 0)
			throw bs_exception ("convert_permeability", "the value of absolute permeability tensor element for node number [d] in Z direction couldn't be less than 0");

		perm[FI_PLANE_ORI_X_IND (i)] = permx[ind];
		perm[FI_PLANE_ORI_Y_IND (i)] = permy[ind];
		perm[FI_PLANE_ORI_Z_IND (i)] = permz[ind];
	}
}

// init data
template <class strategy_t>
int rock_grid<strategy_t>::init_data (index_t cells_count, const index_array_t &index_map, const sp_idata &input_data)
{
	if (input_data->get_rock().size())
	{
		comp_const.assign (input_data->get_rock().begin(), input_data->get_rock().end());
	}
	else
	{
		bs_throw_exception ("Compressibility of rock has not been specified");
	}
	if (input_data->get_p_ref().size())
	{
		comp_ref_pressure.assign (input_data->get_p_ref().begin(), input_data->get_p_ref().end());
	}
	else
	{
		// TODO: BUG: p_ref
		bs_throw_exception ("p_ref of rock has not been specified");
	}

	if (!input_data->d_pool()->contain("PORO"))
		convert_arrays(cells_count, index_map, porosity_p_ref, input_data->carray_fp("PORO"));

	// initialize net to gross
	if (!input_data->d_pool()->contain("NTG"))
	{
		net_to_gros.assign(cells_count, 1.0);
	}
	else
	{
		convert_arrays (cells_count, index_map, net_to_gros, input_data->carray_fp("NTG"));
	}

	if (!input_data->d_pool()->contain ("PERMX") || !input_data->d_pool()->contain ("PERMY") || !input_data->d_pool()->contain ("PERMZ"))
	{
		bs_throw_exception ("Input data doesn't contain permx or permy or permz");
	}
	convert_permeability (cells_count, index_map, permeability, input_data->carray_fp("PERMX"), input_data->carray_fp("PERMY"), input_data->carray_fp("PERMZ"));

	if (input_data->d_pool()->contain("MULTPV"))
	{
		convert_arrays (cells_count, index_map, multpv, input_data->carray_fp("MULTPV"));
	}

	return 0;
}

/*!
 * \brief calculate geometric plane transmissibility (checked by Borschuk Oleg)
 *
 * \param msh                   -- mesh pointer
 * \param ts_params             -- pointer to ts_params
 * \param internal_constants    -- constants
 * \author Shabalin Maxim
 *
 * \return 0 if success
 */
template <class strategy_t>
int rock_grid<strategy_t>::init_planes_trans (index_t cells_count, const item_array_t &mesh_volumes, const sp_fi_params &/*ts_params*/, physical_constants &/*internal_constants*/)
{
	if (net_to_gros.empty ())
	{
		bs_throw_exception ("net_to_gros is empty");
	}
	if (mesh_volumes.empty ())
	{
		bs_throw_exception ("mesh_volumes is empty");
	}

	volume.assign (mesh_volumes.begin (), mesh_volumes.end ());
	if (!multpv.empty ())
	{
		if (multpv.size () == volume.size ())
		{
			for (index_t i = 0; i < cells_count; ++i)
			{
				volume[i] *= multpv[i];
			}
		}
		else
		{
			bs_throw_exception (boost::format ("Size of multpv not equal with size of volume (multpv: %ld, volume: %ld)") 
					% multpv.size () % volume.size ());
		}
	}

	return 0;
}

BLUE_SKY_TYPE_STD_CREATE_T_DEF(rock_grid,(class))
BLUE_SKY_TYPE_STD_COPY_T_DEF(rock_grid,(class))
BLUE_SKY_TYPE_IMPL_T_EXT(1, (rock_grid<base_strategy_fif>), 1, (objbase), "rock_grid_fi", "rock_grid_fi_i class", "rock grid fi i", false)
BLUE_SKY_TYPE_IMPL_T_EXT(1, (rock_grid<base_strategy_did>), 1, (objbase), "rock_grid_di", "rock_grid_di class", "rock grid di i", false)
BLUE_SKY_TYPE_IMPL_T_EXT(1, (rock_grid<base_strategy_dif>), 1, (objbase), "rock_grid_mixi", "rock_grid_di class", "rock grid mixi i", false)

}

