#include "bs_bos_core_data_storage_stdafx.h"

#include "rock_grid.h"
#include "arrays.h"
#include "data_class.h"

#include "plane_orientation.h"
#include "write_time_to_log.h"


#ifndef PI
#define PI      3.141592653589793238463
#endif

namespace blue_sky
  {
  // Special blue-sky objbase class implementations
  rock_grid::rock_grid (bs_type_ctor_param /*param*/)
      : bs_refcounter(), objbase()
  {
    n_elements = 0;
    n_pvt_regions = 0;
    //n_planes = 0;
    //n_boundary_planes = 0;
  }

  rock_grid::rock_grid (const rock_grid& prop)
      : bs_refcounter(), objbase(prop)
  {
    if (&prop != this)
      *this = prop;
  }

  rock_grid::~rock_grid ()
  {
  }

  // init rock_grid
  
  void rock_grid::init (const sp_idata &input_data, t_long n_els, t_int n_pvt_regs)
  {
    n_elements = n_els;
    n_pvt_regions = n_pvt_regs;

    net_to_gros.resize(n_elements);
    permeability.resize(n_elements * PLANE_ORIENTATION_TOTAL);
    porosity_p_ref.resize(n_elements);
    comp_const.resize(n_pvt_regions);
    comp_ref_pressure.resize(n_pvt_regions);

    BS_ASSERT (input_data->contains_fp_array ("MULTPV"));
    multpv.resize(n_elements);
  }

  /*!
  \brief convert  -- convert 'carray' to UfaSolver format and write it to 'array'

  \param msh    -- mesh
  \param array  -- pointer to destination array
  \param carray -- pointer to source array
  */

  void
    convert_permeability (t_long cells_count, const spv_long index_map, stdv_float perm, spv_float permx, spv_float permy, spv_float permz)
  {
    t_long *index_map_data = &(*index_map)[0];
    t_float *permx_data = &(*permx)[0];
    t_float *permy_data = &(*permy)[0];
    t_float *permz_data = &(*permz)[0];
    
    for (t_long i = 0; i < cells_count; ++i)
      {
        t_long ind = index_map_data [i];

        if (permx_data[ind] < 0)
          throw bs_exception ("convert_permeability", "the value of absolute permeability tensor element for node number [d] in X direction couldn't be less than 0");

        if (permy_data[ind] < 0)
          throw bs_exception ("convert_permeability", "the value of absolute permeability tensor element for node number [d] in Y direction couldn't be less than 0");

        if (permz_data[ind] < 0)
          throw bs_exception ("convert_permeability", "the value of absolute permeability tensor element for node number [d] in Z direction couldn't be less than 0");

        perm[FI_PLANE_ORI_X_IND (i)] = permx_data[ind];
        perm[FI_PLANE_ORI_Y_IND (i)] = permy_data[ind];
        perm[FI_PLANE_ORI_Z_IND (i)] = permz_data[ind];
      }
  }

  // init data
  
  int
  rock_grid::init_data (t_long cells_count, const spv_long index_map, const sp_idata &input_data)
  {
    if (input_data->get_rock()->size())
      {
        comp_const.assign (input_data->get_rock()->begin(), input_data->get_rock()->end());
      }
    else
      {
        bs_throw_exception ("Compressibility of rock has not been specified");
      }
    if (input_data->get_p_ref()->size())
      {
        comp_ref_pressure.assign (input_data->get_p_ref()->begin(), input_data->get_p_ref()->end());
      }
    else
      {
        // TODO: BUG: p_ref
        bs_throw_exception ("p_ref of rock has not been specified");
      }

    convert_arrays(cells_count, index_map, porosity_p_ref, input_data->get_fp_array("PORO"));
    convert_arrays (cells_count, index_map, net_to_gros, input_data->get_fp_array ("NTG"));
    convert_arrays (cells_count, index_map, multpv, input_data->get_fp_array ("MULTPV"));

    convert_permeability (cells_count, index_map, permeability, 
                          input_data->get_fp_array("PERMX"), 
                          input_data->get_fp_array("PERMY"), 
                          input_data->get_fp_array("PERMZ"));

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
  
  int
  rock_grid::init_planes_trans (t_long cells_count, const spv_float mesh_volumes, const sp_fi_params &/*ts_params*/, physical_constants &/*internal_constants*/)
  {
    // FIXME: net_to_gros don't used in this function
    if (net_to_gros.empty ())
      {
        bs_throw_exception ("net_to_gros is empty");
      }
    if (!mesh_volumes->size ())
      {
        bs_throw_exception ("mesh_volumes is empty");
      }

    // OPENMP
    volume.assign (mesh_volumes->begin (), mesh_volumes->end ());
    if (multpv.size () == volume.size ())
      {
        for (t_long i = 0; i < cells_count; ++i)
          {
            volume[i] *= multpv[i];
          }
      }
    else
      {
        bs_throw_exception (boost::format ("Size of multpv not equal with size of volume (multpv: %ld, volume: %ld)") 
                            % multpv.size () % volume.size ());
      }

    return 0;
  }


  BLUE_SKY_TYPE_STD_CREATE(rock_grid)
  BLUE_SKY_TYPE_STD_COPY(rock_grid)
  BLUE_SKY_TYPE_IMPL(rock_grid, objbase, "rock_grid", "rock_grid class", "rock grid")
}
