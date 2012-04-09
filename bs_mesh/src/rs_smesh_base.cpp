/*!
	\file rs_smesh_base.cpp
	\brief This file implement interface class #rs_smesh_base which transfers  mesh data to the reservoir simulation process
	\author Iskhakov Ruslan
	\date 2008-05-20
 */
#include "bs_mesh_stdafx.h"

#include "rs_smesh_base.h"


rs_smesh_base ::rs_smesh_base ()
{
}



void
rs_smesh_base ::init_props (const sp_hdm_t hdm)
{
  base_t::init_props (hdm);
#ifndef PURE_MESH  
  
  nx = hdm->get_prop ()->get_i("nx");
  ny = hdm->get_prop ()->get_i("ny");
  nz = hdm->get_prop ()->get_i("nz");
  
  permx_array = hdm->get_pool ()->get_fp_data("PERMX");
  permy_array = hdm->get_pool ()->get_fp_data("PERMY");
  permz_array = hdm->get_pool ()->get_fp_data("PERMZ");
  multx_array = hdm->get_pool ()->get_fp_data("MULTX");
  multy_array = hdm->get_pool ()->get_fp_data("MULTY");
  multz_array = hdm->get_pool ()->get_fp_data("MULTZ");
#else
  nx = hdm->nx;
  ny = hdm->ny;
  nz = hdm->nz;
  
  permx_array = hdm->permx_array;
  permy_array = hdm->permy_array;
  permz_array = hdm->permz_array;
  multx_array = hdm->multx_array;
  multy_array = hdm->multy_array;
  multz_array = hdm->multz_array;
#endif
}


void rs_smesh_base::inside_to_XYZ(const t_long index, t_long &i1,t_long &j1, t_long &k1) const
  {
    
#ifndef PURE_MESH
    t_long r_index = (*base_t::base_t::int_to_ext)[index];
#else
    t_long r_index = base_t::base_t::int_to_ext[index];
#endif
    k1 = r_index / (nx*ny);
    j1 = (r_index - k1*nx*ny) / nx;
    i1 = r_index - k1*nx*ny - j1*nx;
  }

/*

int rs_smesh_base::init_int_to_ext()
{
  
  if (base_t::base_t::ext_to_int.size() == 0)
    return -1;
    
  base_t::base_t::int_to_ext.resize (base_t::base_t::n_active_elements, 0);

  for (size_t i = 0; i < base_t::base_t::ext_to_int.size(); i++)
    {
      if (base_t::base_t::ext_to_int[i] != -1)
        {
          base_t::base_t::int_to_ext [base_t::base_t::ext_to_int[i]] = (t_long)i;
        }
    }
  return 0;  
}
*/


void rs_smesh_base::check_data() const
{
  base_t::check_data ();

#ifndef PURE_MESH  
  if (nx <= 0)
    bs_throw_exception (boost::format ("nx = %d is out of range")% nx);
  if (ny <= 0)
    bs_throw_exception (boost::format ("ny = %d is out of range")% ny);
  if (nz <= 0)
    bs_throw_exception (boost::format ("nz = %d is out of range")% nz);
#endif

  // FIXME: init_props will raise exceptions if no array in hdm
  if (!permx_array)
    bs_throw_exception ("PERMX array is not initialized");
  if (!permy_array)
    bs_throw_exception ("PERMY array is not initialized");
  if (!permz_array)
    bs_throw_exception ("PERMZ array is not initialized");
}


int rs_smesh_base::get_elems_n_in_layers(const direction d_dir, stdv_int &elem_in_layers) const
{
  t_long i_index;
#ifndef PURE_MESH
  t_int const *actnum = base_t::actnum_array->data ();
#else
  t_int const *actnum = base_t::actnum_array;
#endif

  if (d_dir == along_dim3)
    {
      elem_in_layers.resize(nz, 0);
      i_index = 0;
      for (int i = 0; i < nz; ++i)
        {
          for (int j = 0; j < nx*ny; ++j, ++i_index)
            {
              if (actnum[i_index])
                ++elem_in_layers[i];
            }
        }
    }
  else if (d_dir == along_dim2)
    {
      elem_in_layers.resize(ny,0);
      for (int i = 0; i < ny; ++i)
        {
          for (int j = 0; j < nz; ++j)
            {
              i_index = nx*j+i*nx*nz;
              for (int k = 0; k < nx; ++k, ++i_index)
                {
                  if (actnum[i_index])
                    ++elem_in_layers[i];
                }
            }
        }
    }
  else //d_dir == along_dim1
    {
      elem_in_layers.resize(nx,0);
      for (int i = 0; i < nx; ++i)
        {
          for (int j = 0; j < nz; ++j)
            {
              i_index = ny*j+i*ny*nz;
              for (int k = 0; k < ny; ++k, ++i_index)
                {
                  if (actnum[i_index])
                    ++elem_in_layers[i];
                }
            }
        }
    }
  return (int) elem_in_layers.size();
}

//BS_INST_STRAT(rs_smesh_base);
