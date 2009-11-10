/*!
	\file mesh_rs.cpp
	\brief This file implement interface class #mesh_rs which transfers  mesh data to the reservoir simulation process
	\author Iskhakov Ruslan
	\date 2008-05-20
 */
#include "bs_mesh_stdafx.h"

#include "rs_smesh_base.h"

template <typename strategy_t>
void
mesh_rs <strategy_t>::init_props (const sp_idata_t &idata)
{
  base_t::init_props (idata);
  
  nx = idata->nx;
  ny = idata->ny;
  nz = idata->nz;
  
  sp_permx = idata->get_float_non_empty_array("PERMX");
  sp_permy = idata->get_float_non_empty_array("PERMY");
  sp_permz = idata->get_float_non_empty_array("PERMZ");
  sp_multx = idata->get_float_array("MULTX");
  sp_multy = idata->get_float_array("MULTY");
  sp_multz = idata->get_float_array("MULTZ");
}

template<class strategy_t>
void mesh_rs<strategy_t>::inside_to_XYZ(const index_t index, index_t &i1,index_t &j1, index_t &k1) const
  {
    index_t r_index = base_t::base_t::int_to_ext[index];
    k1 = r_index / (nx*ny);
    j1 = (r_index - k1*nx*ny) / nx;
    i1 = r_index - k1*nx*ny - j1*nx;
  }

template<class strategy_t>
int mesh_rs<strategy_t>::init_int_to_ext()
{
  if (base_t::base_t::ext_to_int.size() == 0)
    return -1;
    
  base_t::base_t::int_to_ext.resize (base_t::base_t::n_active_elements, 0);

  for (size_t i = 0; i < base_t::base_t::ext_to_int.size(); i++)
    {
      if (base_t::base_t::ext_to_int[i] != -1)
        {
          base_t::base_t::int_to_ext [base_t::base_t::ext_to_int[i]] = (index_t)i;
        }
    }
  return 0;  
}

template<class strategy_t>
void mesh_rs<strategy_t>::check_data() const
{
  base_t::check_data ();
  
  if (nx <= 0)
    bs_throw_exception (boost::format ("nx = %d is out of range")% nx);
  if (ny <= 0)
    bs_throw_exception (boost::format ("ny = %d is out of range")% ny);
  if (nz <= 0)
    bs_throw_exception (boost::format ("nz = %d is out of range")% nz);
    
  if (!sp_permx.size ())
    bs_throw_exception ("PERMX array is not initialized");
  if (!sp_permy.size ())
    bs_throw_exception ("PERMY array is not initialized");
  if (!sp_permz.size ())
    bs_throw_exception ("PERMZ array is not initialized");
}

template<class strategy_t>
int mesh_rs<strategy_t>::get_elems_n_in_layers(const direction d_dir, index_array_t &elem_in_layers) const
{
  index_t i_index;
  if (d_dir == along_dim3)
    {
      elem_in_layers.resize(nz, 0);
      i_index = 0;
      for (int i = 0; i < nz; ++i)
        {
          for (int j = 0; j < nx*ny; ++j, ++i_index)
            {
              if (base_t::sp_actnum[i_index])
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
                  if (base_t::sp_actnum[i_index])
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
                  if (base_t::sp_actnum[i_index])
                    ++elem_in_layers[i];
                }
            }
        }
    }
  return (int) elem_in_layers.size();
}

BS_INST_STRAT(mesh_rs);
