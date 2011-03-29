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
  permx_array = 0;
  permy_array = 0;
  permz_array = 0;
  multx_array = 0;
  multy_array = 0;
  multz_array = 0;
}



void
rs_smesh_base ::init_props (const sp_idata_t &idata)
{
  base_t::init_props (idata);
  
  nx = idata->props->get_i("nx");
  ny = idata->props->get_i("ny");
  nz = idata->props->get_i("nz");
  
  spv_float data_array;
  
  data_array = idata->get_fp_array("PERMX");
  if (data_array->size()) permx_array = &(*data_array)[0];
  
  data_array = idata->get_fp_array("PERMY");
  if (data_array->size()) permy_array = &(*data_array)[0];
  
  data_array = idata->get_fp_array("PERMZ");
  if (data_array->size()) permz_array = &(*data_array)[0];
  
  data_array = idata->get_fp_array("MULTX");
  if (data_array && data_array->size()) multx_array = &(*data_array)[0];
  
  data_array = idata->get_fp_array("MULTY");
  if (data_array && data_array->size()) multy_array = &(*data_array)[0];
  
  data_array = idata->get_fp_array("MULTZ");
  if (data_array && data_array->size()) multz_array = &(*data_array)[0];
}


void rs_smesh_base::inside_to_XYZ(const t_long index, t_long &i1,t_long &j1, t_long &k1) const
  {
    t_long r_index = (*base_t::base_t::int_to_ext)[index];
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
  
  if (nx <= 0)
    bs_throw_exception (boost::format ("nx = %d is out of range")% nx);
  if (ny <= 0)
    bs_throw_exception (boost::format ("ny = %d is out of range")% ny);
  if (nz <= 0)
    bs_throw_exception (boost::format ("nz = %d is out of range")% nz);
    
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
  if (d_dir == along_dim3)
    {
      elem_in_layers.resize(nz, 0);
      i_index = 0;
      for (int i = 0; i < nz; ++i)
        {
          for (int j = 0; j < nx*ny; ++j, ++i_index)
            {
              if (base_t::actnum_array[i_index])
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
                  if (base_t::actnum_array[i_index])
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
                  if (base_t::actnum_array[i_index])
                    ++elem_in_layers[i];
                }
            }
        }
    }
  return (int) elem_in_layers.size();
}

//BS_INST_STRAT(rs_smesh_base);
