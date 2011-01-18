#include "bs_mesh_stdafx.h"

#include "mesh_ijk.h"
//#include "mesh_grdecl.h"
#include <stack>
#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bos_report.h"
#include "arrays.h"
#include "naive_file_reader.h"
#include BS_STOP_PLUGIN_IMPORT ()


using namespace std;
using namespace blue_sky;

#define BLOCK_NUM(i, j, k, nx, ny) (i)+(j)*(nx)+(k)*(nx)*(ny)

template <typename strategy_t>
mesh_ijk <strategy_t>::mesh_ijk ()
{
  dx_array   = 0;
  dy_array   = 0;
  dz_array   = 0;
  tops_array = 0;
}


template <typename strategy_t>
void
mesh_ijk <strategy_t>::init_props (const sp_idata_t &idata)
{
  sp_fp_storage_array_t data_array;
  
  data_array = idata->get_fp_non_empty_array("DX");
  if (data_array->size()) dx_array = &(*data_array)[0];
  
  data_array = idata->get_fp_non_empty_array("DY");
  if (data_array->size()) dy_array = &(*data_array)[0];
  
  data_array = idata->get_fp_non_empty_array("DZ");
  if (data_array->size()) dz_array = &(*data_array)[0];
  
  data_array = idata->get_fp_non_empty_array("TOPS");
  if (data_array->size()) tops_array = &(*data_array)[0];
  
  base_t::init_props (idata);
}

template<class strategy_t>
int mesh_ijk<strategy_t>::init_ext_to_int()
{
  i_type_t *ext_to_int_data, *int_to_ext_data;
  calc_shift_arrays();

  item_array_t volumes_temp (n_elements);
  int splicing_num = 0;//splicing(volumes_temp);

  //make proxy array
  ext_to_int->resize (n_elements);
  ext_to_int->assign(0);
  ext_to_int_data = &(*ext_to_int)[0];
  
  
  size_t n_count = 0;

  i_type_t nn_active = 0, i_index; //number of non-active previous cells
  for (i_type_t i = 0; i < nz; ++i)
    {
      for (i_type_t j = 0; j < ny; ++j)
        for (i_type_t k = 0; k < nx; ++k, ++n_count)
          {
            i_index = BLOCK_NUM (k, j, i, nx, ny);

            volumes_temp[i_index] = dx_array[i_index] * dy_array[i_index] * dz_array[i_index];

            if (!actnum_array[i_index])
              {
                nn_active++;
                ext_to_int_data[n_count] = -1;
              }
            else
              ext_to_int_data[n_count] = i_index - nn_active;
          }
    }
  init_int_to_ext();
  int_to_ext_data = &(*int_to_ext)[0];
  
  //fill volume array (except non-active block and using proxy array)
  volumes->resize(n_active_elements);
  fp_type_t *volumes_data = &(*volumes)[0];
  
  for (int i = 0; i < n_active_elements; ++i)
    volumes_data[i] = volumes_temp[int_to_ext_data[i]];

  calc_depths();
  return splicing_num;
}

template<class strategy_t>
void mesh_ijk<strategy_t>::check_data() const
{
  base_t::check_data ();

  if (!dx_array)
    bs_throw_exception ("DX array is not initialized");
  if (!dy_array)
    bs_throw_exception ("DY array is not initialized");
  if (!dz_array)
    bs_throw_exception ("DZ array is not initialized");
  if (!tops_array)
    bs_throw_exception ("TOPS array is not initialized");
}


#if 1
/*!
 * \brief  return coords of block vertexes
 * \param  i,j,k      - IJK index of block
 * \output cube_vertex - array of 3d-points - 8 block vertexes and block center coordinates
 */
template<class strategy_t>
grd_ecl::fpoint3d_vector
mesh_ijk<strategy_t>::calc_element (i_type_t index) const
  {
    grd_ecl::fpoint3d_vector cube_vertex;
    /*
     *                             X
     *                  0+-------+1
     *                  /|      /|
     *                 / |     / |
     *               2+--|----+3 |
     *              Y |  4+-- |--+5
     *                | /Z    | /
     *                |/      |/
     *             6  +-------+7
     */

    // upper
    cube_vertex[0] = fpoint3d (dx_shift_array[index],                dy_shift_array[index],                dz_shift_array[index]);
    cube_vertex[1] = fpoint3d (dx_shift_array[index] + dx_array[index], dy_shift_array[index],                dz_shift_array[index]);
    cube_vertex[2] = fpoint3d (dx_shift_array[index],                dy_shift_array[index] + dy_array[index], dz_shift_array[index]);
    cube_vertex[3] = fpoint3d (dx_shift_array[index] + dx_array[index], dy_shift_array[index] + dy_array[index], dz_shift_array[index]);
    // lower
    cube_vertex[4] = fpoint3d (dx_shift_array[index],                dy_shift_array[index],                dz_shift_array[index] + dz_array[index]);
    cube_vertex[5] = fpoint3d (dx_shift_array[index] + dx_array[index], dy_shift_array[index],                dz_shift_array[index] + dz_array[index]);
    cube_vertex[6] = fpoint3d (dx_shift_array[index],                dy_shift_array[index] + dy_array[index], dz_shift_array[index] + dz_array[index]);
    cube_vertex[7] = fpoint3d (dx_shift_array[index] + dx_array[index], dy_shift_array[index] + dy_array[index], dz_shift_array[index] + dz_array[index]);
    // center
    //cube_vertex[8] = fpoint3d (dx_shift_array[index] + dx_array[index]/2., dy_shift_array[index] + dy_array[index]/2., dz_shift_array[index] + dz_array[index]/2.));

    return cube_vertex;
  }

#endif

/*!
 * \brief  return coords of block vertexes
 * \param  index       - index of block in mesh
 * \output cube_vertex - array of 3d-points - 8 block vertexes and block center coordinates
 */
template<class strategy_t>
grd_ecl::fpoint3d_vector
mesh_ijk<strategy_t>::calc_element (const i_type_t i, const i_type_t j, const i_type_t k) const
  {
    i_type_t index = XYZ_to_inside (i, j, k);
    return calc_element (index);
  }

template <typename strategy_t>
typename mesh_ijk<strategy_t>::center_t
mesh_ijk<strategy_t>::get_center (i_type_t n_block) const
{
  BS_ASSERT (n_block != -1) (n_block);
  center_t center;
  
  center[0] = dx_shift_array[n_block] + dx_array[n_block] / 2;
  center[1] = dy_shift_array[n_block] + dy_array[n_block] / 2;
  center[2] = dz_shift_array[n_block] + dz_array[n_block] / 2;

  return center;
}  

template <typename strategy_t>
typename mesh_ijk<strategy_t>::center_t
mesh_ijk<strategy_t>::get_center (i_type_t i, i_type_t j, i_type_t k) const
{
  i_type_t n_block = BLOCK_NUM (i, j, k, nx, ny);
  return get_center (n_block);
}

template<class strategy_t>
int mesh_ijk<strategy_t>::splicing(item_array_t& volumes_temp)
{
  int splicing_num = 0;
  BS_ASSERT (false && "NOT IMPL YET");
  //for (i_type_t i = 0; i < nx; ++i)
  //{
  //  for (i_type_t j = 0; j < ny; ++j)
  //  {
  //    for (i_type_t k = 0; k < nz; ++k)
  //    {
  //      int index = i + j*nx + (k+1)*nx*ny;
  //      if ((*actnum_array)[index])
  //      {
  //        fp_type_t vol = dx_array[index]*dy_array[index]*dz_array[index];
  //        vol *= (*ntg_array)[index];
  //        if (vol*(*poro_array)[index] < minpv)
  //        {
  //          (*actnum_array)[index] = 0;
  //          splicing_num++;
  //        }
  //        volumes_temp[index] = vol;
  //      }
  //    }
  //  }
  //}
  return splicing_num;
}

template<class strategy_t>
int mesh_ijk<strategy_t>::build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn,
    sp_i_array_t boundary_array)
{
  n_connections = 0;
  sp_bcsr_t conn_trans;
  i_type_t n_non_zeros;
  i_type_t block_idx_ext, block_idx, next_block_idx_ext, conn_idx = 0;
  i_type_t i, j, k;

  jacobian->init_struct (n_active_elements, n_active_elements, n_active_elements);

  i_type_t *rows_ptr = &(*jacobian->get_rows_ptr())[0];
  i_type_t *ext_to_int_data = &(*ext_to_int)[0];
  
  for (i = 0; i < n_active_elements + 1; ++i)
    {
      rows_ptr[i] = 0;
    }

  boundary_array->clear();

  //all blocks are butting
  //first step - define and fill rows_ptr

  for (i = 0; i < nx; ++i)
    {
      for (j = 0; j < ny; ++j)
        {
          for (k = 0; k < nz; ++k)
            {
              block_idx_ext = BLOCK_NUM(i, j, k, nx, ny);
              if (!actnum_array[block_idx_ext])//skip non-active cells
                continue;

              //look only 3 positive-direction side
              next_block_idx_ext = BLOCK_NUM (i + 1, j, k, nx, ny);
              if ((i + 1 < nx) && actnum_array[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int_data[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int_data[next_block_idx_ext] + 1]++;
                }

              next_block_idx_ext = BLOCK_NUM (i, j + 1, k, nx, ny);
              if ((j + 1 < ny) && actnum_array[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int_data[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int_data[next_block_idx_ext] + 1]++;
                }

              next_block_idx_ext = BLOCK_NUM (i, j, k + 1, nx, ny);
              if ((k + 1 < nz) && actnum_array[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int_data[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int_data[next_block_idx_ext] + 1]++;
                }
            }
        }
    }

  //jacobian

  //sum rows_ptr
  for (i = 1; i < n_active_elements + 1; ++i)
    {
      rows_ptr[i]++; // diagonal
      rows_ptr[i] += rows_ptr[i - 1];
    }
  n_non_zeros = rows_ptr[n_active_elements];

  //create cols_ind
  jacobian->get_cols_ind()->resize (n_non_zeros);
  i_type_t *cols_ind = &(*jacobian->get_cols_ind())[0];

  //tran

  n_connections = (n_non_zeros - n_active_elements) / 2; //connection number

  if (n_connections == 0)
    {
      bs_throw_exception ("Mesh has no connections!");
    }

  conn_trans = flux_conn->get_conn_trans();
  conn_trans->init (n_connections, n_active_elements, 1, 2 * n_connections);

  i_type_t *rows_ptr_tran = &(*conn_trans->get_rows_ptr())[0];
  i_type_t *cols_ind_tran = &(*conn_trans->get_cols_ind())[0];
  fp_storage_type_t *values_tran =  &(*conn_trans->get_values())[0];
  
  flux_conn->get_matrix_block_idx_plus ()->resize (n_connections * 2);
  flux_conn->get_matrix_block_idx_plus ()->resize (n_connections * 2);
  
  i_type_t *matrix_block_idx_plus = &(*flux_conn->get_matrix_block_idx_plus ())[0];
  i_type_t *matrix_block_idx_minus = &(*flux_conn->get_matrix_block_idx_minus ())[0];

 

  for (i = 0; i < n_connections + 1; ++i)
    rows_ptr_tran[i] = i * 2;

  //additional array for current index in rows_ptr
  index_array_t tmp_rows_ptr;
  tmp_rows_ptr.resize (n_active_elements);

  // tmp_rows_ptr point to second value in each row
  // first value is always diagonal
  for (i = 0; i < n_active_elements; ++i)
    tmp_rows_ptr[i] = rows_ptr[i] + 1;

  for (i = 0; i < nx; ++i)
    {
      for (j = 0; j < ny; ++j)
        {
          for (k = 0; k < nz; ++k)
            {
              block_idx_ext = BLOCK_NUM (i, j, k, nx, ny);

              if (!actnum_array[block_idx_ext])//skip-non-active cells
                continue;

              block_idx = ext_to_int_data[block_idx_ext];

              //set diagonal value first in the row
              cols_ind[rows_ptr[block_idx]] = block_idx;

              //look only 3 positive-direction side active blocks
              next_block_idx_ext = BLOCK_NUM (i + 1, j, k, nx, ny);
              if ((i + 1 < nx) && actnum_array[next_block_idx_ext])
                {
                  set_neigbour_data (block_idx, block_idx_ext, next_block_idx_ext, conn_idx,
                                     rows_ptr, cols_ind, tmp_rows_ptr,
                                     matrix_block_idx_minus, matrix_block_idx_plus,
                                     cols_ind_tran, values_tran, along_dim1);
                }

              next_block_idx_ext = BLOCK_NUM (i, j + 1, k, nx, ny);
              if ((j+1 < ny) && actnum_array[next_block_idx_ext])//skip non-active
                {
                  set_neigbour_data (block_idx, block_idx_ext, next_block_idx_ext, conn_idx,
                                    rows_ptr, cols_ind, tmp_rows_ptr,
                                    matrix_block_idx_minus, matrix_block_idx_plus,
                                    cols_ind_tran, values_tran, along_dim2);
                }

              next_block_idx_ext = BLOCK_NUM (i, j, k + 1, nx, ny);
              if ((k+1 < nz) && actnum_array[next_block_idx_ext])//skip non-active
                {
                  set_neigbour_data (block_idx, block_idx_ext, next_block_idx_ext, conn_idx,
                                    rows_ptr, cols_ind, tmp_rows_ptr,
                                    matrix_block_idx_minus, matrix_block_idx_plus,
                                    cols_ind_tran, values_tran, along_dim3);
                }
            }
        }
    }
  return 0;
}

template<class strategy_t>
void mesh_ijk<strategy_t>::set_neigbour_data (const i_type_t index1, const i_type_t index1_ext, const i_type_t index2_ext, i_type_t &conn_idx,
                                              const i_type_t *rows_ptr, i_type_t *cols_ind, index_array_t &tmp_rows_ptr,
                                              i_type_t *m_memory, i_type_t *p_memory,
                                              i_type_t *cols_ind_tran, fp_storage_type_t *values_tran, direction dir)
{
  //change jacobian
  i_type_t index2 = (*ext_to_int)[index2_ext];

  cols_ind[tmp_rows_ptr[index1]] = index2;
  cols_ind[tmp_rows_ptr[index2]] = index1;

  m_memory[conn_idx] = rows_ptr[index1];
  m_memory[conn_idx + 1] = tmp_rows_ptr[index1];

  p_memory[conn_idx] = tmp_rows_ptr[index2];
  p_memory[conn_idx + 1] = rows_ptr[index2];

  tmp_rows_ptr[index1]++;
  tmp_rows_ptr[index2]++;

  //change flux_connection
  fp_type_t tran = calculate_tran(index1_ext, index2_ext, dir);

  cols_ind_tran[conn_idx] = index1;
  cols_ind_tran[conn_idx + 1] = index2;
  values_tran[conn_idx] = tran;
  values_tran[conn_idx + 1] = -tran;
  conn_idx += 2;
}

template<class strategy_t>
int mesh_ijk<strategy_t>::find_neighbours(sp_bcsr_t neig_matrix)
{
  return 0;
}

// calculating method have been taken from td eclipse (page 893)
template<class strategy_t>
typename mesh_ijk<strategy_t>::fp_type_t mesh_ijk<strategy_t>::calculate_tran(const i_type_t i, const i_type_t j, const  direction d_dir) const
  {
    fp_type_t tran;
    fp_type_t *depths_data = &(*depths)[0];
    i_type_t *ext_to_int_data = &(*ext_to_int)[0];

    fp_type_t A; //area between i and j block
    fp_type_t DIPC; //correction of inclination
    fp_type_t B,DHS,DVS; //additional variable

    fp_type_t ntg_i = 1;
    fp_type_t ntg_j = 1;
    if (!ntg_array)
      {
        ntg_i = ntg_array[i];
        ntg_j = ntg_array[j];
      }

    if (d_dir == along_dim1) //lengthwise OX
      {
        A = (dx_array[j] * dy_array[i] * dz_array[i] * ntg_i + dx_array[i] * dy_array[j] * dz_array[j] * ntg_j) / (dx_array[i] + dx_array[j]);
        B = (dx_array[i] / permx_array[i] + dx_array[j] / permx_array[j]) / 2;

        DHS = ((dx_array[i] + dx_array[j])/2) * ((dx_array[i] + dx_array[j])/2);
        DVS = (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]) * (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]);

        DIPC = DHS / (DHS + DVS);
        tran = darcy_constant * A * DIPC / B;
        if (multx_array)
          tran *= multx_array[i];
      }
    else if (d_dir == along_dim2) //lengthwise OY
      {
        A = (dy_array[j] * dx_array[i] * dz_array[i] * ntg_i + dy_array[i] * dx_array[j] * dz_array[j] * ntg_j) / (dy_array[i] + dy_array[j]);
        B = (dy_array[i] / permy_array[i] + dy_array[j] / permy_array[j]) / 2;

        DHS = ((dy_array[i] + dy_array[j]) / 2) * ((dy_array[i] + dy_array[j]) / 2);
        DVS = (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]) * (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]);

        DIPC = DHS / (DHS + DVS);
        tran = darcy_constant * A * DIPC / B;
        if (multy_array)
          tran *= multy_array[i];
      }
    else //lengthwise OZ
      {
        A = (dz_array[j] * dx_array[i] * dy_array[i] + dz_array[i] * dx_array[j] * dy_array[j]) / (dz_array[i] + dz_array[j]);
        B = (dz_array[i] / permz_array[i] + dz_array[j] / permz_array[j]) / 2;

        tran = darcy_constant * A / B;
        if (multz_array)
          tran *= multz_array[i];
      }

    return tran;
  }


template<class strategy_t>
void mesh_ijk<strategy_t>::get_block_dx_dy_dz(i_type_t n_elem, fp_type_t &dx, fp_type_t &dy, fp_type_t &dz) const
  {
    dx = dx_array[n_elem];
    dy = dy_array[n_elem];
    dz = dz_array[n_elem];
  }

template<class strategy_t>
int mesh_ijk<strategy_t>::calc_depths ()
{
  depths->resize (n_active_elements);
  fp_type_t *depths_data = &(*depths)[0];
  i_type_t index; //index of current block
  i_type_t *ext_to_int_data = &(*ext_to_int)[0];

  for (i_type_t i = 0; i < nx; ++i)
    for (i_type_t j = 0; j < ny; ++j)
      {
        fp_type_t current_tops = tops_array[i+j*nx]; //tops of current block
        for (i_type_t k = 0; k < nz; ++k)
          {
            index = i + j * nx + k * nx * ny;
            if (actnum_array[index])
              depths_data[ext_to_int_data[index]] = current_tops + dz_array[index]/2;
            current_tops += dz_array[index];
          }
      }
  return 0;
}

template<class strategy_t>
int mesh_ijk<strategy_t>::calc_shift_arrays()
{
  BS_ASSERT (nx * ny * nz) (nx) (ny) (nz);
  dx_shift_array.assign (nx * ny * nz, 0);
  dy_shift_array.assign (nx * ny * nz, 0);

  for (int k = 0; k < nz; ++k)
    {
      //dx_shift_array calculating
      for (int j = 0; j < ny; ++j)
        {
          size_t index = j*nx + k*nx*ny;
          dx_shift_array[index] =  0;
          index++;
          for (int i = 1; i < nx; ++i, ++index)
            dx_shift_array[index] = dx_shift_array[index-1] + dx_array[index];
        }
      //dy_shift_array calculating
      for (int i = 0; i < nx; ++i)
        {
          size_t index = i + k*nx*ny;
          dy_shift_array[index] = 0;
          index += nx;
          for (int j = 1; j < ny; ++j, index += nx)
            dy_shift_array[index] = dy_shift_array[index-nx] + dy_array[index];
        }
    }
  dz_shift_array.assign (nx * ny * nz, 0);
  int layer_size = nx*ny;
  for (int i = 0; i < nx; ++i)
    {
      for (int j = 0; j < ny; ++j)
        {
          size_t index = i+j*nx;
          dz_shift_array[index] = 0;
          index += layer_size;
          for (int k = 1; k < nz; ++k, index += layer_size)
            dz_shift_array[index] = dz_shift_array[index-layer_size] + dz_array[index];
        }
    }
  return 0;
}


BS_INST_STRAT(mesh_ijk);
