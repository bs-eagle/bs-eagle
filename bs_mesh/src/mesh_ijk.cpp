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
void
mesh_ijk <strategy_t>::init_props (const sp_idata_t &data)
{
  sp_dx   = data->get_float_non_empty_array("DX");
  sp_dy   = data->get_float_non_empty_array("DY");
  sp_dz   = data->get_float_non_empty_array("DZ");
  sp_tops = data->get_float_non_empty_array("TOPS");

  base_t::init_props (data);
}

template<class strategy_t>
int mesh_ijk<strategy_t>::init_ext_to_int()
{

  calc_shift_arrays();

  item_array_t volumes_temp (n_elements);
  int splicing_num = 0;//splicing(volumes_temp);

  //make proxy array
  ext_to_int.resize (n_elements, 0);
  size_t n_count = 0;

  index_t nn_active = 0, i_index; //number of non-active previous cells
  for (index_t i = 0; i < nz; ++i)
    {
      for (index_t j = 0; j < ny; ++j)
        for (index_t k = 0; k < nx; ++k, ++n_count)
          {
            i_index = BLOCK_NUM (k, j, i, nx, ny);

            volumes_temp[i_index] = sp_dx[i_index] * sp_dy[i_index] * sp_dz[i_index];

            if (!sp_actnum[i_index])
              {
                nn_active++;
                ext_to_int[n_count] = -1;
              }
            else
              ext_to_int[n_count] = i_index - nn_active;
          }
    }
  init_int_to_ext();

  //fill volume array (except non-active block and using proxy array)
  volumes.resize(n_active_elements);
  for (int i = 0; i < n_active_elements; ++i)
    volumes[i] = volumes_temp[int_to_ext[i]];

  calc_depths();
  return splicing_num;
}

template<class strategy_t>
void mesh_ijk<strategy_t>::check_data() const
{
  base_t::check_data ();

  if (!sp_dx.size ())
    bs_throw_exception ("DX array is not initialized");
  if (!sp_dy.size ())
    bs_throw_exception ("DY array is not initialized");
  if (!sp_dz.size ())
    bs_throw_exception ("DZ array is not initialized");
  if (!sp_tops.size ())
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
mesh_ijk<strategy_t>::calc_element (index_t index) const
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
    cube_vertex[1] = fpoint3d (dx_shift_array[index] + sp_dx[index], dy_shift_array[index],                dz_shift_array[index]);
    cube_vertex[2] = fpoint3d (dx_shift_array[index],                dy_shift_array[index] + sp_dy[index], dz_shift_array[index]);
    cube_vertex[3] = fpoint3d (dx_shift_array[index] + sp_dx[index], dy_shift_array[index] + sp_dy[index], dz_shift_array[index]);
    // lower
    cube_vertex[4] = fpoint3d (dx_shift_array[index],                dy_shift_array[index],                dz_shift_array[index] + sp_dz[index]);
    cube_vertex[5] = fpoint3d (dx_shift_array[index] + sp_dx[index], dy_shift_array[index],                dz_shift_array[index] + sp_dz[index]);
    cube_vertex[6] = fpoint3d (dx_shift_array[index],                dy_shift_array[index] + sp_dy[index], dz_shift_array[index] + sp_dz[index]);
    cube_vertex[7] = fpoint3d (dx_shift_array[index] + sp_dx[index], dy_shift_array[index] + sp_dy[index], dz_shift_array[index] + sp_dz[index]);
    // center
    //cube_vertex[8] = fpoint3d (dx_shift_array[index] + sp_dx[index]/2., dy_shift_array[index] + sp_dy[index]/2., dz_shift_array[index] + sp_dz[index]/2.));

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
mesh_ijk<strategy_t>::calc_element (const index_t i, const index_t j, const index_t k) const
  {
    index_t index = XYZ_to_inside (i, j, k);
    return calc_element (index);
  }

template <typename strategy_t>
typename mesh_ijk<strategy_t>::center_t
mesh_ijk<strategy_t>::get_center (index_t n_block) const
{
  BS_ASSERT (n_block != -1) (n_block);
  center_t center;
  
  center[0] = dx_shift_array[n_block] + sp_dx[n_block] / 2;
  center[1] = dy_shift_array[n_block] + sp_dy[n_block] / 2;
  center[2] = dz_shift_array[n_block] + sp_dz[n_block] / 2;

  return center;
}  

template <typename strategy_t>
typename mesh_ijk<strategy_t>::center_t
mesh_ijk<strategy_t>::get_center (index_t i, index_t j, index_t k) const
{
  index_t n_block = BLOCK_NUM (i, j, k, nx, ny);
  return get_center (n_block);
}

template<class strategy_t>
int mesh_ijk<strategy_t>::splicing(item_array_t& volumes_temp)
{
  int splicing_num = 0;
  BS_ASSERT (false && "NOT IMPL YET");
  //for (index_t i = 0; i < nx; ++i)
  //{
  //  for (index_t j = 0; j < ny; ++j)
  //  {
  //    for (index_t k = 0; k < nz; ++k)
  //    {
  //      int index = i + j*nx + (k+1)*nx*ny;
  //      if ((*sp_actnum)[index])
  //      {
  //        item_t vol = sp_dx[index]*sp_dy[index]*sp_dz[index];
  //        vol *= (*sp_ntg)[index];
  //        if (vol*(*sp_poro)[index] < minpv)
  //        {
  //          (*sp_actnum)[index] = 0;
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
int mesh_ijk<strategy_t>::build_jacobian_and_flux_connections (const sp_bcsr_t &jacobian, const sp_flux_conn_iface_t &flux_conn,
    index_array_t &boundary_array)
{
  n_connections = 0;
  sp_bcsr_t conn_trans;
  index_t n_non_zeros;
  index_t block_idx_ext, block_idx, next_block_idx_ext, conn_idx = 0;
  index_t i, j, k;

  jacobian->alloc_rows_ptr (n_active_elements);
  jacobian->n_cols = n_active_elements;

  index_array_t &rows_ptr = jacobian->get_rows_ptr();
  rows_ptr.assign (n_active_elements + 1, 0);

  boundary_array.clear();

  //all blocks are butting
  //first step - define and fill rows_ptr

  for (i = 0; i < nx; ++i)
    {
      for (j = 0; j < ny; ++j)
        {
          for (k = 0; k < nz; ++k)
            {
              block_idx_ext = BLOCK_NUM(i, j, k, nx, ny);
              if (!sp_actnum[block_idx_ext])//skip non-active cells
                continue;

              //look only 3 positive-direction side
              next_block_idx_ext = BLOCK_NUM (i + 1, j, k, nx, ny);
              if ((i + 1 < nx) && sp_actnum[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int[next_block_idx_ext] + 1]++;
                }

              next_block_idx_ext = BLOCK_NUM (i, j + 1, k, nx, ny);
              if ((j + 1 < ny) && sp_actnum[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int[next_block_idx_ext] + 1]++;
                }

              next_block_idx_ext = BLOCK_NUM (i, j, k + 1, nx, ny);
              if ((k + 1 < nz) && sp_actnum[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int[next_block_idx_ext] + 1]++;
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
  jacobian->alloc_cols_ind (n_non_zeros);
  index_array_t &cols_ind = jacobian->get_cols_ind();

  //tran

  n_connections = (n_non_zeros - n_active_elements) / 2; //connection number

  if (n_connections == 0)
    {
      bs_throw_exception ("Mesh has no connections!");
    }

  conn_trans = flux_conn->get_conn_trans();
  conn_trans->init (n_connections, n_active_elements, 1, 2 * n_connections);

  index_array_t &rows_ptr_tran = conn_trans->get_rows_ptr();
  index_array_t &cols_ind_tran = conn_trans->get_cols_ind();
  rhs_item_array_t &values_tran = conn_trans->get_values();
  index_array_t &matrix_block_idx_plus = flux_conn->get_matrix_block_idx_plus ();
  index_array_t &matrix_block_idx_minus = flux_conn->get_matrix_block_idx_minus ();

  matrix_block_idx_minus.resize (n_connections * 2, -1);
  matrix_block_idx_plus.resize (n_connections * 2, -1);

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

              if (!sp_actnum[block_idx_ext])//skip-non-active cells
                continue;

              block_idx = ext_to_int[block_idx_ext];

              //set diagonal value first in the row
              cols_ind[rows_ptr[block_idx]] = block_idx;

              //look only 3 positive-direction side active blocks
              next_block_idx_ext = BLOCK_NUM (i + 1, j, k, nx, ny);
              if ((i + 1 < nx) && sp_actnum[next_block_idx_ext])
                {
                  set_neigbour_data (block_idx, block_idx_ext, next_block_idx_ext, conn_idx,
                                     rows_ptr, cols_ind, tmp_rows_ptr,
                                     matrix_block_idx_minus, matrix_block_idx_plus,
                                     cols_ind_tran, values_tran, along_dim1);
                }

              next_block_idx_ext = BLOCK_NUM (i, j + 1, k, nx, ny);
              if ((j+1 < ny) && sp_actnum[next_block_idx_ext])//skip non-active
                {
                  set_neigbour_data (block_idx, block_idx_ext, next_block_idx_ext, conn_idx,
                                    rows_ptr, cols_ind, tmp_rows_ptr,
                                    matrix_block_idx_minus, matrix_block_idx_plus,
                                    cols_ind_tran, values_tran, along_dim2);
                }

              next_block_idx_ext = BLOCK_NUM (i, j, k + 1, nx, ny);
              if ((k+1 < nz) && sp_actnum[next_block_idx_ext])//skip non-active
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
void mesh_ijk<strategy_t>::set_neigbour_data (const index_t index1, const index_t index1_ext, const index_t index2_ext, index_t &conn_idx,
                                                   const index_array_t &rows_ptr, index_array_t &cols_ind, index_array_t& tmp_rows_ptr,
                                                   index_array_t& m_memory, index_array_t& p_memory,
                                                   index_array_t &cols_ind_tran, rhs_item_array_t &values_tran, direction dir)
{
  //change jacobian
  index_t index2 = ext_to_int[index2_ext];

  cols_ind[tmp_rows_ptr[index1]] = index2;
  cols_ind[tmp_rows_ptr[index2]] = index1;

  m_memory[conn_idx] = rows_ptr[index1];
  m_memory[conn_idx + 1] = tmp_rows_ptr[index1];

  p_memory[conn_idx] = tmp_rows_ptr[index2];
  p_memory[conn_idx + 1] = rows_ptr[index2];

  tmp_rows_ptr[index1]++;
  tmp_rows_ptr[index2]++;

  //change flux_connection
  item_t tran = calculate_tran(index1_ext, index2_ext, dir);

  cols_ind_tran[conn_idx] = index1;
  cols_ind_tran[conn_idx + 1] = index2;
  values_tran[conn_idx] = tran;
  values_tran[conn_idx + 1] = -tran;
  conn_idx += 2;
}

template<class strategy_t>
int mesh_ijk<strategy_t>::find_neighbours(sp_bcsr_t &neig_matrix)
{
  return 0;
}

// calculating method have been taken from td eclipse (page 893)
template<class strategy_t>
typename mesh_ijk<strategy_t>::item_t mesh_ijk<strategy_t>::calculate_tran(const index_t i, const index_t j, const  direction d_dir) const
  {
    item_t tran;

    item_t A; //area between i and j block
    item_t DIPC; //correction of inclination
    item_t B,DHS,DVS; //additional variable

    item_t ntg_i = 1;
    item_t ntg_j = 1;
    if (!sp_ntg.empty ())
      {
        ntg_i = sp_ntg[i];
        ntg_j = sp_ntg[j];
      }

    if (d_dir == along_dim1) //lengthwise OX
      {
        A = (sp_dx[j] * sp_dy[i] * sp_dz[i] * ntg_i + sp_dx[i] * sp_dy[j] * sp_dz[j] * ntg_j) / (sp_dx[i] + sp_dx[j]);
        B = (sp_dx[i] / sp_permx[i] + sp_dx[j] / sp_permx[j]) / 2;

        DHS = ((sp_dx[i] + sp_dx[j])/2) * ((sp_dx[i] + sp_dx[j])/2);
        DVS = (depths[ext_to_int[i]] - depths[ext_to_int[j]]) * (depths[ext_to_int[i]] - depths[ext_to_int[j]]);

        DIPC = DHS / (DHS + DVS);
        tran = darcy_constant * A * DIPC / B;
        if (!sp_multx.empty ())
          tran *= sp_multx[i];
      }
    else if (d_dir == along_dim2) //lengthwise OY
      {
        A = (sp_dy[j] * sp_dx[i] * sp_dz[i] * ntg_i + sp_dy[i] * sp_dx[j] * sp_dz[j] * ntg_j) / (sp_dy[i] + sp_dy[j]);
        B = (sp_dy[i] / sp_permy[i] + sp_dy[j] / sp_permy[j]) / 2;

        DHS = ((sp_dy[i] + sp_dy[j]) / 2) * ((sp_dy[i] + sp_dy[j]) / 2);
        DVS = (depths[ext_to_int[i]] - depths[ext_to_int[j]]) * (depths[ext_to_int[i]] - depths[ext_to_int[j]]);

        DIPC = DHS / (DHS + DVS);
        tran = darcy_constant * A * DIPC / B;
        if (!sp_multy.empty ())
          tran *= sp_multy[i];
      }
    else //lengthwise OZ
      {
        A = (sp_dz[j] * sp_dx[i] * sp_dy[i] + sp_dz[i] * sp_dx[j] * sp_dy[j]) / (sp_dz[i] + sp_dz[j]);
        B = (sp_dz[i] / sp_permz[i] + sp_dz[j] / sp_permz[j]) / 2;

        tran = darcy_constant * A / B;
        if (!sp_multz.empty ())
          tran *= sp_multz[i];
      }

    return tran;
  }


template<class strategy_t>
void mesh_ijk<strategy_t>::get_block_dx_dy_dz(index_t n_elem, item_t &dx, item_t &dy, item_t &dz) const
  {
    dx = sp_dx[n_elem];
    dy = sp_dy[n_elem];
    dz = sp_dz[n_elem];
  }

template<class strategy_t>
int mesh_ijk<strategy_t>::calc_depths ()
{
  depths.resize (n_active_elements,0);
  index_t index; //index of current block
  for (index_t i = 0; i < nx; ++i)
    for (index_t j = 0; j < ny; ++j)
      {
        item_t current_tops = sp_tops[i+j*nx]; //tops of current block
        for (index_t k = 0; k < nz; ++k)
          {
            index = i + j * nx + k * nx * ny;
            if (sp_actnum[index])
              depths[ext_to_int[index]] = current_tops + sp_dz[index]/2;
            current_tops += sp_dz[index];
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
            dx_shift_array[index] = dx_shift_array[index-1] + sp_dx[index];
        }
      //dy_shift_array calculating
      for (int i = 0; i < nx; ++i)
        {
          size_t index = i + k*nx*ny;
          dy_shift_array[index] = 0;
          index += nx;
          for (int j = 1; j < ny; ++j, index += nx)
            dy_shift_array[index] = dy_shift_array[index-nx] + sp_dy[index];
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
            dz_shift_array[index] = dz_shift_array[index-layer_size] + sp_dz[index];
        }
    }
  return 0;
}


BS_INST_STRAT(mesh_ijk);
