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
#include BS_STOP_PLUGIN_IMPORT ()


using namespace std;
using namespace blue_sky;

#define BLOCK_NUM(i, j, k, nx, ny) (i)+(j)*(nx)+(k)*(nx)*(ny)


mesh_ijk ::mesh_ijk ()
{
}



void
mesh_ijk ::init_props (const sp_hdm_t hdm)
{
  spv_float data_array;
  
  dx_array = hdm->get_pool ()->get_fp_data("DX");
  dy_array = hdm->get_pool ()->get_fp_data("DY");
  dz_array = hdm->get_pool ()->get_fp_data("DZ");
  tops_array = hdm->get_pool ()->get_fp_data("TOPS");
  
  base_t::init_props (hdm);
}


int mesh_ijk::init_ext_to_int()
{
  t_long *ext_to_int_data, *int_to_ext_data;
  calc_shift_arrays();

  // FIXME: splicing
  stdv_double volumes_temp (n_elements);
  int splicing_num = 0;//splicing(volumes_temp);

  //make proxy array
  ext_to_int->resize (n_elements);
  ext_to_int->assign(0);
  ext_to_int_data = ext_to_int->data ();
  
  
  size_t n_count = 0;

  t_long nn_active = 0, i_index; //number of non-active previous cells
  t_float const *dx = dx_array->data ();
  t_float const *dy = dy_array->data ();
  t_float const *dz = dz_array->data ();
  t_int const *actnum = actnum_array->data ();
  for (t_long i = 0; i < nz; ++i)
    {
      for (t_long j = 0; j < ny; ++j)
        for (t_long k = 0; k < nx; ++k, ++n_count)
          {
            i_index = BLOCK_NUM (k, j, i, nx, ny);

            volumes_temp[i_index] = dx[i_index] * dy[i_index] * dz[i_index];

            if (!actnum[i_index])
              {
                nn_active++;
                ext_to_int_data[n_count] = -1;
              }
            else
              ext_to_int_data[n_count] = i_index - nn_active;
          }
    }
  init_int_to_ext();
  int_to_ext_data = int_to_ext->data ();
  
  //fill volume array (except non-active block and using proxy array)
  volumes->resize(n_active_elements);
  t_double *volumes_data = &(*volumes)[0];
  
  for (int i = 0; i < n_active_elements; ++i)
    volumes_data[i] = volumes_temp[int_to_ext_data[i]];

  calc_depths();
  return splicing_num;
}


void mesh_ijk::check_data() const
{
  base_t::check_data ();
}


#if 1
/*!
 * \brief  return coords of block vertexes
 * \param  i,j,k      - IJK index of block
 * \output cube_vertex - array of 3d-points - 8 block vertexes and block center coordinates
 */

grd_ecl::fpoint3d_vector
mesh_ijk::calc_element (t_long index) const
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

    t_float const *dx = dx_array->data ();
    t_float const *dy = dy_array->data ();
    t_float const *dz = dz_array->data ();
    // upper
    cube_vertex[0] = fpoint3d (dx_shift_array[index],             dy_shift_array[index],             dz_shift_array[index]);
    cube_vertex[1] = fpoint3d (dx_shift_array[index] + dx[index], dy_shift_array[index],             dz_shift_array[index]);
    cube_vertex[2] = fpoint3d (dx_shift_array[index],             dy_shift_array[index] + dy[index], dz_shift_array[index]);
    cube_vertex[3] = fpoint3d (dx_shift_array[index] + dx[index], dy_shift_array[index] + dy[index], dz_shift_array[index]);
    // lower
    cube_vertex[4] = fpoint3d (dx_shift_array[index],             dy_shift_array[index],             dz_shift_array[index] + dz[index]);
    cube_vertex[5] = fpoint3d (dx_shift_array[index] + dx[index], dy_shift_array[index],             dz_shift_array[index] + dz[index]);
    cube_vertex[6] = fpoint3d (dx_shift_array[index],             dy_shift_array[index] + dy[index], dz_shift_array[index] + dz[index]);
    cube_vertex[7] = fpoint3d (dx_shift_array[index] + dx[index], dy_shift_array[index] + dy[index], dz_shift_array[index] + dz[index]);
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

grd_ecl::fpoint3d_vector
mesh_ijk::calc_element (const t_long i, const t_long j, const t_long k) const
  {
    t_long index = XYZ_to_inside (i, j, k);
    return calc_element (index);
  }


mesh_ijk::center_t
mesh_ijk::get_center (t_long n_block) const
{
  BS_ASSERT (n_block != -1) (n_block);
  center_t center;
  
  t_float const *dx = dx_array->data ();
  t_float const *dy = dy_array->data ();
  t_float const *dz = dz_array->data ();
  center[0] = dx_shift_array[n_block] + dx[n_block] / 2;
  center[1] = dy_shift_array[n_block] + dy[n_block] / 2;
  center[2] = dz_shift_array[n_block] + dz[n_block] / 2;

  return center;
}  


mesh_ijk::center_t
mesh_ijk::get_center (t_long i, t_long j, t_long k) const
{
  t_long n_block = BLOCK_NUM (i, j, k, nx, ny);
  return get_center (n_block);
}


int mesh_ijk::splicing(stdv_double& /*volumes_temp*/)
{
  int splicing_num = 0;
  BS_ASSERT (false && "NOT IMPL YET");
  //for (t_long i = 0; i < nx; ++i)
  //{
  //  for (t_long j = 0; j < ny; ++j)
  //  {
  //    for (t_long k = 0; k < nz; ++k)
  //    {
  //      int index = i + j*nx + (k+1)*nx*ny;
  //      if ((*actnum_array)[index])
  //      {
  //        t_double vol = dx_array[index]*dy_array[index]*dz_array[index];
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


int mesh_ijk::build_jacobian_and_flux_connections (const sp_bcsr_t jacobian, const sp_flux_conn_iface_t flux_conn,
    spv_long boundary_array)
{
  n_connections = 0;
  sp_bcsr_t conn_trans;
  t_long n_non_zeros;
  t_long block_idx_ext, block_idx, next_block_idx_ext, conn_idx = 0;
  t_long i, j, k;

  jacobian->init_struct (n_active_elements, n_active_elements, n_active_elements);

  t_long *rows_ptr = &(*jacobian->get_rows_ptr())[0];
  t_long *ext_to_int_data = &(*ext_to_int)[0];
  
  for (i = 0; i < n_active_elements + 1; ++i)
    {
      rows_ptr[i] = 0;
    }

  boundary_array->clear();

  //all blocks are butting
  //first step - define and fill rows_ptr

  t_int const *actnum = actnum_array->data ();
  for (i = 0; i < nx; ++i)
    {
      for (j = 0; j < ny; ++j)
        {
          for (k = 0; k < nz; ++k)
            {
              block_idx_ext = BLOCK_NUM(i, j, k, nx, ny);
              if (!actnum[block_idx_ext])//skip non-active cells
                continue;

              //look only 3 positive-direction side
              next_block_idx_ext = BLOCK_NUM (i + 1, j, k, nx, ny);
              if ((i + 1 < nx) && actnum[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int_data[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int_data[next_block_idx_ext] + 1]++;
                }

              next_block_idx_ext = BLOCK_NUM (i, j + 1, k, nx, ny);
              if ((j + 1 < ny) && actnum[next_block_idx_ext])
                {
                  rows_ptr[ext_to_int_data[block_idx_ext] + 1]++;
                  rows_ptr[ext_to_int_data[next_block_idx_ext] + 1]++;
                }

              next_block_idx_ext = BLOCK_NUM (i, j, k + 1, nx, ny);
              if ((k + 1 < nz) && actnum[next_block_idx_ext])
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
  t_long *cols_ind = &(*jacobian->get_cols_ind())[0];

  //tran

  n_connections = (n_non_zeros - n_active_elements) / 2; //connection number

  if (n_connections == 0)
    {
      bs_throw_exception ("Mesh has no connections!");
    }

  conn_trans = flux_conn->get_conn_trans();
  conn_trans->init (n_connections, n_active_elements, 1, 2 * n_connections);

  t_long *rows_ptr_tran = &(*conn_trans->get_rows_ptr())[0];
  t_long *cols_ind_tran = &(*conn_trans->get_cols_ind())[0];
  t_double *values_tran =  &(*conn_trans->get_values())[0];
  
  flux_conn->get_matrix_block_idx_plus ()->resize (n_connections * 2);
  flux_conn->get_matrix_block_idx_plus ()->resize (n_connections * 2);
  
  t_long *matrix_block_idx_plus = &(*flux_conn->get_matrix_block_idx_plus ())[0];
  t_long *matrix_block_idx_minus = &(*flux_conn->get_matrix_block_idx_minus ())[0];

 

  for (i = 0; i < n_connections + 1; ++i)
    rows_ptr_tran[i] = i * 2;

  //additional array for current index in rows_ptr
  stdv_long tmp_rows_ptr;
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

              if (!actnum[block_idx_ext])//skip-non-active cells
                continue;

              block_idx = ext_to_int_data[block_idx_ext];

              //set diagonal value first in the row
              cols_ind[rows_ptr[block_idx]] = block_idx;

              //look only 3 positive-direction side active blocks
              next_block_idx_ext = BLOCK_NUM (i + 1, j, k, nx, ny);
              if ((i + 1 < nx) && actnum[next_block_idx_ext])
                {
                  set_neigbour_data (block_idx, block_idx_ext, next_block_idx_ext, conn_idx,
                                     rows_ptr, cols_ind, tmp_rows_ptr,
                                     matrix_block_idx_minus, matrix_block_idx_plus,
                                     cols_ind_tran, values_tran, along_dim1);
                }

              next_block_idx_ext = BLOCK_NUM (i, j + 1, k, nx, ny);
              if ((j+1 < ny) && actnum[next_block_idx_ext])//skip non-active
                {
                  set_neigbour_data (block_idx, block_idx_ext, next_block_idx_ext, conn_idx,
                                    rows_ptr, cols_ind, tmp_rows_ptr,
                                    matrix_block_idx_minus, matrix_block_idx_plus,
                                    cols_ind_tran, values_tran, along_dim2);
                }

              next_block_idx_ext = BLOCK_NUM (i, j, k + 1, nx, ny);
              if ((k+1 < nz) && actnum[next_block_idx_ext])//skip non-active
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


void mesh_ijk::set_neigbour_data (const t_long index1, const t_long index1_ext, const t_long index2_ext, t_long &conn_idx,
                                              const t_long *rows_ptr, t_long *cols_ind, stdv_long &tmp_rows_ptr,
                                              t_long *m_memory, t_long *p_memory,
                                              t_long *cols_ind_tran, t_double *values_tran, direction dir)
{
  //change jacobian
  t_long index2 = (*ext_to_int)[index2_ext];

  cols_ind[tmp_rows_ptr[index1]] = index2;
  cols_ind[tmp_rows_ptr[index2]] = index1;

  m_memory[conn_idx] = rows_ptr[index1];
  m_memory[conn_idx + 1] = tmp_rows_ptr[index1];

  p_memory[conn_idx] = tmp_rows_ptr[index2];
  p_memory[conn_idx + 1] = rows_ptr[index2];

  tmp_rows_ptr[index1]++;
  tmp_rows_ptr[index2]++;

  //change flux_connection
  t_double tran = calculate_tran(index1_ext, index2_ext, dir);

  cols_ind_tran[conn_idx] = index1;
  cols_ind_tran[conn_idx + 1] = index2;
  values_tran[conn_idx] = tran;
  values_tran[conn_idx + 1] = -tran;
  conn_idx += 2;
}


int mesh_ijk::find_neighbours(sp_bcsr_t /*neig_matrix*/)
{
  return 0;
}

// calculating method have been taken from td eclipse (page 893)

t_double mesh_ijk::calculate_tran(const t_long i, const t_long j, const  direction d_dir) const
  {
    t_double tran;
    t_double *depths_data = depths->data ();
    t_long *ext_to_int_data = ext_to_int->data ();

    t_double A; //area between i and j block
    t_double DIPC; //correction of inclination
    t_double B,DHS,DVS; //additional variable

    t_float const *ntg = ntg_array->data ();
    t_float const *dx = dx_array->data ();
    t_float const *dy = dy_array->data ();
    t_float const *dz = dz_array->data ();
    t_float const *permx = permx_array->data ();
    t_float const *permy = permy_array->data ();
    t_float const *permz = permz_array->data ();
    t_float const *multx = multx_array->data ();
    t_float const *multy = multy_array->data ();
    t_float const *multz = multz_array->data ();

    t_double ntg_i = 1;
    t_double ntg_j = 1;
    if (!ntg_array)
      {
        ntg_i = ntg[i];
        ntg_j = ntg[j];
      }

    if (d_dir == along_dim1) //lengthwise OX
      {
        A = (dx[j] * dy[i] * dz[i] * ntg_i + dx[i] * dy[j] * dz[j] * ntg_j) / (dx[i] + dx[j]);
        B = (dx[i] / permx[i] + dx[j] / permx[j]) / 2;

        DHS = ((dx[i] + dx[j])/2) * ((dx[i] + dx[j])/2);
        DVS = (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]) * (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]);

        DIPC = DHS / (DHS + DVS);
        tran = darcy_constant * A * DIPC / B;
        if (multx_array)
          tran *= multx[i];
      }
    else if (d_dir == along_dim2) //lengthwise OY
      {
        A = (dy[j] * dx[i] * dz[i] * ntg_i + dy[i] * dx[j] * dz[j] * ntg_j) / (dy[i] + dy[j]);
        B = (dy[i] / permy[i] + dy[j] / permy[j]) / 2;

        DHS = ((dy[i] + dy[j]) / 2) * ((dy[i] + dy[j]) / 2);
        DVS = (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]) * (depths_data[ext_to_int_data[i]] - depths_data[ext_to_int_data[j]]);

        DIPC = DHS / (DHS + DVS);
        tran = darcy_constant * A * DIPC / B;
        if (multy_array)
          tran *= multy[i];
      }
    else //lengthwise OZ
      {
        A = (dz[j] * dx[i] * dy[i] + dz[i] * dx[j] * dy[j]) / (dz[i] + dz[j]);
        B = (dz[i] / permz[i] + dz[j] / permz[j]) / 2;

        tran = darcy_constant * A / B;
        if (multz_array)
          tran *= multz[i];
      }

    return tran;
  }



void mesh_ijk::get_block_dx_dy_dz(t_long n_elem, t_double &dx, t_double &dy, t_double &dz) const
  {
    dx = dx_array->data ()[n_elem];
    dy = dy_array->data ()[n_elem];
    dz = dz_array->data ()[n_elem];
  }


int mesh_ijk::calc_depths ()
{
  depths->resize (n_active_elements);
  t_double *depths_data = &(*depths)[0];
  t_long index; //index of current block
  t_long *ext_to_int_data = &(*ext_to_int)[0];

  t_float const *tops = tops_array->data ();
  t_int const *actnum = actnum_array->data ();
  t_float const *dz = dz_array->data ();
  for (t_long i = 0; i < nx; ++i)
    for (t_long j = 0; j < ny; ++j)
      {
        t_float current_tops = tops[i+j*nx]; //tops of current block
        for (t_long k = 0; k < nz; ++k)
          {
            index = i + j * nx + k * nx * ny;
            if (actnum[index])
              depths_data[ext_to_int_data[index]] = current_tops + dz[index]/2;
            current_tops += dz[index];
          }
      }
  return 0;
}


int mesh_ijk::calc_shift_arrays()
{
  BS_ASSERT (nx * ny * nz) (nx) (ny) (nz);
  dx_shift_array.assign (nx * ny * nz, 0);
  dy_shift_array.assign (nx * ny * nz, 0);

  t_float const *dx = dx_array->data ();
  t_float const *dy = dy_array->data ();
  t_float const *dz = dz_array->data ();

  for (int k = 0; k < nz; ++k)
    {
      //dx_shift_array calculating
      for (int j = 0; j < ny; ++j)
        {
          size_t index = j*nx + k*nx*ny;
          dx_shift_array[index] =  0;
          index++;
          for (int i = 1; i < nx; ++i, ++index)
            dx_shift_array[index] = dx_shift_array[index-1] + dx[index];
        }
      //dy_shift_array calculating
      for (int i = 0; i < nx; ++i)
        {
          size_t index = i + k*nx*ny;
          dy_shift_array[index] = 0;
          index += nx;
          for (int j = 1; j < ny; ++j, index += nx)
            dy_shift_array[index] = dy_shift_array[index-nx] + dy[index];
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
            dz_shift_array[index] = dz_shift_array[index-layer_size] + dz[index];
        }
    }
  return 0;
}


//BS_INST_STRAT(mesh_ijk);
