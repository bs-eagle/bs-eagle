/*! \file mesh_grdecl.cpp
\brief This file implement class for working with grid_ecllipse
\author Iskhakov Ruslan
\date 2008-05-20 */
#include "bs_mesh_stdafx.h"
#include "mesh_grdecl.h"

using namespace grd_ecl;
using namespace blue_sky;

const char filename_hdf5[] = "grid_swap.h5";

#ifdef BS_MESH_WRITE_TRANSMISS_MATRIX  
  FILE*  fp;
#endif BS_MESH_WRITE_TRANSMISS_MATRIX 


template<class strategy_t>
inline void
mesh_grdecl<strategy_t>::calc_corner_point(const pool_item_t z, const pool_item_t *coord, fpoint3d_t &p)const
  {
    p.z = z;
    /*
    float temp = (p.z-m_Coord.pStart.z)/(m_Coord.pEnd.z-m_Coord.pStart.z);
    p.x = temp *(m_Coord.pEnd.x-m_Coord.pStart.x)+m_Coord.pStart.x;
    p.y = temp *(m_Coord.pEnd.y-m_Coord.pStart.y)+m_Coord.pStart.y;
    */
    float temp = (p.z - coord[2]) / (coord[5] - coord[2]);
    p.x = temp *(coord[3] - coord[0]) + coord[0];
    p.y = temp *(coord[4] - coord[1]) + coord[1];
  }
  
//! get element corners index in zcorn_array of block[i,j,k]
 /* nodes layout
     *                             X
     *                    0+-------+1
     *                    /|     / |
     *                  /  |   /   |
     *               2/-------+3   |
     *              Y |   4+--|----+5
     *                |   /Z  |   /
     *                | /     | /
     *              6 /-------/7
     */

template<class strategy_t>
inline void
typename mesh_grdecl<strategy_t>::get_element_zcorn_index(index_t i, index_t j, index_t k, element_zcorn_index_t& element)  const
{
  //typename mesh_grdecl<strategy_t>::element_zcorn_index_t element;
  
  index_t index1 = i * 2 + j * 4 * nx + k * 8 * nx * ny;
  index_t index2 = index1 + 4 * nx * ny;

  element[0] = index1;
  element[1] = index1 + 1;
  element[2] = index1 + 2 * nx;
  element[3] = index1 + 2 * nx + 1;

  element[4] = index2;
  element[5] = index2 + 1;
  element[6] = index2 + 2 * nx;
  element[7] = index2 + 2 * nx + 1;
}

//! get element corners index in zcorn_array of block[i,j,k]
template<class strategy_t>
typename mesh_grdecl<strategy_t>::plane_zcorn_index_t
typename mesh_grdecl<strategy_t>::get_plane_zcorn_index (element_zcorn_index_t element, 
                                                         element_plane_orientation_t orientation)  const
{
  typename mesh_grdecl<strategy_t>::plane_zcorn_index_t plane;
  switch (orientation)
    {
      case x_axis_minus:  //left_cross
        plane[0] = element[0];
        plane[1] = element[2];
        plane[2] = element[4];
        plane[3] = element[6];
        break;
      case x_axis_plus:   //right_cross
        plane[0] = element[1];
        plane[1] = element[3];
        plane[2] = element[5];
        plane[3] = element[7];
        break;
      case y_axis_minus:  //top_cross
        plane[0] = element[0];
        plane[1] = element[1];
        plane[2] = element[4];
        plane[3] = element[5];
        break;
      case y_axis_plus:   //bottom_cross
        plane[0] = element[2];
        plane[1] = element[3];
        plane[2] = element[6];
        plane[3] = element[7];
        break;
      case z_axis_minus:  //lower_cross
        plane[0] = element[0];
        plane[1] = element[2];
        plane[2] = element[1];
        plane[3] = element[3];
        break;
      case z_axis_plus:   //upper_cross
        plane[0] = element[4];
        plane[1] = element[6];
        plane[2] = element[5];
        plane[3] = element[7];
        break;
      
      default:
        bs_throw_exception ("Invalid orientation!");;
    }
}

template<class strategy_t>
typename mesh_grdecl<strategy_t>::element_t
mesh_grdecl<strategy_t>::calc_element (const index_t index) const
  {
    index_t i, j, k;
    element_t element;
    
    inside_to_XYZ (index, i, j, k);
    calc_element (i, j, k, element);
    return element;
  }

template<class strategy_t>
void
mesh_grdecl<strategy_t>::calc_element (const index_t i, const index_t j, const index_t k, element_t &element) const
  {
    corners_t corners;

#ifdef _DEBUG
    BS_ASSERT ( (i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz));
#endif

    /* nodes layout
     *                             X
     *                    0+-------+1
     *                    /|     / |
     *                  /  |   /   |
     *               2/-------+3   |
     *              Y |   4+--|----+5
     *                |   /Z  |   /
     *                | /     | /
     *              6 /-------/7
     */

    //define index

    index_t index1 = i * 2 + j * 4 * nx + k * 8 * nx * ny;//upper side
    index_t index2 = index1 + 4 * nx * ny;//lower side
    index_t iCOORD = i + j * (nx + 1);

    calc_corner_point (zcorn_array[index1], &coord_array[iCOORD * 6], corners[0]);
    calc_corner_point (zcorn_array[index1 + 1], &coord_array[(iCOORD + 1) * 6], corners[1]);
    calc_corner_point (zcorn_array[index1 + 2 * nx], &coord_array[(iCOORD + (nx + 1)) * 6], corners[2]);
    calc_corner_point (zcorn_array[index1 + 2 * nx + 1], &coord_array[(iCOORD + (nx + 1) + 1) * 6], corners[3]);

    calc_corner_point (zcorn_array[index2], &coord_array[iCOORD * 6], corners[4]);
    calc_corner_point (zcorn_array[index2 + 1], &coord_array[(iCOORD + 1) * 6], corners[5]);
    calc_corner_point (zcorn_array[index2 + 2 * nx], &coord_array[(iCOORD + (nx + 1)) * 6], corners[6]);
    calc_corner_point (zcorn_array[index2 + 2 * nx + 1], &coord_array[(iCOORD + (nx + 1) + 1) * 6], corners[7]);
    
    element.init (corners);
  }
  
template<class strategy_t>
bool mesh_grdecl<strategy_t>::is_small(const index_t i, const index_t j, const index_t k, item_t eps)  const
  {
    if (k >= nz)
      return false;

    item_t dz1, dz2, dz3, dz4; //height for each coord
    //define index
    index_t index1 = i*2+j*4*nx+k*8*nx*ny;	//lower side
    index_t index2 = index1+4*nx*ny;			//upper side
    dz1 = zcorn_array[index2]-zcorn_array[index1];
    dz2 = zcorn_array[index2+1]-zcorn_array[index1+1];
    dz3 = zcorn_array[index2+2*nx]-zcorn_array[index1+2*nx];
    dz4 = zcorn_array[index2+2*nx+1]-zcorn_array[index1+2*nx+1];

    if (dz1 <= eps && dz2 <= eps && dz3 <= eps && dz4 <= eps)
      return true;
    return false;
  }

template<class strategy_t>
void
mesh_grdecl<strategy_t>::get_plane_crossing_projection_on_all_axis(const plane_t &plane1, const plane_t &plane2, fpoint3d_t &A)const
  {
    quadrangle_t quad1, quad2;

    //project on YoZ
    for (size_t i = 0; i < N_PLANE_CORNERS; i++)
      {
        quad1[i] = fpoint2d(plane1[i].y, plane1[i].z);
        quad2[i] = fpoint2d(plane2[i].y, plane2[i].z);
      }

    A.x = get_quadrangle_crossing_area (quad1, quad2, EPS_SQ);

    //project on XoZ
    for (size_t i = 0; i < N_PLANE_CORNERS; i++)
      {
        quad1[i] = fpoint2d(plane1[i].x, plane1[i].z);
        quad2[i] = fpoint2d(plane2[i].x, plane2[i].z);
      }
    A.y = get_quadrangle_crossing_area (quad1, quad2, EPS_SQ);

    if ((A.x + A.y )!= 0.0)
      {
        //project on XoY
        for (size_t i = 0; i < N_PLANE_CORNERS; i++)
          {
            quad1[i] = fpoint2d (plane1[i].x, plane1[i].y);
            quad2[i] = fpoint2d (plane2[i].x, plane2[i].y);
          }
        A.z = get_quadrangle_crossing_area (quad1, quad2, EPS_SQ);
     }
  }

template<class strategy_t>
typename mesh_grdecl<strategy_t>::item_t 
mesh_grdecl<strategy_t>::get_center_zcorn(const element_zcorn_index_t &element)const
  {
    item_t res = 0.0;
    size_t i, n = element.size();
    
    for (i = 0; i < n; i++)
      res += zcorn_array[element[i]];
      
    return res/n;
  }



template<class strategy_t>
int mesh_grdecl<strategy_t>::init_ext_to_int()
{
 
  write_time_to_log init_time ("Mesh initialization", ""); 
  item_array_t volumes_temp(n_elements);

  //tools::save_seq_vector ("actnum.bs.txt").save sp_actnum;
  
  int splicing_num = splicing(volumes_temp);
  
  //check_adjacency (1);
  //tools::save_seq_vector ("active_blocks.bs.txt").save sp_actnum;
  
  check_adjacency ();
  
  //make proxy array
  ext_to_int.resize(n_elements,0);
  size_t n_count = 0;

  index_t nn_active = 0, i_index; //number of non-active previous cells
  for (index_t i = 0; i < nz; ++i)
    {
      for (index_t j = 0; j < ny; ++j)
        for (index_t k = 0; k < nx; ++k, ++n_count)
          {
            i_index = k + (nx * j) + (i * nx * ny);
            
            if (!sp_actnum[i_index])
              {
                nn_active++;
                ext_to_int[n_count] = -1;
              }
            else
              ext_to_int[n_count] = i_index-nn_active;
          }
    }

  //tools::save_seq_vector ("ext_to_int.bs.txt").save (ext_to_int);

  init_int_to_ext();

  //fill volume array (except non-active block and using proxy array)
  volumes.resize(n_active_elements);
  for (int i = 0; i < n_active_elements; ++i)
    volumes[i] = volumes_temp[int_to_ext[i]];
    
  calc_depths();
  

  //bs_throw_exception ("NOT IMPL YET");
  return n_active_elements;
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::splice_two_blocks (const index_t i, const index_t j, const index_t k, const index_t k1)
{
  index_t index, index1, index2;

  index1 = i * 2 + j * 4 * nx + k1 * 8 * nx * ny; //upper side of [i,j,k1]
  index2 = i * 2 + j * 4 * nx + (k1 * 8 + 4) * nx * ny; //lower side of [i,j,k1]
  if (k > k1)
    {
      index = i * 2 + j * 4 * nx + k * 8 * nx * ny; //upper side of [i,j,k]

      zcorn_array[index] = zcorn_array[index1];
      zcorn_array[index + 1] = zcorn_array[index1 + 1];
      zcorn_array[index + 2 * nx] = zcorn_array[index1 + 2 * nx];
      zcorn_array[index + 2 * nx + 1] = zcorn_array[index1 + 2 * nx + 1];

      zcorn_array[index2] = zcorn_array[index1];
      zcorn_array[index2 + 1] = zcorn_array[index1 + 1];
      zcorn_array[index2 + 2 * nx] = zcorn_array[index1 + 2 * nx];
      zcorn_array[index2 + 2 * nx + 1] = zcorn_array[index1 + 2 * nx + 1];
    }
  else
    {
      index = i * 2 + j * 4 * nx + (k * 8 + 4) * nx * ny; //lower side of [i,j,k]

      zcorn_array[index] = zcorn_array[index2];
      zcorn_array[index + 1] = zcorn_array[index2 + 1];
      zcorn_array[index + 2 * nx] = zcorn_array[index2 + 2 * nx];
      zcorn_array[index + 2 * nx + 1] = zcorn_array[index2 + 2 * nx + 1];

      zcorn_array[index1] = zcorn_array[index2];
      zcorn_array[index1 + 1] = zcorn_array[index2 + 1];
      zcorn_array[index1 + 2 * nx] = zcorn_array[index2 + 2 * nx];
      zcorn_array[index1 + 2 * nx + 1] = zcorn_array[index2 + 2 * nx + 1];
    }
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::are_two_blocks_close (const index_t i, const index_t j, const index_t k, const index_t k1)
{
  index_t index, index1;
  if (k > k1)
    {
      index = i * 2 + j * 4 * nx + k * 8 *nx * ny; //upper side of [i,j,k]
      index1 = i * 2 + j * 4 * nx + (k1 * 8 + 4) * nx * ny; //lower side of [i,j,k1]
    }
  else
    {
      index = i * 2 + j * 4 * nx + (k * 8 + 4) * nx * ny; //lower side of [i,j,k]
      index1 = i * 2 + j * 4 * nx + k1 * 8 * nx * ny; //upper side of [i,j,k1]
    }

  if ((fabs (zcorn_array[index1] - zcorn_array[index]) < max_thickness) &&
      (fabs (zcorn_array[index1 + 1] - zcorn_array[index + 1]) < max_thickness) &&
      (fabs (zcorn_array[index1 + 2 * nx] - zcorn_array[index + 2 * nx]) < max_thickness) &&
      (fabs (zcorn_array[index1 + 2 * nx + 1] - zcorn_array[index + 2 * nx + 1]) < max_thickness))
    {
      return true;
    }
  else
    {
      return false;
    }
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::check_adjacency(int shift_zcorn)
{
   index_t i, j, k;
   index_t index, zindex, zindex1;
   index_t n_adjacent = 0;
   array_uint8_t &actnum     = sp_actnum;
   
  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j)
      for (k = 0; k < nz; ++k)
        {
          
          index = i + j * nx + k * nx * ny;
          
          
          // miss inactive blocks
          if (!actnum[index])
            {
              continue;
            }
                       
          // check blocks adjacency
          zindex = i * 2 + j * 4 * nx + k * 8 *nx * ny;
          zindex1 = zindex + 2; // check next by x
         
          if (i + 1 == nx ||
              zcorn_array[zindex + 1] == zcorn_array[zindex1] &&
              zcorn_array[zindex +  2 * nx + 1] == zcorn_array[zindex1 +  2 * nx] &&
              zcorn_array[zindex + 4 * nx * ny + 1] == zcorn_array[zindex1 + 4 * nx * ny ] &&
              zcorn_array[zindex + 4 * nx * ny + 2 * nx + 1] == zcorn_array[zindex1 + 4 * nx * ny + 2 * nx])
            {
              
              zindex1 = zindex + 4 * nx; // check next by y
              if (j + 1 == ny ||
                  zcorn_array[zindex + 2 * nx] == zcorn_array[zindex1] &&
                  zcorn_array[zindex + 2 * nx + 1] == zcorn_array[zindex1 +  1] &&
                  zcorn_array[zindex + 4 * nx * ny + 2 * nx] == zcorn_array[zindex1 + 4 * nx * ny] &&
                  zcorn_array[zindex + 4 * nx * ny + 2 * nx + 1] == zcorn_array[zindex1+ 4 * nx * ny + 1])
                  {
                    n_adjacent++;
                  }
            }
            if (shift_zcorn && (i + j) % 2 == 1)
              {
                zcorn_array[zindex] += 2;
                zcorn_array[zindex + 1] += 2;
                zcorn_array[zindex +  2 * nx] += 2;
                zcorn_array[zindex +  2 * nx + 1] += 2;
                zcorn_array[zindex + 4 * nx * ny] += 2;
                zcorn_array[zindex + 4 * nx * ny + 1] += 2;
                zcorn_array[zindex + 4 * nx * ny +  2 * nx] += 2;
                zcorn_array[zindex + 4 * nx * ny +  2 * nx + 1] += 2;
              }
       }
  
  BOSOUT (section::mesh, level::medium) << "  adjacent cells:"<< n_adjacent <<" ("<< n_adjacent * 100 / (n_active_elements)<<"% active)" << bs_end;
  return (n_adjacent == n_active_elements);
}

template<class strategy_t>
int mesh_grdecl<strategy_t>::splicing(item_array_t& volumes_temp)
{
  index_t nCount = 0;
  item_t vol_sum, vol_block;
  index_t i, j, k, k1;
  index_t small_block_top, big_block_top;
  index_t n_inactive_orig, n_inactive_vol, n_incative_splice;
  index_t index;
  element_t element;

#ifdef _DEBUG
    BS_ASSERT (zcorn_array.size() && coord_array.size());
#endif

  array_uint8_t &actnum     = sp_actnum;
  const array_float16_t &poro = sp_poro;
  

  n_inactive_orig = 0;
  n_inactive_vol = 0;
  n_incative_splice = 0;

  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j)
      {
        small_block_top = -1;
        big_block_top = -1;
        vol_sum = 0.0;
        for (k = 0; k < nz; ++k)
          {
            index = i + j * nx + k * nx * ny;
            
            // miss inactive blocks
            if (!actnum[index])
              {
                // blocks can`t be spliced across inactive blocks
                big_block_top = -1;
                small_block_top = -1;
                vol_sum = 0.0;
                n_inactive_orig++;
                continue;
              }
            calc_element (i, j, k, element);
            vol_block = element.calc_volume ();
            item_t vol_block_poro = vol_block * poro[index];

            if (vol_block_poro <= minpv)
              {
                // block is too small, set as inactive
                actnum[index] = 0;
                ++nCount;
                n_inactive_vol++;
                // blocks can`t be spliced across inactive blocks
                vol_sum = 0.0;
                big_block_top = -1;
                small_block_top = -1;
                continue;
              }
            else if (vol_block_poro < minsv)
              {
                // this is small block

                volumes_temp[index] = vol_block;

                if (big_block_top != -1)
                  {
                    // this block is absorbed by bigger block above
                    splice_two_blocks (i, j, big_block_top, k);
                    // add volume this small block to the volume of big block
                    volumes_temp[i + j * nx + big_block_top * nx * ny] += vol_block;

                    // make block inactive
                    actnum[index] = 0;
                    ++nCount;
                    n_incative_splice++;
                    small_block_top = -1;
                  }

                // check if this block is close enough to next block
                if ((k < (nz - 1)) && (are_two_blocks_close (i, j, k, k + 1)))
                  {
                    if (big_block_top == -1)
                      {
                        if (small_block_top == -1)
                          {
                            // if already not spliced, can be spliced with lower block
                            small_block_top = k;
                          }

                        // aggregate volume of small blocks in case they`ll be spliced with lower block
                        vol_sum += vol_block;
                      }
                  }
                else
                  {
                    // this block and next block are not close, can`t be coupled
                    small_block_top = -1;
                    big_block_top = -1;
                    vol_sum = 0.0;
                    // TODO: make multperm = 0
                  }

              }
            else
              {
                // this is big block

                volumes_temp[index] = vol_block;

                if (small_block_top != -1)
                  {
                    // this block absorbes all smaller blocks above
                    for (k1 = k - 1; k1 >= small_block_top; --k1)
                      {
                        splice_two_blocks (i, j, k, k1);
                        n_incative_splice++;
                        // make small block inactive
                        actnum[i + j * nx + k1 * nx * ny] = 0;
                        ++nCount;
                      }

                    // add volume all small blocks above to the volume of this big block
                    volumes_temp[index] += vol_sum;
                    vol_sum = 0.0;
                    small_block_top = -1;
                  }

                // check if this block is close enough to next block
                if ((k < (nz - 1)) && (are_two_blocks_close (i, j, k, k + 1)))
                  {
                    // this block can absorb lower small blocks
                    big_block_top = k;
                  }
                else
                  {
                    // this block and next block are not close, can`t be coupled
                    small_block_top = -1;
                    big_block_top = -1;
                    // TODO: make multperm = 0
                  }
              }
           }
        }
  
  
  index_t n_total = n_active_elements + n_inactive_orig;
  /*
  if (n_total != nx * ny * nz)
  
    {
      BOSOUT (section::mesh, level::error) << "MESH_GRDECL total cells assert failed!" << bs_end;
      return -1;
    }  
  */      
  
  BOSOUT (section::mesh, level::medium) << "Mesh cells info:" << bs_end;
  BOSOUT (section::mesh, level::medium) << "  total: "<< n_total << bs_end; 
  BOSOUT (section::mesh, level::medium) << "  initial active: "<< n_active_elements <<" ("<< n_active_elements * 100 / (n_total)<<"%)" << bs_end;
  BOSOUT (section::mesh, level::medium) << "  marked inactive: "<< nCount << " (" << n_inactive_vol << " by volume, " << n_incative_splice << " by splice)" << bs_end;
  n_active_elements -= nCount;
  BOSOUT (section::mesh, level::medium) << "  total active: "<< n_active_elements <<" ("<< n_active_elements * 100 / (n_total)<<"%)" << bs_end;
  
  return nCount;
}


template<class strategy_t>
int mesh_grdecl<strategy_t >::calc_depths ()
{
  depths.resize(n_active_elements,0);
  index_t index = 0;
  element_zcorn_index_t element;
  for (index_t k = 0; k < nz; ++k)
    for (index_t j = 0; j < ny; ++j)
      for (index_t i = 0; i < nx; ++i, ++index)
        {
          if (sp_actnum[index])
            {
              get_element_zcorn_index(i, j, k, element);
              depths[ext_to_int[index]] = get_center_zcorn(element);
            }
        }
  return 0;
}


static int n_tran_calc = 0;

// calculating method have been taken from td eclipse (page 896)
// calc transmissibility between fully adjacent cells index1 and index2
template<class strategy_t>
typename mesh_grdecl<strategy_t>::item_t 
mesh_grdecl<strategy_t>::calc_tran(const index_t ext_index1, const index_t ext_index2, const plane_t &plane1, 
                                       const fpoint3d_t &center1, const fpoint3d_t &center2, direction d_dir, plane_t* plane2) const
{
  item_t tran;
  fpoint3d_t D1, D2, A, plane_contact_center;
  
  n_tran_calc ++;
  
  // if (plane2 != 0) then column is not adjacent - but there still can be adjacent cells, so check it
  if (plane2 == 0 || ((plane1[0].z == (*plane2)[0].z) && (plane1[1].z == (*plane2)[1].z) && (plane1[2].z == (*plane2)[2].z) && (plane1[3].z == (*plane2)[3].z)))
  //if (plane2 == 0)
    {
      get_plane_center(plane1, plane_contact_center);
      plane_contact_center.distance_to_point (center1, D1);
      plane_contact_center.distance_to_point (center2, D2);
      A = get_projection_on_all_axis_for_one_side(plane1);
    }
  else 
    {
      fpoint3d_t plane1_center, plane2_center;
      
      get_plane_center(plane1, plane1_center);
      get_plane_center(*plane2, plane2_center);
      
      plane1_center.distance_to_point (center1, D1);
      plane2_center.distance_to_point (center2, D2);
      
      // check simple cases of plane intersection  
      
      // 1. plane1 is completely inside plane2
      if ((plane1[0].z >= (*plane2)[0].z) && (plane1[1].z >= (*plane2)[1].z) && (plane1[2].z <= (*plane2)[2].z) && (plane1[3].z <= (*plane2)[3].z))
        A = get_projection_on_all_axis_for_one_side(plane1);
      else 
      // 2. plane2 is completely inside plane1  
      if ((plane1[0].z <= (*plane2)[0].z) && (plane1[1].z <= (*plane2)[1].z) && (plane1[2].z >= (*plane2)[2].z) && (plane1[3].z >= (*plane2)[3].z))
        A = get_projection_on_all_axis_for_one_side(*plane2);
      else
      // 3. plane1 is lower than plane2
      if ((plane1[0].z >= (*plane2)[0].z) && (plane1[1].z >= (*plane2)[1].z) && (plane1[0].z <= (*plane2)[2].z) && (plane1[1].z <= (*plane2)[3].z) && (plane1[2].z >= (*plane2)[2].z) && (plane1[3].z >= (*plane2)[3].z))
        {
          (*plane2)[0] = plane1[0];
          (*plane2)[1] = plane1[1];
          A = get_projection_on_all_axis_for_one_side(*plane2);
        }
      else    
       // 4. plane2 is lower than plane1
      if ((plane1[0].z <= (*plane2)[0].z) && (plane1[1].z <= (*plane2)[1].z) && (plane1[2].z >= (*plane2)[0].z) && (plane1[3].z >= (*plane2)[1].z) && (plane1[2].z <= (*plane2)[2].z) && (plane1[3].z <= (*plane2)[3].z))
        {
          (*plane2)[2] = plane1[2];
          (*plane2)[3] = plane1[3];
          A = get_projection_on_all_axis_for_one_side(*plane2);
        }
      else    
        get_plane_crossing_projection_on_all_axis(plane1, *plane2, A);
    }
  


  float koef1, koef2 ; //koef = (A,Di)/(Di,Di)
  
  koef1 = A*D1 / (D1*D1);
  koef2 = A*D2 / (D2*D2);
  
  #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index1 < 1000)
          {
             /*
            fprintf (fp, "[A: %lf; %lf; %lf K1: %lf, K2: %lf]\n", A.x, A.y, A.z, koef1, koef2);
           
            fprintf (fp, "[D1: %lf; %lf; %lf, D2: %lf; %lf; %lf, plane_center: %lf; %lf; %lf, center1: %lf; %lf; %lf, center2: %lf; %lf; %lf]", \
            D1.x, D1.y, D1.z, D2.x, D2.y, D2.z,  plane_contact_center.x, plane_contact_center.y, plane_contact_center.z,  
            center1.x, center1.y, center1.z, center2.x, center2.y, center2.z);
            fprintf (fp, "[D1: %lf; %lf; %lf, D2: %lf; %lf; %lf, plane_center: %lf; %lf; %lf, center1: %lf; %lf; %lf, center2: %lf; %lf; %lf]", \
            D1.x, D1.y, D1.z, D2.x, D2.y, D2.z,  plane_contact_center.x, plane_contact_center.y, plane_contact_center.z,  
            center1.x, center1.y, center1.z, center2.x, center2.y, center2.z);
            */
          }
  #endif //BS_MESH_WRITE_TRANSMISS_MATRIX   

  if (koef1 < 10e-16 || koef2 < 10e-16)
    {
      /*
      BOSWARN (section::mesh, level::warning) 
        << boost::format ("For indexes (%d, %d) transmissibility will be set to 0 because koef1 = 0 (%f) or koef2 = 0 (%f)")
        % ext_index1 % ext_index2 % koef1 % koef2 
        << bs_end;
      */
      return 0;
    }

  item_t Ti, Tj;

  item_t ntg_index1 = 1;
  item_t ntg_index2 = 1;
  if (!sp_ntg.empty ())
    {
      ntg_index1 = sp_ntg[ext_index1];
      ntg_index2 = sp_ntg[ext_index2];
    }

  if (d_dir == along_dim1) //lengthwise OX
    {
      Ti = sp_permx[ext_index1]*ntg_index1*koef1;
      Tj = sp_permx[ext_index2]*ntg_index2*koef2;
      tran = darcy_constant / (1 / Ti + 1 / Tj);
      if (!sp_multx.empty ())
        tran *= sp_multx[ext_index1];
    }
  else if (d_dir == along_dim2) //lengthwise OY
    {
      Ti = sp_permy[ext_index1]*ntg_index1*koef1;
      Tj = sp_permy[ext_index2]*ntg_index2*koef2;
      tran = darcy_constant / (1 / Ti + 1 / Tj);
      if (!sp_multy.empty ())
        tran *= sp_multy[ext_index1];
    }
  else //lengthwise OZ
    {
      Ti = sp_permz[ext_index1]*koef1;
      Tj = sp_permz[ext_index2]*koef2;
      tran = darcy_constant / (1 / Ti + 1 / Tj);
      if (!sp_multz.empty ())
        tran *= sp_multz[ext_index1];
    }
  
  /*
  #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
    if (ext_index1 < 1000)
      fprintf (fp, "[Ti: %lf; Tj: %lf]", Ti, Tj);
  #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
  */
  return tran;
}



template<class strategy_t>
int mesh_grdecl<strategy_t>::build_jacobian_and_flux_connections (const sp_bcsr_t &jacobian,const sp_flux_conn_iface_t & flux_conn,
                                                                  index_array_t &boundary_array)

{
  return build_jacobian_and_flux_connections_add_boundary (jacobian, flux_conn, boundary_array);
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::init_props(const sp_idata_t &idata)
{
  base_t::init_props (idata);

  // init ZCORN
  zcorn_array = idata->get_float_non_empty_array("ZCORN");
  min_z = *(std::min_element(zcorn_array.begin(),zcorn_array.end()));
  max_z = *(std::max_element(zcorn_array.begin(),zcorn_array.end()));

  // init COORD
  coord_array = idata->get_float_non_empty_array("COORD");

  max_x = min_x = coord_array[0];
  max_y = min_y = coord_array[1];

  array_float16_t::iterator it, e = coord_array.end ();

  for (it = coord_array.begin (); it != e; it += 6)
    {
      // move matching points apart
      if (it[2] == it[5])
        it[5] += 1.0f;

      //looking for max&min coordinates
      if (min_x > it[0]) min_x = it[0];
      if (min_x > it[3]) min_x = it[3];

      if (min_y > it[1]) min_y = it[1];
      if (min_y > it[4]) min_y = it[4];

      if (max_x < it[0]) max_x = it[0];
      if (max_x < it[3]) max_x = it[3];

      if (max_y < it[1]) max_y = it[1];
      if (max_y < it[4]) max_y = it[4];

    }

}

template<class strategy_t>
void mesh_grdecl<strategy_t>::check_data() const
{
  base_t::check_data ();

  if (min_x < 0)
    bs_throw_exception (boost::format ("min_x = %d is out of range")% min_x);
  if (min_y < 0)
    bs_throw_exception (boost::format ("min_y = %d is out of range")% min_y);
  if (min_z < 0)
    bs_throw_exception (boost::format ("min_z = %d is out of range")% min_z);

  if (!coord_array.size ())
    bs_throw_exception ("COORD array is not initialized");
  if (!zcorn_array.size ())
    bs_throw_exception ("ZCORN array is not initialized");
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::get_block_dx_dy_dz (index_t n_elem, item_t &dx, item_t &dy, item_t &dz) const
  {
    element_t &elem = calc_element(n_elem);
    point3d_t sizes = elem.get_dx_dy_dz (); 
    dx = sizes[0];
    dy = sizes[1];
    dz = sizes[2];
  }

template<class strategy_t>
float mesh_grdecl<strategy_t>:: get_block_dx(index_t n_elem) const
  {
    return calc_element(n_elem).get_dx();
  }

template<class strategy_t>
float mesh_grdecl<strategy_t>:: get_block_dy(index_t n_elem) const
  {
    return calc_element(n_elem).get_dy();
  }

template<class strategy_t>
float mesh_grdecl<strategy_t>:: get_block_dz(index_t n_elem) const
  {
    return calc_element(n_elem).get_dz();
  }

template<class strategy_t>
float  mesh_grdecl<strategy_t>::get_dtop(index_t n_elem) const
{
  element_t elem;
  
  elem = calc_element (n_elem);

  return elem.get_center().z - elem.get_dz();
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::generate_array()
{
#if 0
  index_t n_size = n_elements;
  sp_poro->clear();
  sp_ntg->clear();

  sp_permx->clear();
  sp_permy->clear();
  sp_permz->clear();

  sp_multx->clear();
  sp_multy->clear();
  sp_multz->clear();


  for (index_t i =0; i < n_size; ++i)
    {
      sp_poro->push_back(0.2f);
      sp_ntg->push_back(0.4f);

      sp_permx->push_back(10.0f);
      sp_permy->push_back(10.0f);
      sp_permz->push_back(10.0f);

      sp_multx->push_back(1.0f);
      sp_multy->push_back(1.0f);
      sp_multz->push_back(1.0f);
    }
#endif
}


template <typename strategy_t, typename loop_t>
struct build_jacobian_rows_class
{
  typedef typename strategy_t::index_t          index_t;
  typedef typename strategy_t::item_t           item_t;
  typedef typename strategy_t::index_array_t    index_array_t;
  typedef typename strategy_t::item_array_t     item_array_t;

  typedef mesh_grdecl <strategy_t>                mesh_t;
  typedef typename mesh_t::plane_t                plane_t;
  typedef typename mesh_t::element_zcorn_index_t  element_zcorn_index_t;

  build_jacobian_rows_class (mesh_grdecl <strategy_t> *mesh, loop_t *loop, std::set <index_t, std::less <index_t> > &boundary_set, index_array_t *rows_ptr)
  : mesh (mesh)
  , loop (loop)
  , boundary_set (boundary_set)
  , rows_ptr (rows_ptr)
  , nx (mesh->nx)
  , ny (mesh->ny)
  , nz (mesh->nz)
  {
  }

  void
  prepare (index_t i, index_t j, index_t k)
  {
    ext_index  = i + j * nx + k * nx * ny;
  }

  // check if (i, j) column of cells is adjacent to neigbours
  // that is true, if every cell of a column is fully adjacent to X and Y neighbour
  
  bool
  check_column_adjacency (index_t i, index_t j)
  {
    bool flag = true;
    index_t k, zindex, zindex1, index;
         
    for (k = 0; k < nz; ++k)
      {
        index = i + j * nx + k * nx * ny;
         
        // miss inactive blocks
        if (!mesh->sp_actnum[index])
          {
            continue;
          }
          
        // check blocks adjacency
        zindex = i * 2 + j * 4 * nx + k * 8 *nx * ny;
        zindex1 = zindex + 2; // check next by x
       
        if (i + 1 == nx ||
            mesh->zcorn_array[zindex + 1] == mesh->zcorn_array[zindex1] &&
            mesh->zcorn_array[zindex +  2 * nx + 1] == mesh->zcorn_array[zindex1 +  2 * nx] &&
            mesh->zcorn_array[zindex + 4 * nx * ny + 1] == mesh->zcorn_array[zindex1 + 4 * nx * ny ] &&
            mesh->zcorn_array[zindex + 4 * nx * ny + 2 * nx + 1] == mesh->zcorn_array[zindex1 + 4 * nx * ny + 2 * nx])
          {
            
            zindex1 = zindex + 4 * nx; // check next by y
            if (j + 1 == ny ||
                mesh->zcorn_array[zindex + 2 * nx] == mesh->zcorn_array[zindex1] &&
                mesh->zcorn_array[zindex + 2 * nx + 1] == mesh->zcorn_array[zindex1 +  1] &&
                mesh->zcorn_array[zindex + 4 * nx * ny + 2 * nx] == mesh->zcorn_array[zindex1 + 4 * nx * ny] &&
                mesh->zcorn_array[zindex + 4 * nx * ny + 2 * nx + 1] == mesh->zcorn_array[zindex1+ 4 * nx * ny + 1])
                {
                  // cell is adjacent;
                }
             else
              {
                flag = false;
                break;
              }
          }
         else
          {
            flag = false;
            break;
          }  
       }
    loop->is_column_adjacent[j * nx + i] = flag;
    return flag;
  }

  void
  change_by_x (index_t i, index_t j, index_t k, index_t ext_index2, bool is_adjacent)
  {
    (*rows_ptr)[mesh->convert_ext_to_int (ext_index) + 1]++;
    (*rows_ptr)[mesh->convert_ext_to_int (ext_index2) + 1]++;
  }

  void
  change_by_y (index_t i, index_t j, index_t k, index_t ext_index2, bool is_adjacent)
  {
    (*rows_ptr)[mesh->convert_ext_to_int (ext_index) + 1]++;
    (*rows_ptr)[mesh->convert_ext_to_int (ext_index2) + 1]++;
  }

  void
  change_by_z (index_t i, index_t j, index_t k, index_t ext_index2, bool is_adjacent)
  {
    (*rows_ptr)[mesh->convert_ext_to_int (ext_index) + 1]++;
    (*rows_ptr)[mesh->convert_ext_to_int (ext_index2) + 1]++;
  }

  void
  add_boundary (index_t external_cell_index)
  {
    boundary_set.insert (external_cell_index);
  }

  mesh_grdecl <strategy_t>                  *mesh;
  loop_t                                    *loop;
  std::set <index_t, std::less <index_t> >  &boundary_set;
  index_array_t                             *rows_ptr;
  index_t                                   nx;
  index_t                                   ny;
  index_t                                   nz;
  index_t                                   ext_index;
};

template <typename T, typename L, typename BS, typename RP>
build_jacobian_rows_class <T, L>
build_jacobian_rows (mesh_grdecl <T> *mesh, L *l, BS &bs, RP *rp)
{
  return build_jacobian_rows_class <T, L> (mesh, l, bs, rp);
}

template <typename strategy_t, typename loop_t>
struct build_jacobian_cols_class
{
  typedef typename strategy_t::index_t          index_t;
  typedef typename strategy_t::item_t           item_t;
  typedef typename strategy_t::index_array_t    index_array_t;
  typedef typename strategy_t::item_array_t     item_array_t;
  typedef typename strategy_t::rhs_item_array_t rhs_item_array_t;

  typedef mesh_grdecl <strategy_t>                mesh_t;
  typedef typename mesh_t::element_t              element_t;
  typedef typename mesh_t::plane_t                plane_t;
  typedef typename mesh_t::element_zcorn_index_t  element_zcorn_index_t;

  build_jacobian_cols_class (mesh_t *mesh, loop_t *loop, index_array_t *rows_ptr, index_array_t *cols_ind,
    index_array_t &cols_ind_transmis, rhs_item_array_t &values_transmis,
    index_array_t &matrix_block_idx_minus, index_array_t &matrix_block_idx_plus)
  : mesh (mesh)
  , loop (loop)
  , rows_ptr (rows_ptr)
  , cols_ind (cols_ind)
  , cols_ind_transmis (cols_ind_transmis)
  , values_transmis (values_transmis)
  , matrix_block_idx_minus (matrix_block_idx_minus)
  , matrix_block_idx_plus (matrix_block_idx_plus)
  , nx (mesh->nx)
  , ny (mesh->ny)
  , nz (mesh->nz)
  {
    //curIndex.assign (rows_ptr->begin (), rows_ptr->end ());
    index_t i;
    index_t n = mesh->n_active_elements;
    
    
    rows_ptr_tmp.resize (n);
    for (i = 0; i < n; ++i)
      {
        // let it point to first non-diagonal column in each row
        rows_ptr_tmp[i] = (*rows_ptr)[i] + 1;
        (*cols_ind)[(*rows_ptr)[i]]= i;
      }
  }
  
  void
  prepare (index_t i, index_t j, index_t k)
  {
    ext_index  = i + j * nx + k * nx * ny;
    int_index  = mesh->convert_ext_to_int (ext_index);

    mesh->calc_element(i, j, k, element);
    center1  = element.get_center();
  }
  
    
  void
  change_jac_and_flux_conn( const index_t ext_index1, const index_t ext_index2, item_t tran)
  {
    index_t index1 = mesh->convert_ext_to_int (ext_index);
    index_t index2 = mesh->convert_ext_to_int (ext_index2);
    
    (*cols_ind)[rows_ptr_tmp[index1]] = index2;
    (*cols_ind)[rows_ptr_tmp[index2]] = index1;

    matrix_block_idx_minus[loop->con_num] = (*rows_ptr)[index1];
    matrix_block_idx_minus[loop->con_num + 1] = rows_ptr_tmp[index1];
    matrix_block_idx_plus[loop->con_num + 1] = (*rows_ptr)[index2];
    matrix_block_idx_plus[loop->con_num] = rows_ptr_tmp[index2];
    
    ++rows_ptr_tmp[index1];
    ++rows_ptr_tmp[index2];
    
    cols_ind_transmis[loop->con_num++] = index1;
    cols_ind_transmis[loop->con_num++] = index2;
    
    values_transmis.push_back(tran);
    values_transmis.push_back(-tran);
  }

  bool
  check_column_adjacency (index_t i, index_t j)
  {
    return loop->is_column_adjacent[i + j * nx];
  }
  
  void
  add_boundary (index_t)
  {
  }

  void
  change_by_x (index_t i, index_t j, index_t k, index_t ext_index2, bool is_adjacent)
  {
    item_t tran;
    plane_t plane1;
    element_t element2;
    fpoint3d center2  = element2.get_center ();
    
    element.get_plane (x_axis_plus, plane1);
    mesh->calc_element(i, j, k, element2);
     
    if (is_adjacent)
      {
        tran = mesh->calc_tran (ext_index, ext_index2, plane1, center1, center2, along_dim1);
      }
    else
      {
        plane_t plane2;
        element2.get_plane (x_axis_minus, plane2);
        tran = mesh->calc_tran (ext_index, ext_index2, plane1, center1, center2, along_dim1, &plane2);
      }
    
    if (tran != 0)
      {
        change_jac_and_flux_conn (ext_index, ext_index2, tran);
      }  
    /*  
    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index < 1000)
          fprintf (fp, " %d [%d;%d;%d] (%lf)", ext_index2, i, j, k, tran);
    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX 
    */
  }

  void
  change_by_y (index_t i, index_t j, index_t k, index_t ext_index2, bool is_adjacent)
  {
    item_t tran;
    
    plane_t plane1;
    element_t element2;
    fpoint3d center2  = element2.get_center ();
    
    element.get_plane (y_axis_plus, plane1);
    mesh->calc_element(i, j, k, element2);
     
    if (is_adjacent)
      {
        tran = mesh->calc_tran (ext_index, ext_index2, plane1, center1, center2, along_dim2);
      }
    else
      {
        plane_t plane2;
        element2.get_plane (y_axis_minus, plane2);
        tran = mesh->calc_tran (ext_index, ext_index2, plane1, center1, center2, along_dim2, &plane2);
      }
     
    if (tran != 0)
      {
        change_jac_and_flux_conn (ext_index, ext_index2, tran);
      }
     /* 
    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index < 1000)
          fprintf (fp, " %d [%d;%d;%d] (%lf)", ext_index2, i, j, k, tran);
    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX 
    */
  }

  void
  change_by_z (index_t i, index_t j, index_t k, index_t ext_index2, bool is_adjacent)
  {
    item_t tran;
    
    plane_t plane1;
    element_t element2;
    fpoint3d center2  = element2.get_center ();
    
    element.get_plane (z_axis_plus, plane1);
    mesh->calc_element(i, j, k, element2);
    
    tran = mesh->calc_tran (ext_index, ext_index2, plane1, center1, center2, along_dim3);
         
    if (tran != 0)
      {
        change_jac_and_flux_conn (ext_index, ext_index2, tran);
      }
    /*  
    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index < 1000)
          fprintf (fp, " %d [%d;%d;%d] (%lf)", ext_index2, i, j, k, tran);
    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
    */
  }

  mesh_t              *mesh;
  loop_t              *loop;
  index_array_t       *rows_ptr;
  index_array_t       *cols_ind;
  index_array_t       &cols_ind_transmis;
  rhs_item_array_t    &values_transmis;
  index_array_t       &matrix_block_idx_minus;
  index_array_t       &matrix_block_idx_plus;

  index_t             nx;
  index_t             ny;
  index_t             nz;

  // points to current column position for each row
  // while adding new columns
  index_array_t       rows_ptr_tmp; 

  index_t             ext_index; // index1
  index_t             int_index; // index_ijk

  element_t           element;
  fpoint3d            center1;

};

template <typename M, typename L, typename RP, typename CI, typename CT, typename FC>
build_jacobian_cols_class <M, L>
build_jacobian_cols (mesh_grdecl <M> *m, L *l, RP *rp, CI *ci, CT &conn_trans, FC &flux_conn)
{
  return build_jacobian_cols_class <M, L> (m, l, rp, ci,
    conn_trans->get_cols_ind (), conn_trans->get_values (),
    flux_conn->get_matrix_block_idx_minus (), flux_conn->get_matrix_block_idx_plus ());
}

template <typename strategy_t>
struct build_jacobian_and_flux : boost::noncopyable
{
  typedef typename strategy_t::index_t          index_t;
  typedef typename strategy_t::item_t           item_t;
  typedef typename strategy_t::index_array_t    index_array_t;
  typedef typename strategy_t::item_array_t     item_array_t;

  typedef mesh_grdecl <strategy_t>                mesh_t;
  typedef typename mesh_t::plane_t                plane_t;
  typedef typename mesh_t::element_zcorn_index_t  element_zcorn_index_t;

  build_jacobian_and_flux (mesh_grdecl <strategy_t> *mesh)
  : mesh (mesh)
  , nx (mesh->nx)
  , ny (mesh->ny)
  , nz (mesh->nz)
  , con_num (0)
  {
    is_column_adjacent.assign (nx * ny, true);
  }

  template <typename loop_body_t>
  void
  cell_loop (loop_body_t loop_body)
  {
    index_t ext_index1, ext_index2;
    index_t last_k_x, last_k_y, k_x, k_y;
    bool is_adjacent;
    index_t n_adj_elems = 0, n_non_adj_elems = 0;
    
    element_zcorn_index_t zcorn_index1, zcorn_index2; 
    array_float16_t &zcorn = mesh->zcorn_array;
    for (index_t i = 0; i < nx; ++i)
      {
        for (index_t j = 0; j < ny; ++j)
          {
            is_adjacent = loop_body.check_column_adjacency (i, j);
            
            if (is_adjacent)
              {   
                
                // simple loop
                
                for (index_t k = 0; k < nz; ++k)
                  {
                    ext_index1  = i + j * nx + k * nx * ny;
                    
                    //skip non-active cells
                    if (!mesh->sp_actnum[ext_index1])
                      continue;
                      
                    n_adj_elems++;
                    
                    /*
                    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
                    if (ext_index1 < 1000)
                      fprintf (fp, "\n%d [%d;%d;%d]", ext_index1, i, j, k);
                    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
                    */
                      
                    loop_body.prepare (i, j, k);
                    
                    if (i+1 < nx && mesh->sp_actnum[ext_index1 + 1])
                      {
                        loop_body.change_by_x (i + 1, j, k, ext_index1 + 1, true);
                      }
                      
                    if (j+1 < ny && mesh->sp_actnum[ext_index1 + nx])
                      {
                        loop_body.change_by_y (i, j + 1, k, ext_index1 + nx, true);
                      }
                     
                    if (k+1 < nz && mesh->sp_actnum[ext_index1 + nx * ny])
                      {
                        loop_body.change_by_z (i, j, k + 1, ext_index1 + nx * ny, true);
                      }
                  }
              }
            else
              {
                // complicated loop
                
                last_k_x = 0;
                last_k_y = 0;
                
                 /*                             X
                 *                    0+-------+1
                 *                    /|     / |
                 *                  /  |   /   |
                 *               2/-------+3   |
                 *              Y |   4+--|----+5
                 *                |   /Z  |   /
                 *                | /     | /
                 *              6 /-------/7
                 */
                            
                for (index_t k = 0; k < nz; ++k)
                  {
                    ext_index1  = i + j * nx + k * nx * ny;
                    
                    //skip non-active cells
                    if (!mesh->sp_actnum[ext_index1])
                      continue;
                    
                    n_non_adj_elems++;
                    
                    /*
                    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
                    if (ext_index1 < 1000)
                      fprintf (fp, "\n%d [%d;%d;%d]", ext_index1, i, j, k);
                    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
                    */
                    
                    mesh->get_element_zcorn_index (i, j, k, zcorn_index1);
                    
                    loop_body.prepare (i, j, k);
                    
                    if (i + 1 < nx) 
                      {
                        k_x = last_k_x - 1;
                        
                        // search first possible neighbour
                        do
                          {
                            k_x++;
                            mesh->get_element_zcorn_index (i + 1, j, k_x, zcorn_index2);
                          }
                        while ((k_x < nz) && ((zcorn[zcorn_index1[1]] >= zcorn[zcorn_index2[4]]) && (zcorn[zcorn_index1[3]] >= zcorn[zcorn_index2[6]])));
                        
                        // calc all neihbours
                        while ((k_x < nz) && ((zcorn[zcorn_index1[5]] > zcorn[zcorn_index2[0]]) || (zcorn[zcorn_index1[7]] > zcorn[zcorn_index2[2]])))
                          {
                            if ((zcorn[zcorn_index1[5]] >= zcorn[zcorn_index2[4]]) && (zcorn[zcorn_index1[7]] >= zcorn[zcorn_index2[6]]))
                              {
                                // this (i + 1, j, k_x) neigbour won`t touch next (i, j, k + 1) element
                                last_k_x = k_x + 1;
                              }
                              
                            ext_index2 = ext_index1 + 1 + (k_x - k) * nx * ny;
                            
                            if (mesh->sp_actnum[ext_index2])
                              {
                                loop_body.change_by_x (i + 1, j, k_x, ext_index2, false);
                              }
                            k_x++;
                            mesh->get_element_zcorn_index (i + 1, j, k_x, zcorn_index2);
                          }
                      }
                      
                    if (j + 1 < ny) 
                      {
                        k_y = last_k_y - 1;
                        
                        // search first possible neighbour
                        do
                          {
                            k_y++;
                            mesh->get_element_zcorn_index (i, j + 1, k_y, zcorn_index2);
                          }
                        while ((k_y < nz) && ((zcorn[zcorn_index1[2]] >= zcorn[zcorn_index2[4]]) && (zcorn[zcorn_index1[3]] >= zcorn[zcorn_index2[5]])));
                        
                        // calc all neighbours
                        while ((k_y < nz) && ((zcorn[zcorn_index1[6]] > zcorn[zcorn_index2[0]]) || (zcorn[zcorn_index1[7]] > zcorn[zcorn_index2[1]])))
                          {
                            if ((zcorn[zcorn_index1[6]] >= zcorn[zcorn_index2[4]]) && (zcorn[zcorn_index1[7]] >= zcorn[zcorn_index2[5]]))
                              {
                                // this (i, j + 1, k_y) neigbour won`t touch next (i, j, k + 1) element
                                last_k_y = k_y + 1;
                              }
                              
                            ext_index2 = ext_index1 + ny + (k_y - k) * nx * ny;
                            
                            if (mesh->sp_actnum[ext_index2])
                              {
                                loop_body.change_by_y (i, j + 1, k_y, ext_index2, false);
                              }
                            k_y++;
                            mesh->get_element_zcorn_index (i, j + 1, k_y, zcorn_index2);
                          }
                      }
                      
                    if (k + 1 < nz && mesh->sp_actnum[ext_index1 + nx * ny])
                      {
                        loop_body.change_by_z (i, j, k + 1, ext_index1 + nx * ny, true);
                      }    
                  }
              }
          }
      }
    BOSWARN (section::mesh, level::warning)<< boost::format ("MESH_GRDECL: elements adj %d, non-adj %d, total %d") \
            % n_adj_elems % n_non_adj_elems % (n_adj_elems + n_non_adj_elems) << bs_end;  
    BOSWARN (section::mesh, level::warning)<< boost::format ("MESH_GRDECL: number of tran calcs is %d") % n_tran_calc << bs_end;  
  }

  mesh_grdecl <strategy_t>    *mesh;
  index_t                     nx;
  index_t                     ny;
  index_t                     nz;
  shared_vector <bool>        is_column_adjacent;
  index_t                     con_num;
};


template<class strategy_t>
int mesh_grdecl<strategy_t>::build_jacobian_and_flux_connections_add_boundary (const sp_bcsr_t &jacobian,
                                                                               const sp_flux_conn_iface_t &flux_conn,
                                                                               index_array_t &boundary_array)
{
  write_time_to_log init_time ("Mesh transmissibility calculation", ""); 
  
  n_connections = 0;
  
  index_t max_size = n_active_elements;
  jacobian->get_cols_ind().clear();
  jacobian->get_rows_ptr().clear();
  jacobian->init_struct(max_size,max_size, max_size);
  index_array_t* rows_ptr = &jacobian->get_rows_ptr();

  sp_bcsr_t conn_trans;

  (*rows_ptr)[0] = 0;

  std::vector<bool> is_butting(nx*ny,false);

  std::set<index_t, std::less<index_t> > boundary_set;

  build_jacobian_and_flux <strategy_t> build_jacobian (this);
  
 #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX     
  fp = fopen ("transmiss.out", "w");
 #endif BS_MESH_WRITE_TRANSMISS_MATRIX   

  //first step - define and fill - rows_ptr (jacobian)
  build_jacobian.cell_loop (build_jacobian_rows (this, &build_jacobian, boundary_set, rows_ptr));

  //////jacobian//////////////////////
  //sum rows_ptr
  for (size_t i = 1; i < rows_ptr->size(); ++i)
    {
      (*rows_ptr)[i]++;
      (*rows_ptr)[i] += (*rows_ptr)[i-1];
    }
  //create cols_ind
  index_array_t* cols_ind = &jacobian->get_cols_ind();

  cols_ind->resize((*rows_ptr)[rows_ptr->size()-1],-1);


  ////////transmis/////////////////////////
  index_t cols_ind_n = (index_t)cols_ind->size();


  index_t con_num = (cols_ind_n-max_size)/2;//connection number
  n_connections = con_num;
  conn_trans = flux_conn->get_conn_trans();
  conn_trans->init_struct(con_num,2*con_num,2*con_num);

  index_array_t *rows_ptr_transmis = &conn_trans->get_rows_ptr();
  index_array_t &matrix_block_idx_minus = flux_conn->get_matrix_block_idx_minus ();
  index_array_t &matrix_block_idx_plus = flux_conn->get_matrix_block_idx_plus ();

  matrix_block_idx_minus.resize(con_num*2,-1);
  matrix_block_idx_plus.resize(con_num*2,-1);

  if (!con_num)
    {
      for (int i = 0; i < cols_ind_n; ++i)
        (*cols_ind)[i] = i;
      return 0;
    }

  for (index_t i = 0; i < con_num+1; ++i)
    (*rows_ptr_transmis)[i] = i*2;

  //second step - fill and define cols_ind
  build_jacobian.con_num = 0;
  build_jacobian.cell_loop (build_jacobian_cols (this, &build_jacobian, rows_ptr, cols_ind, conn_trans, flux_conn));


  boundary_array.assign(boundary_set.begin(), boundary_set.end());
  
  #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX  
    fflush (fp);
    fclose (fp);
  #endif BS_MESH_WRITE_TRANSMISS_MATRIX 

  return (int) boundary_array.size();
}

#ifdef _HDF5_MY
template<class strategy_t>
int mesh_grdecl<strategy_t>::create_array_hdf5(const char *dataset_name, H5::H5File &file_hdf5, H5::DataSet **dataset)
{
  // creating dataset_coords
  hsize_t dims[] = {0};
  hsize_t dims_max[] = {H5S_UNLIMITED};
  H5::DataSpace dataspace_coords(1, dims, dims_max);
  // set the dataset to be chunked
  H5::DSetCreatPropList cparms;
  hsize_t chunk_dims[] ={1};
  cparms.setChunk(1, chunk_dims);
  *dataset = new H5::DataSet(file_hdf5.createDataSet(dataset_name, H5::PredType::NATIVE_FLOAT, dataspace_coords, cparms));
  return 0;
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_activs_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt)
{
  return true;
  /*
  // turn off error printing
  H5E_auto2_t old_func;
  void *old_client_data;
  H5::Exception::getAutoPrint(old_func, &old_client_data);
  H5::Exception::dontPrint(); // turn off error printing

  // check if file exists
  hid_t file_hid;
  file_hid = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_hid < 0)
  return false;
  H5Fclose(file_hid);

  H5::H5File file(file_name, H5F_ACC_RDWR);

  try
  {
  actnum_array.resize((it - is + 1) * (jt - js  + 1) * (kt - ks + 1));
  H5::DataSet dataset = file.openDataSet("activ");
  hsize_t dims_memory[] = {((it-is+1))*(jt - js + 1)};
  H5::DataSpace dataspace_memory(1, dims_memory);
  H5::DataSpace dataspace_file = dataset.getSpace();
  hsize_t count[] = {jt - js + 1};
  hsize_t stride[] = {nx};
  hsize_t block[] = {it - is + 1};
  for (int k = 0; k < kt - ks + 1; k++)
  {
  hsize_t start[] = {is + js * nx + (k + ks) * nx * ny};
  dataspace_file.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
  dataset.read(&actnum_array[k * count[0] * block[0]], H5::PredType::NATIVE_INT, dataspace_memory, dataspace_file);
  }
  }
  catch (H5::FileIException exception)
  {
  exception.printError();
  return false;
  }

  H5::Exception::setAutoPrint(old_func, old_client_data);
  n_active_elements = accumulate(actnum_array.begin(),actnum_array.end(),0);
  return true;
  */
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_cube_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt)
{
  // turn off error printing
  H5E_auto2_t old_func;
  void *old_client_data;
  H5::Exception::getAutoPrint(old_func, &old_client_data);
  H5::Exception::dontPrint();

  // check if file exists
  hid_t file_hid;
  file_hid = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_hid < 0)
    return false;
  H5Fclose(file_hid);

  H5::H5File file(file_name, H5F_ACC_RDONLY);

  try
    {
      // read "coord" array
      H5::DataSet dataset_coord = file.openDataSet("coord");
      hsize_t dims_memory_coord[] = {6 * (it - is + 2) * (jt - js  + 2)};
      H5::DataSpace dataspace_memory_coord(1, dims_memory_coord);
      H5::DataSpace dataspace_file_coord = dataset_coord.getSpace();
      coord_array.clear();
      vector<item_t> buf_array(6 * (it - is + 2) * (jt - js  + 2));
      // hyperslab settings
      hsize_t count_coord[] = {jt - js + 2};
      hsize_t start_coord[] = {(is + js * (nx + 1)) * 6};
      hsize_t stride_coord[] = {6 * (nx + 1)};
      hsize_t block_coord[] = {(it - is + 2) * 6};
      // we read "coord" array by 1 read operation
      dataspace_file_coord.selectHyperslab(H5S_SELECT_SET, count_coord, start_coord, stride_coord, block_coord);
      dataset_coord.read(&buf_array[0], H5::PredType::NATIVE_FLOAT, dataspace_memory_coord, dataspace_file_coord);
      for (int i = 0; i < (it - is + 2) * (jt - js  + 2); i++)
        coord_array.push_back(coordElem (fpoint3d(buf_array[i*6],buf_array[i*6+1],buf_array[i*6+2]),
                                         fpoint3d(buf_array[i*6+3],buf_array[i*6+4],buf_array[i*6+5])));

      // read "zcorn" array
      H5::DataSet dataset_zcorn = file.openDataSet("zcorn");
      hsize_t dims_memory_zcorn[] = {8 * (it - is + 1) * (jt - js + 1)};
      H5::DataSpace dataspace_memory_zcorn(1, dims_memory_zcorn);
      H5::DataSpace dataspace_file_zcorn = dataset_zcorn.getSpace();
      zcorn_array.resize(8 * (it - is + 1) * (jt - js  + 1) * (kt - ks + 1));
      // hyperslab settings (common for each plane)
      hsize_t count_zcorn[] = {4 * (jt - js + 1)};
      hsize_t stride_zcorn[] = {2 * nx};
      hsize_t block_zcorn[] = {2 * (it - is + 1)};
      // we read "zcorn" array by planes - kt-ks+1 read operations
      for (int k = 0; k < kt - ks + 1; k++)
        {
          // determine start array for each plane individually
          hsize_t start_zcorn[] = {2 * is + 4 * js * nx + 8 * k * nx * ny};
          dataspace_file_zcorn.selectHyperslab(H5S_SELECT_SET, count_zcorn, start_zcorn, stride_zcorn, block_zcorn);
          dataset_zcorn.read(&zcorn_array[k * count_zcorn[0] * block_zcorn[0]], H5::PredType::NATIVE_FLOAT, dataspace_memory_zcorn, dataspace_file_zcorn);
        }
    }
  catch (H5::Exception exception)
    {
      exception.printError();
      coord_array.clear();
      zcorn_array.clear();
      return false;
    }

  // enable error printing
  H5::Exception::setAutoPrint(old_func, old_client_data);

  return true;
}

template<class strategy_t>
int mesh_grdecl<strategy_t>::append_array_hdf5(const item_t *arr, size_t arr_length, H5::DataSet *dataset)
{
  // determine new dimensions of dataset
  hsize_t dims_old[1];
  hsize_t dims_memory[1] = {arr_length};
  H5::DataSpace dataspace_file = dataset->getSpace();
  dataspace_file.getSimpleExtentDims(dims_old);
  hsize_t dims_new[1] = {dims_old[0] + dims_memory[0]};
  // extend dataset
  dataset->extend(dims_new);
  dataspace_file = dataset->getSpace();
  H5::DataSpace dataspace_memory(1, dims_memory);
  // select hyperslab to write
  hsize_t count[1] = {dims_memory[0]};
  hsize_t start[1] = {dims_old[0]};
  dataspace_file.selectHyperslab(H5S_SELECT_SET, count, start);
  // write array
  dataset->write(arr, H5::PredType::NATIVE_FLOAT, dataspace_memory, dataspace_file);

  return 0;
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_cube_with_hdf5_swap(const char* file_name)
{
  try
    {
      fstream file(file_name, ios::in);
      if (!file.is_open())
        return false;

      H5::H5File file_hdf5(filename_hdf5, H5F_ACC_TRUNC);
      H5::DataSet *dataset_coords = 0;
      H5::DataSet *dataset_zcorn = 0;

      char buf[100];
      item_array_t buf_array;
      bool array_end;

      /////////////////////////////////
      // read "coord"
      /////////////////////////////////
      while (!file.eof())
        {
          file >> buf;
          if (strcmp(buf, "COORD") == 0)
            break;
        }
      create_array_hdf5("coord", file_hdf5, &dataset_coords);
      array_end = false;
      while (!file.eof())
        {
          // check if it is array end
          file >> buf;
          if (buf[0] == '/')
            array_end = true;

          if (array_end != true)
            {
              buf_array.push_back((float) atof(buf));
            }

          // if maximum space is used
          // or if it is array end
          // write "coords" to the swap file
          if ((buf_array.size() >= HDF5_MAX_SPACE)  ||
              ((array_end == true) && (buf_array.size() > 0))
             )
            {
              // save coord_array to hdf5 swap file
              append_array_hdf5(&buf_array[0], buf_array.size(), dataset_coords);
              buf_array.clear();
            }

          if (array_end == true)
            break;
        }
      delete dataset_coords;

      /////////////////////////////////
      // read "zcorn"
      /////////////////////////////////
      while (!file.eof())
        {
          file >> buf;

          if (strcmp(buf, "ZCORN") == 0)
            break;
        }
      // creating dataset_zcorn
      create_array_hdf5("zcorn", file_hdf5, &dataset_zcorn);
      array_end = false;
      while (!file.eof())
        {
          // check if it is array end
          file >> buf;
          if (buf[0] == '/')
            array_end = true;

          if (array_end != true)
            {
              buf_array.push_back((float) atof(buf));
            }

          // if maximum space is used
          // or if it is array end
          // write "coords" to the swap file
          if ((buf_array.size() >= HDF5_MAX_SPACE)  ||
              ((array_end == true) && (buf_array.size() > 0))
             )
            {
              // save coord_array to hdf5 swap file
              append_array_hdf5(&buf_array[0], buf_array.size(), dataset_zcorn);
              buf_array.clear();
            }

          if (array_end == true)
            break;
        }
      delete dataset_zcorn;
      return true;
    }
  catch (H5::Exception exception)
    {
      exception.printError();
      return false;
    }
}
#endif

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_actnum(const char* file_name)
{
#if 0
  fstream file(file_name,  ios::in);
  char buf[255];
  char *start_ptr,*end_ptr;
  if (!file.is_open())
    {
      n_active_elements =  accumulate(sp_actnum.begin(),sp_actnum.end(),0);
      return false;
    }
  sp_actnum.clear();

  while (!file.eof())
    {
      file >> buf;
      start_ptr = buf;
      if (strcmp (buf, "/") != 0)
        sp_actnum.push_back((int)strtod (start_ptr, &end_ptr));
      else
        break;
    }
  n_active_elements =  accumulate(sp_actnum.begin(),sp_actnum.end(),0);
#endif
  return true;
}


template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_cube(const char* file_name)
{
  /*
  using namespace std;
  fstream file(file_name,  ios::in);
  char buf[100];
  char *start_ptr, *end_ptr;
  if (!file.is_open())
    return false;

  while (!file.eof())
    {
      file >> buf;
      if (  strcmp (buf, "COORD") == 0)
        break;
    }

  max_x = max_y = -100000;
  min_x = min_y = 1000000;

  index_t i_count = 0;
  const int max_count = 6; //buffer_size
  vector<float> buf_vec(max_count);
  while (!file.eof())
    {
      //coordElem curElem;
      file >> buf;

      if (buf[0] != '/')
        {
          buf_vec[i_count++] = (float)atof(buf);
          if (i_count == max_count) //swap_buffer
            {
              i_count = 0;
              coordElem curElem = coordElem(fpoint3d(buf_vec[0],buf_vec[1],buf_vec[2]),
                                            fpoint3d(buf_vec[3],buf_vec[4],buf_vec[5]));
              coord_array.push_back(curElem);

              //looking for max&min coordinates
              if (curElem.pStart.x > max_x)   max_x = curElem.pStart.x;
              if (curElem.pEnd.x > max_x)     max_x = curElem.pEnd.x;

              if (curElem.pStart.y > max_y)   max_y = curElem.pStart.y;
              if (curElem.pEnd.y > max_y)     max_y = curElem.pEnd.y;

              if (curElem.pStart.x < min_x)   min_x = curElem.pStart.x;
              if (curElem.pEnd.x < min_x)     min_x = curElem.pEnd.x;

              if (curElem.pStart.y < min_y)   min_y = curElem.pStart.y;
              if (curElem.pEnd.y < min_y)     min_y = curElem.pEnd.y;
            }
        }
      else
        break;
    }

  while (!file.eof())
    {
      file >> buf;
      if (  strcmp (buf, "ZCORN") == 0)
        break;
    }

  while (!file.eof())
    {
      file >> buf;
      start_ptr = buf;
      if (strcmp (buf, "/") != 0)
        zcorn_array.push_back((float)strtod (start_ptr, &end_ptr));
      else
        break;
    }
  min_z = *(std::min_element(zcorn_array.begin(),zcorn_array.end()));
  max_z = *(std::max_element(zcorn_array.begin(),zcorn_array.end()));
  return true;
  */
  return false;
}


BS_INST_STRAT(mesh_grdecl);