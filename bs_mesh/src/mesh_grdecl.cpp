/*! \file mesh_grdecl.cpp
\brief This file implement class for working with grid_ecllipse
\author Iskhakov Ruslan
\date 2008-05-20 */
#include "bs_mesh_stdafx.h"
#include "mesh_grdecl.h"
#include <iterator>

using namespace grd_ecl;
using namespace blue_sky;

const char filename_hdf5[] = "grid_swap.h5";

#ifdef BS_MESH_WRITE_TRANSMISS_MATRIX
  FILE*  fp;
#endif //BS_MESH_WRITE_TRANSMISS_MATRIX 

struct mesh_grdecl::inner {
    // shorter aliases
    typedef t_long int_t;
    typedef t_double fp_t;
    typedef t_float fp_stor_t;

    void init_minmax(mesh_grdecl& m) const {
        int_t i, n;
        fp_stor_t *it;

        // init ZCORN
        m.min_z = *(std::min_element(zcorn_->begin(), zcorn_->end()));
        m.max_z = *(std::max_element(zcorn_->begin(), zcorn_->end()));

        n = (int_t) coord_->size();

        m.max_x = m.min_x = m.coord_array[0];
        m.max_y = m.min_y = m.coord_array[1];

        for (i = 0; i < n; i += 6) {
            it = m.coord_array + i;
            // move matching points apart
            if (it[2] == it[5])
                it[5] += 1.0f;

            //looking for max&min coordinates
            if (m.min_x > it[0]) m.min_x = it[0];
            if (m.min_x > it[3]) m.min_x = it[3];

            if (m.min_y > it[1]) m.min_y = it[1];
            if (m.min_y > it[4]) m.min_y = it[4];

            if (m.max_x < it[0]) m.max_x = it[0];
            if (m.max_x < it[3]) m.max_x = it[3];

            if (m.max_y < it[1]) m.max_y = it[1];
            if (m.max_y < it[4]) m.max_y = it[4];
        }
    }

    // hold reference to coord and czron arrays if generated internally
    spv_float coord_;
    spv_float zcorn_;
};

mesh_grdecl::mesh_grdecl ()
    : pinner_(new inner), coord_array(0), zcorn_array(0)
{}

void mesh_grdecl::init_props(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
    if(coord && coord->size()) {
        pinner_->coord_ = coord;
        coord_array = &pinner_->coord_->ss(0);
    }
    if(zcorn && zcorn->size()) {
        pinner_->zcorn_ = zcorn;
        zcorn_array = &pinner_->zcorn_->ss(0);
    }
    this->nx = nx;
    this->ny = ny;
    this->nz = (zcorn->size() / nx / ny) >> 3;
    this->n_elements = nx * ny * nz;

    // postinit
    pinner_->init_minmax(*this);
}

void mesh_grdecl::clear() {
    // free coord
    if(pinner_->coord_ && pinner_->coord_->size()) {
        pinner_->coord_->clear();
        pinner_->coord_.release();
        coord_array = NULL;
    }

    // free zcorn
    if(pinner_->zcorn_ && pinner_->zcorn_->size()) {
        pinner_->zcorn_->clear();
        pinner_->zcorn_.release();
        zcorn_array = NULL;
    }

    this->nx = 0;
    this->ny = 0;
    this->nz = 0;
    this->n_elements = 0;

    // zero bounds
    this->max_x = this->max_y = this->max_z = 0;
    this->min_x = this->min_y = this->min_z = 0;
}

void mesh_grdecl::init_props(t_long nx, t_long ny, t_long nz, spv_float dx, spv_float dy, spv_float dz) {
    // generate COORD & ZCORN
    std::pair< spv_float, spv_float > cz = gen_coord_zcorn(nx, ny, nz, dx, dy, dz);
    init_props(nx, ny, cz.first, cz.second);
}

void mesh_grdecl::init_props(const sp_hdm_t hdm)
{
  base_t::init_props (hdm);
  spv_float data_array;
  //t_long i, n;
  //t_float *it;
  
  // init ZCORN
  data_array = hdm->get_pool ()->get_fp_data("ZCORN");
  if (data_array->size()) {
      pinner_->zcorn_ = data_array;
      zcorn_array = &(*data_array)[0];
  }
  // init COORD
  data_array = hdm->get_pool ()->get_fp_data("COORD");
  if (data_array->size()) {
      pinner_->coord_ = data_array;
      coord_array = &(*data_array)[0];
  }

  // postinit
  pinner_->init_minmax(*this);
  hdm->get_prop()->add_property_f(0, "min_x", std::string ("Minimum X coordinate for current mesh"));
  hdm->get_prop()->add_property_f(0, "min_y", std::string ("Minimum Y coordinate for current mesh"));
  hdm->get_prop()->add_property_f(0, "min_z", std::string ("Minimum Z coordinate for current mesh"));

  hdm->get_prop()->add_property_f(0, "max_x", std::string ("Maximum X coordinate for current mesh"));
  hdm->get_prop()->add_property_f(0, "max_y", std::string ("Maximum Y coordinate for current mesh"));
  hdm->get_prop()->add_property_f(0, "max_z", std::string ("Maximum Z coordinate for current mesh"));

  hdm->get_prop()->set_f("min_x", min_x);
  hdm->get_prop()->set_f("min_y", min_y);
  hdm->get_prop()->set_f("min_z", min_z);

  hdm->get_prop()->set_f("max_x", max_x);
  hdm->get_prop()->set_f("max_y", max_y);
  hdm->get_prop()->set_f("max_z", max_z);
}


void mesh_grdecl::check_data() const
{
  write_time_to_log init_time ("Mesh check", ""); 
  base_t::check_data ();

  if (min_x < 0)
    bs_throw_exception (boost::format ("min_x = %d is out of range")% min_x);
  if (min_y < 0)
    bs_throw_exception (boost::format ("min_y = %d is out of range")% min_y);
  if (min_z < 0)
    bs_throw_exception (boost::format ("min_z = %d is out of range")% min_z);

  if (!coord_array)
    bs_throw_exception ("COORD array is not initialized");
  if (!zcorn_array)
    bs_throw_exception ("ZCORN array is not initialized");

  // 1. check if all cells are convex

  element_t element;
  t_long wrong_cells = 0;
  t_int *actnum = actnum_array->data ();

  for (t_long i = 0; i < nx; ++i)
    for (t_long j = 0; j < ny; ++j)
      {
        for (t_long k = 0; k < nz; ++k)
          {
            t_long index = k + j * nz + i * ny * nz;

            // miss inactive blocks
            if (actnum[index] || k == 0 || k == nz - 1)
              {
                calc_element (i, j, k, element);
                mesh_element3d::corners_t corns = element.get_corners();
                /*
                if (k == 0)
                  {
                    t_long cindex = 6 * (i + j * (nx + 1));

                    coord_array[cindex] = corns[0].x;
                    coord_array[cindex + 1] = corns[0].y;
                    coord_array[cindex + 2] = corns[0].z;

                    cindex += 6;


                    calc_corner_point (zcorn_array[index1 + 1], &coord_array[(iCOORD + 1) * 6], corners[1]);
                    calc_corner_point (zcorn_array[index1 + 2 * nx], &coord_array[(iCOORD + (nx + 1)) * 6], corners[2]);
                    calc_corner_point (zcorn_array[index1 + 2 * nx + 1], &coord_array[(iCOORD + (nx + 1) + 1) * 6], corners[3]);

                */
                if (actnum[index])
                  {
                    // check X
                    if (((corns[1].x - corns[0].x) * (corns[3].x - corns[2].x) < 0) ||
                        ((corns[1].x - corns[0].x) * (corns[5].x - corns[4].x) < 0) ||
                        ((corns[1].x - corns[0].x) * (corns[7].x - corns[6].x) < 0))
                      {
                        actnum[index] = 0;
                        wrong_cells ++;
                        continue;
                      }

                    // check Y
                    if (((corns[2].y - corns[0].y) * (corns[3].y - corns[1].y) < 0) ||
                        ((corns[2].y - corns[0].y) * (corns[7].y - corns[5].y) < 0) ||
                        ((corns[2].y - corns[0].y) * (corns[6].y - corns[4].y) < 0))
                      {
                        actnum[index] = 0;
                        wrong_cells ++;
                        continue;
                      }
                
                
                    // check Z
                    if (((corns[4].z - corns[0].z) * (corns[5].z - corns[1].z) < 0) ||
                        ((corns[4].z - corns[0].z) * (corns[6].z - corns[2].z) < 0) ||
                        ((corns[4].z - corns[0].z) * (corns[7].z - corns[3].z) < 0))
                      {
                        actnum[index] = 0;
                        wrong_cells ++;
                        continue;
                      }
                  }
              }
        }
    }

  if (wrong_cells)
    BOSOUT (section::mesh, level::medium) << "% wrong (nonconvex) cells found! Marked inactive." << wrong_cells << bs_end;

  // 2. check for ZCORN intersections (and correct them)
/*
  for (t_long i = 0; i < 2 * nx; ++i)
    for (t_long j = 0; j < 2 * ny; ++j)
      for (t_long k = 1; k < nz; ++k)
        {
          t_long z_index = i + 2 * nx * j + 4 * nx * ny * (2 * k - 1);
          if (zcorn_array[z_index] > zcorn_array[z_index + 4 * nx * ny])
            {
               t_long index1 = k + std::ceil(j / 2) * nz + std::ceil(i / 2) * ny * nz - 1;
               t_long index2 = index1 + 1;

               if (actnum[index1] == 0)
                 zcorn_array[z_index] = zcorn_array[z_index + 4 * nx * ny];
               else if (actnum[index2] == 0)
                 zcorn_array[z_index + 4 * nx * ny] = zcorn_array[z_index]
               else
                 {
                   t_float middle = (zcorn_array[z_index + 4 * nx * ny] + zcorn_array[z_index]) / 2;
                   zcorn_array[z_index + 4 * nx * ny] = zcorn_array[z_index] = middle;
                 }
            }
        };
  */

}

inline void
mesh_grdecl::calc_corner_point(const t_float z, const t_float *coord, fpoint3d_t &p)const
  {
    p.z = z;
    /*
   t_double temp = (p.z-m_Coord.pStart.z)/(m_Coord.pEnd.z-m_Coord.pStart.z);
    p.x = temp *(m_Coord.pEnd.x-m_Coord.pStart.x)+m_Coord.pStart.x;
    p.y = temp *(m_Coord.pEnd.y-m_Coord.pStart.y)+m_Coord.pStart.y;
    */
   t_double temp = (p.z - coord[2]) / (coord[5] - coord[2]);
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


inline void
mesh_grdecl::get_element_zcorn_index (const t_long i, const t_long j, const t_long k, element_zcorn_t_long& element)  const
{
  // mesh_grdecl::element_zcorn_t_long element;
  
  t_long index1 = i * 2 + j * 4 * nx + k * 8 * nx * ny;
  t_long index2 = index1 + 4 * nx * ny;

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

mesh_grdecl::plane_zcorn_t_long
mesh_grdecl::get_plane_zcorn_index (element_zcorn_t_long &element, 
                                                         element_plane_orientation_t orientation)  const
{
   mesh_grdecl::plane_zcorn_t_long plane;
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
  return plane;  
}


 mesh_grdecl::element_t
mesh_grdecl::calc_element (const t_long index) const
  {
    t_long i, j, k;
    element_t element;
    
    inside_to_XYZ (index, i, j, k);
    calc_element (i, j, k, element);
    return element;
  }


void
mesh_grdecl::calc_element (const t_long i, const t_long j, const t_long k, element_t &element) const
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

    t_long index1 = i * 2 + j * 4 * nx + k * 8 * nx * ny;//upper side
    t_long index2 = index1 + 4 * nx * ny;//lower side
    t_long iCOORD = i + j * (nx + 1);

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
  

bool mesh_grdecl::is_small(const t_long i, const t_long j, const t_long k, t_double eps)  const
  {
    if (k >= nz)
      return false;

    t_double dz1, dz2, dz3, dz4; //height for each coord
    //define index
    t_long index1 = i*2+j*4*nx+k*8*nx*ny;   //lower side
    t_long index2 = index1+4*nx*ny;         //upper side
    dz1 = zcorn_array[index2]-zcorn_array[index1];
    dz2 = zcorn_array[index2+1]-zcorn_array[index1+1];
    dz3 = zcorn_array[index2+2*nx]-zcorn_array[index1+2*nx];
    dz4 = zcorn_array[index2+2*nx+1]-zcorn_array[index1+2*nx+1];

    if (dz1 <= eps && dz2 <= eps && dz3 <= eps && dz4 <= eps)
      return true;
    return false;
  }


void
mesh_grdecl::get_plane_crossing_projection_on_all_axis(const plane_t &plane1, const plane_t &plane2, fpoint3d_t &A)const
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


t_double 
mesh_grdecl::get_center_zcorn(const element_zcorn_t_long &element)const
  {
    t_double res = 0.0;
    size_t i, n = element.size();
    
    for (i = 0; i < n; i++)
      res += zcorn_array[element[i]];
      
    return res/n;
  }




int mesh_grdecl::init_ext_to_int()
{
  write_time_to_log init_time ("Mesh initialization", ""); 
  stdv_float volumes_temp(n_elements);

  //tools::save_seq_vector ("actnum.bs.txt").save actnum_array;
  
  // FIXME: was commented
  int splicing_num = splicing(volumes_temp);
  splicing_num;
  
  //check_adjacency (1);
  //tools::save_seq_vector ("active_blocks.bs.txt").save actnum_array;
  
  check_adjacency ();
  
  //make proxy array
  ext_to_int->resize(n_elements);
  ext_to_int->assign(0);
  t_long *ext_to_int_data = ext_to_int->data ();
  size_t n_count = 0;

  t_long nn_active = 0; //number of non-active previous cells
  t_int const *actnum = actnum_array->data ();
  for (t_long i = 0; i < nx; ++i)
    {
      for (t_long j = 0; j < ny; ++j)
        for (t_long k = 0; k < nz; ++k, ++n_count)
          {
            t_long i_index = k + (nz * j) + (i * ny * nz);
            
            if (!actnum[i_index])
              {
                nn_active++;
                ext_to_int_data[n_count] = -1;
              }
            else
              ext_to_int_data[n_count] = i_index-nn_active;
          }
    }

  //tools::save_seq_vector ("ext_to_int.bs.txt").save (ext_to_int);

  init_int_to_ext();
  t_long *int_to_ext_data = int_to_ext->data ();
  
  //fill volume array (except non-active block and using proxy array)
  volumes->resize(n_active_elements);
  t_float *volumes_data = volumes->data ();
  
  // FIXME: init volumes_temp
  for (int i = 0; i < n_active_elements; ++i)
    volumes_data[i] = volumes_temp[int_to_ext_data[i]];
    
  calc_depths();
  

  //bs_throw_exception ("NOT IMPL YET");
  return n_active_elements;
}


void mesh_grdecl::splice_two_blocks (const t_long i, const t_long j, const t_long k, const t_long k1)
{
  t_long index, index1, index2;

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


bool mesh_grdecl::are_two_blocks_close (const t_long i, const t_long j, const t_long k, const t_long k1)
{
  t_long index, index1;
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


bool mesh_grdecl::check_adjacency(int shift_zcorn)
{
  t_long n_adjacent = 0;
   
  t_int const *actnum = actnum_array->data ();
  t_float *zcorn = pinner_->zcorn_->data ();
  for (t_long i = 0; i < nx; ++i)
    for (t_long j = 0; j < ny; ++j)
      for (t_long k = 0; k < nz; ++k)
        {
          
          t_long index = k + j * nz + i * ny * nz;
          
          
          // miss inactive blocks
          if (!actnum [index])
            {
              continue;
            }
                       
          // check blocks adjacency
          t_long zindex = i * 2 + j * 4 * nx + k * 8 *nx * ny;
          t_long zindex1 = zindex + 2; // check next by x
         
          if (i + 1 == nx ||
              (zcorn[zindex + 1] == zcorn[zindex1] &&
              zcorn[zindex +  2 * nx + 1] == zcorn[zindex1 +  2 * nx] &&
              zcorn[zindex + 4 * nx * ny + 1] == zcorn[zindex1 + 4 * nx * ny ] &&
              zcorn[zindex + 4 * nx * ny + 2 * nx + 1] == zcorn[zindex1 + 4 * nx * ny + 2 * nx]))
            {
              
              zindex1 = zindex + 4 * nx; // check next by y
              if (j + 1 == ny ||
                  (zcorn[zindex + 2 * nx] == zcorn[zindex1] &&
                  zcorn[zindex + 2 * nx + 1] == zcorn[zindex1 +  1] &&
                  zcorn[zindex + 4 * nx * ny + 2 * nx] == zcorn[zindex1 + 4 * nx * ny] &&
                  zcorn[zindex + 4 * nx * ny + 2 * nx + 1] == zcorn[zindex1+ 4 * nx * ny + 1]))
                  {
                    n_adjacent++;
                  }
            }
            if (shift_zcorn && (i + j) % 2 == 1)
              {
                // FIXME: precalculate
                zcorn[zindex] += 2;
                zcorn[zindex + 1] += 2;
                zcorn[zindex +  2 * nx] += 2;
                zcorn[zindex +  2 * nx + 1] += 2;
                zcorn[zindex + 4 * nx * ny] += 2;
                zcorn[zindex + 4 * nx * ny + 1] += 2;
                zcorn[zindex + 4 * nx * ny +  2 * nx] += 2;
                zcorn[zindex + 4 * nx * ny +  2 * nx + 1] += 2;
              }
       }
  
  BOSOUT (section::mesh, level::medium) << "  adjacent cells:"<< n_adjacent <<" ("<< n_adjacent * 100 / (n_active_elements)<<"% active)" << bs_end;
  return (n_adjacent == n_active_elements);
}

spv_float mesh_grdecl::get_cell_volumes (const t_long Nx, const t_long Ny, const t_long Nz) const
{
    t_long i, j, k, ind;
    element_t element;

    spv_float volumes = BS_KERNEL.create_object(v_float::bs_type());
    volumes->resize(Nx*Ny*Nz);

    // important: XYZ order

    ind = 0;
    for (k = 0; k < Nz; ++k)
        for (j = 0; j < Ny; ++j)
            for (i = 0; i < Nx; ++i)
              {
                calc_element (i, j, k, element);
                (*volumes)[(ind++)] = element.calc_volume ();
              }
    return volumes;
}

int mesh_grdecl::splicing(stdv_float& volumes_temp)
{
#ifdef _DEBUG
    BS_ASSERT (zcorn_array != 0 && coord_array != 0);
#endif


  t_long nCount = 0;
  t_long n_inactive_orig = 0;
  t_long n_inactive_vol = 0;
  t_long n_incative_splice = 0;

  // FIXME: how to check ranges?
  t_int *actnum = actnum_array->data ();
  t_float const *poro = poro_array->data ();
  for (t_long i = 0; i < nx; ++i)
    for (t_long j = 0; j < ny; ++j)
      {
        t_long small_block_top = -1;
        t_long big_block_top = -1;
        t_double vol_sum = 0.0;
        for (t_long k = 0; k < nz; ++k)
          {
            t_long index = k + j * nz + i * ny * nz;
            
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

            element_t element;
            calc_element (i, j, k, element);
            t_double vol_block = element.calc_volume ();
            t_double vol_block_poro = vol_block * poro[index];

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
                    volumes_temp[big_block_top + j * nz + i * ny * nz] += vol_block;

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
                    for (t_long k1 = k - 1; k1 >= small_block_top; --k1)
                      {
                        splice_two_blocks (i, j, k, k1);
                        n_incative_splice++;
                        // make small block inactive
                        actnum[k1 + j * nz + i * ny * nz] = 0;
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
  
  
  t_long n_total = n_active_elements + n_inactive_orig;
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



int mesh_grdecl::calc_depths ()
{
  depths->resize (n_active_elements);

  t_float *depths_data = depths->data ();
  t_long const *ext_to_int_data = ext_to_int->data ();
  t_int const *actnum = actnum_array->data ();

  t_long index = 0;
  for (t_long i = 0; i < nx; ++i)
    for (t_long j = 0; j < ny; ++j)
      for (t_long k = 0; k < nz; ++k, ++index)
        {
          if (actnum[index])
            {
              // FIXME: check index and array length
              element_zcorn_t_long element;
              get_element_zcorn_index(i, j, k, element);
              depths_data[ext_to_int_data[index]] = get_center_zcorn(element);
            }
        }
  return 0;
}


static int n_tran_calc = 0;

// calculating method have been taken from td eclipse (page 896)
// calc transmissibility between fully adjacent cells index1 and index2

// FIXME: ntg_ and etc arrays
t_double 
mesh_grdecl::calc_tran(const t_long ext_index1, const t_long ext_index2, const plane_t &plane1, 
                                       const fpoint3d_t &center1, const fpoint3d_t &center2, direction d_dir, plane_t* plane2) const
{
  t_double tran;
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
  


 t_double koef1, koef2 ; //koef = (A,Di)/(Di,Di)
  
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

  t_double Ti, Tj;

  // FIXME: ntg_array and etc access
  t_double ntg_index1 = 1;
  t_double ntg_index2 = 1;
  if (ntg_array)
    {
      ntg_index1 = ntg_array->data ()[ext_index1];
      ntg_index2 = ntg_array->data ()[ext_index2];
    }

  if (d_dir == along_dim1) //lengthwise OX
    {
      Ti = permx_array->data ()[ext_index1]*ntg_index1*koef1;
      Tj = permx_array->data ()[ext_index2]*ntg_index2*koef2;
      tran = darcy_constant / (1 / Ti + 1 / Tj);
      if (multx_array)
        tran *= multx_array->data ()[ext_index1];
    }
  else if (d_dir == along_dim2) //lengthwise OY
    {
      Ti = permy_array->data ()[ext_index1]*ntg_index1*koef1;
      Tj = permy_array->data ()[ext_index2]*ntg_index2*koef2;
      tran = darcy_constant / (1 / Ti + 1 / Tj);
      if (multy_array)
        tran *= multy_array->data ()[ext_index1];
    }
  else //lengthwise OZ
    {
      Ti = permz_array->data ()[ext_index1]*koef1;
      Tj = permz_array->data ()[ext_index2]*koef2;
      tran = darcy_constant / (1 / Ti + 1 / Tj);
      if (multz_array)
        tran *= multz_array->data ()[ext_index1];
    }
  
  /*
  #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
    if (ext_index1 < 1000)
      fprintf (fp, "[Ti: %lf; Tj: %lf]", Ti, Tj);
  #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
  */
  return tran;
}


t_double 
mesh_grdecl::calc_tran_boundary (const t_long ext_index1, const plane_t &plane1, const fpoint3d_t &center1, direction d_dir) const
{
  t_double tran;
  fpoint3d_t D1, A, plane_contact_center;
  
  get_plane_center(plane1, plane_contact_center);
  plane_contact_center.distance_to_point (center1, D1);
  A = get_projection_on_all_axis_for_one_side(plane1);
  
  t_double koef1; //koef = (A,Di)/(Di,Di)
  
  koef1 = A*D1 / (D1*D1);
  
  t_double Ti;

  // FIXME: ntg_array and etc access
  t_double ntg_index1 = 1;
  if (ntg_array)
    {
      ntg_index1 = ntg_array->data ()[ext_index1];
    }

  if (d_dir == along_dim1) //lengthwise OX
    {
      Ti = permx_array->data ()[ext_index1]*ntg_index1*koef1;
      tran = darcy_constant / (2 / Ti);
      if (multx_array)
        tran *= multx_array->data ()[ext_index1];
    }
  else if (d_dir == along_dim2) //lengthwise OY
    {
      Ti = permy_array->data ()[ext_index1]*ntg_index1*koef1;
      tran = darcy_constant / (2 / Ti);
      if (multy_array)
        tran *= multy_array->data ()[ext_index1];
    }
  else //lengthwise OZ
    {
      Ti = permz_array->data ()[ext_index1]*koef1;
      tran = darcy_constant / (2 / Ti);
      if (multz_array)
        tran *= multz_array->data ()[ext_index1];
    }
  
  return tran;
}




int mesh_grdecl::build_jacobian_and_flux_connections (const sp_bcsr_t jacobian,const sp_flux_conn_iface_t flux_conn,
                                                                  spv_long boundary_array)

{
  return build_jacobian_and_flux_connections_add_boundary (jacobian, flux_conn, boundary_array);
}


void mesh_grdecl::get_block_dx_dy_dz (t_long n_elem, t_double &dx, t_double &dy, t_double &dz) const
  {
    element_t elem = calc_element(n_elem);
    point3d_t sizes = elem.get_dx_dy_dz (); 
    dx = sizes[0];
    dy = sizes[1];
    dz = sizes[2];
  }
  
spv_double mesh_grdecl::get_element_sizes (const t_long n_element) const
  {
      double dx, dy, dz;
      get_block_dx_dy_dz(n_element, dx, dy, dz);
      spv_double sizes = BS_KERNEL.create_object(v_double::bs_type());
      sizes->resize(3);
      (*sizes)[0] = dx;
      (*sizes)[1] = dy;
      (*sizes)[2] = dz;
      return sizes;
  }


t_double mesh_grdecl:: get_block_dx(t_long n_elem) const
  {
    return calc_element(n_elem).get_dx();
  }


t_double mesh_grdecl:: get_block_dy(t_long n_elem) const
  {
    return calc_element(n_elem).get_dy();
  }


t_double mesh_grdecl:: get_block_dz(t_long n_elem) const
  {
    return calc_element(n_elem).get_dz();
  }

t_double mesh_grdecl:: get_block_dz_ext(t_long i, t_long j, t_long k) const
  {
    element_t element;
    calc_element(i, j, k, element);
    return element.get_dz();
  }

t_double  mesh_grdecl::get_dtop(t_long n_elem) const
{
  element_t elem;
  
  elem = calc_element (n_elem);

  return elem.get_center().z - elem.get_dz();
}


boost::python::list mesh_grdecl::calc_element_tops ()
{
  element_t element;
  spv_float tops, prop;
  spv_long indexes;
  t_long i, j, k, c, ind, *indexes_data;
  t_float *tops_data, *prop_data;
  boost::python::list myavi_list;

  tops = give_kernel::Instance().create_object(v_float::bs_type());
  indexes = give_kernel::Instance().create_object(v_long::bs_type());
  prop = give_kernel::Instance().create_object(v_float::bs_type());

  tops->resize (n_elements * 8 * 3);
  indexes->resize (n_elements * 8);
  prop->resize (n_elements);

  indexes_data = &(*indexes)[0];
  tops_data = &(*tops)[0];
  prop_data = &(*prop)[0];

  ind = 0;
   
  t_float const *poro = poro_array->data ();
  for (i = 0; i < nx; ++i)
      for (j = 0; j < ny; ++j)
          for (k = 0; k < nz; ++k, ++ind)
            {
            
                /*
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
    */
              calc_element (i, j, k, element);
              prop_data[ind] = poro[ind];
              for (c = 0; c < 8; ++c)
                {
                  tops_data[8 * 3 * ind + 3 * c] = element.get_corners()[c].x;
                  tops_data[8 * 3 * ind + 3 * c + 1] = element.get_corners()[c].y;
                  tops_data[8 * 3 * ind + 3 * c + 2] = element.get_corners()[c].z * 10;
                }

                
              indexes_data[8 * ind] = 8 * ind;
              indexes_data[8 * ind + 1] = 8 * ind + 1;
              indexes_data[8 * ind + 2] = 8 * ind + 3;
              indexes_data[8 * ind + 3] = 8 * ind + 2;

              indexes_data[8 * ind + 4] = 8 * ind + 4;
              indexes_data[8 * ind + 5] = 8 * ind + 5;
              indexes_data[8 * ind + 6] = 8 * ind + 7;
              indexes_data[8 * ind + 7] = 8 * ind + 6;

/*            
              // planes indexes
              indexes_data[24 * ind] = 8 * ind;
              indexes_data[24 * ind + 1] = 8 * ind + 1;
              indexes_data[24 * ind + 2] = 8 * ind + 3;
              indexes_data[24 * ind + 3] = 8 * ind + 2;

              indexes_data[24 * ind + 4] = 8 * ind + 4;
              indexes_data[24 * ind + 5] = 8 * ind + 5;
              indexes_data[24 * ind + 6] = 8 * ind + 7;
              indexes_data[24 * ind + 7] = 8 * ind + 6;

              indexes_data[24 * ind + 8] = 8 * ind;
              indexes_data[24 * ind + 9] = 8 * ind + 2;
              indexes_data[24 * ind + 10] = 8 * ind + 6;
              indexes_data[24 * ind + 11] = 8 * ind + 4;

              indexes_data[24 * ind + 12] = 8 * ind + 1;
              indexes_data[24 * ind + 13] = 8 * ind + 3;
              indexes_data[24 * ind + 14] = 8 * ind + 7;
              indexes_data[24 * ind + 15] = 8 * ind + 5;

              indexes_data[24 * ind + 16] = 8 * ind ;
              indexes_data[24 * ind + 17] = 8 * ind + 1;
              indexes_data[24 * ind + 18] = 8 * ind + 5;
              indexes_data[24 * ind + 19] = 8 * ind + 4;

              indexes_data[24 * ind + 20] = 8 * ind + 2;
              indexes_data[24 * ind + 21] = 8 * ind + 3;
              indexes_data[24 * ind + 22] = 8 * ind + 7;
              indexes_data[24 * ind + 23] = 8 * ind + 6;
*/
            }

  myavi_list.append(tops);
  myavi_list.append(indexes);
  myavi_list.append(prop);

  return myavi_list;
}

spv_float mesh_grdecl::calc_cells_vertices() {
  spv_float tops = give_kernel::Instance().create_object(v_float::bs_type());
  tops->resize (n_elements * 8 * 3);

  t_float* tops_data = &(*tops)[0];
  t_long ind = 0;

  element_t element;
  for (t_long i = 0; i < nx; ++i)
      for (t_long j = 0; j < ny; ++j)
          for (t_long k = 0; k < nz; ++k, ++ind)
            {
              calc_element (i, j, k, element);
              for (t_long c = 0; c < 8; ++c)
                {
                  const t_long offs = 8 * 3 * ind + 3 * c;
                  const element_t::fpoint3d_t& cur_corner = element.get_corners()[c];
                  tops_data[offs]     = cur_corner.x;
                  tops_data[offs + 1] = cur_corner.y;
                  tops_data[offs + 2] = cur_corner.z;
                }
            }

  return tops;
}

spv_float mesh_grdecl::calc_cells_vertices_xyz() {
  spv_float tops = give_kernel::Instance().create_object(v_float::bs_type());
  tops->resize (n_elements * 8 * 3);

  t_float* tops_data = &(*tops)[0];
  t_long ind = 0;

  element_t element;
  for (t_long k = 0; k < nz; ++k)
    for (t_long j = 0; j < ny; ++j)
      for (t_long i = 0; i < nx; ++i, ++ind)
          {
              calc_element (i, j, k, element);
              for (t_long c = 0; c < 8; ++c)
                {
                  const t_long offs = 8 * 3 * ind + 3 * c;
                  const element_t::fpoint3d_t& cur_corner = element.get_corners()[c];
                  tops_data[offs]     = cur_corner.x;
                  tops_data[offs + 1] = cur_corner.y;
                  tops_data[offs + 2] = cur_corner.z;
                }
            }

  return tops;
}

boost::python::list mesh_grdecl::calc_element_center ()
{
  element_t element;
  spv_float centers, prop;
  t_long i, j, k /*, c */, ind /*, *indexes_data */;
  t_float *centers_data, *prop_data;
  boost::python::list myavi_list;

  centers = give_kernel::Instance().create_object(v_float::bs_type());
  prop = give_kernel::Instance().create_object(v_float::bs_type());

  centers->resize (n_elements * 3);
  prop->resize (n_elements);

  centers_data = &(*centers)[0];
  prop_data = &(*prop)[0];

  ind = 0;
   
  for (i = 0; i < nx; ++i)
      for (j = 0; j < ny; ++j)
          for (k = 0; k < nz; ++k, ++ind)
            {
              calc_element (i, j, k, element);
              centers_data[3 * ind] = element.get_center().x;
              centers_data[3 * ind + 1] = element.get_center().y;
              centers_data[3 * ind + 2] = element.get_center().z;
            }

  myavi_list.append(centers);
  myavi_list.append(prop);

  return myavi_list;
}



void mesh_grdecl::generate_array()
{
#if 0
  t_long n_size = n_elements;
  poro_array->clear();
  ntg_array->clear();

  permx_array->clear();
  permy_array->clear();
  permz_array->clear();

  multx_array->clear();
  multy_array->clear();
  multz_array->clear();


  for (t_long i =0; i < n_size; ++i)
    {
      poro_array->push_back(0.2f);
      ntg_array->push_back(0.4f);

      permx_array->push_back(10.0f);
      permy_array->push_back(10.0f);
      permz_array->push_back(10.0f);

      multx_array->push_back(1.0f);
      multy_array->push_back(1.0f);
      multz_array->push_back(1.0f);
    }
#endif
}


template <typename loop_t>
struct build_jacobian_rows_class
{
  typedef mesh_grdecl                 mesh_t;
  typedef  mesh_t::plane_t                plane_t;
  typedef  mesh_t::element_zcorn_t_long  element_zcorn_t_long;

  build_jacobian_rows_class (mesh_grdecl  *mesh, loop_t *loop, std::set <t_long, std::less <t_long> > &boundary_set, t_long *rows_ptr)
  : mesh (mesh)
  , loop (loop)
  , boundary_set (boundary_set)
  , rows_ptr (rows_ptr)
  , nx (mesh->get_nx())
  , ny (mesh->get_ny())
  , nz (mesh->get_nz())
  {
  }

  void
  prepare (t_long i, t_long j, t_long k)
  {
    ext_index  = k + (nz * j) + (i * ny * nz);
  }

  // check if (i, j) column of cells is adjacent to neigbours
  // that is true, if every cell of a column is fully adjacent to X and Y neighbour
  
  bool
  check_column_adjacency (t_long i, t_long j)
  {
    bool flag = true;
    t_long k, zindex, zindex1, index;
         
    t_int const *actnum = mesh->actnum_array->data ();
    t_float const *zcorn = mesh->pinner_->zcorn_->data ();
    for (k = 0; k < nz; ++k)
      {
        index = k + (nz * j) + (i * ny * nz);
         
        // miss inactive blocks
        if (!actnum[index])
          {
            continue;
          }
          
        // check blocks adjacency
        zindex = i * 2 + j * 4 * nx + k * 8 *nx * ny;
        zindex1 = zindex + 2; // check next by x
       
        if (i + 1 == nx ||
            (zcorn[zindex + 1] == zcorn[zindex1] &&
            zcorn[zindex +  2 * nx + 1] == zcorn[zindex1 +  2 * nx] &&
            zcorn[zindex + 4 * nx * ny + 1] == zcorn[zindex1 + 4 * nx * ny ] &&
            zcorn[zindex + 4 * nx * ny + 2 * nx + 1] == zcorn[zindex1 + 4 * nx * ny + 2 * nx]))
          {
            
            zindex1 = zindex + 4 * nx; // check next by y
            if (j + 1 == ny ||
                (zcorn[zindex + 2 * nx] == zcorn[zindex1] &&
                zcorn[zindex + 2 * nx + 1] == zcorn[zindex1 +  1] &&
                zcorn[zindex + 4 * nx * ny + 2 * nx] == zcorn[zindex1 + 4 * nx * ny] &&
                zcorn[zindex + 4 * nx * ny + 2 * nx + 1] == zcorn[zindex1+ 4 * nx * ny + 1]))
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
  change_by_x (t_long /*i*/, t_long /*j*/, t_long /*k*/, t_long ext_index2, bool /*is_adjacent*/)
  {
    rows_ptr[mesh->convert_ext_to_int (ext_index) + 1]++;
    rows_ptr[mesh->convert_ext_to_int (ext_index2) + 1]++;
  }

  void
  change_by_y (t_long /*i*/, t_long /*j*/, t_long /*k*/, t_long ext_index2, bool /*is_adjacent*/)
  {
    rows_ptr[mesh->convert_ext_to_int (ext_index) + 1]++;
    rows_ptr[mesh->convert_ext_to_int (ext_index2) + 1]++;
  }

  void
  change_by_z (t_long /*i*/, t_long /*j*/, t_long /*k*/, t_long ext_index2, bool /*is_adjacent*/)
  {
    rows_ptr[mesh->convert_ext_to_int (ext_index) + 1]++;
    rows_ptr[mesh->convert_ext_to_int (ext_index2) + 1]++;
  }

  void
  add_boundary (t_long external_cell_index)
  {
    boundary_set.insert (external_cell_index);
  }

  mesh_grdecl                   *mesh;
  loop_t                                    *loop;
  std::set <t_long, std::less <t_long> >  &boundary_set;
  t_long                                   *rows_ptr;
  t_long                                   nx;
  t_long                                   ny;
  t_long                                   nz;
  t_long                                   ext_index;
};

template <typename L, typename BS, typename RP>
build_jacobian_rows_class <L>
build_jacobian_rows (mesh_grdecl *mesh, L *l, BS &bs, RP *rp)
{
  return build_jacobian_rows_class <L> (mesh, l, bs, rp);
}

template <typename loop_t>
struct build_jacobian_cols_class
{
  typedef mesh_grdecl                 mesh_t;
  typedef  mesh_t::element_t              element_t;
  typedef  mesh_t::plane_t                plane_t;
  typedef  mesh_t::element_zcorn_t_long  element_zcorn_t_long;

  build_jacobian_cols_class (mesh_t *mesh, loop_t *loop, t_long *rows_ptr, t_long *cols_ind,
    t_long *cols_ind_transmis, t_float *values_transmis,
    t_long *matrix_block_idx_minus, t_long *matrix_block_idx_plus)
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
    t_long i;
    t_long n = mesh->n_active_elements;
    
    
    rows_ptr_tmp.resize (n);
    for (i = 0; i < n; ++i)
      {
        // let it point to first non-diagonal column in each row
        rows_ptr_tmp[i] = rows_ptr[i] + 1;
        cols_ind[rows_ptr[i]]= i;
      }
  }
  
  void
  prepare (t_long i, t_long j, t_long k)
  {
    ext_index  = k + (nz * j) + (i * ny * nz);
    int_index  = mesh->convert_ext_to_int (ext_index);

    mesh->calc_element(i, j, k, element);
    center1  = element.get_center();
  }
  
    
  void
  change_jac_and_flux_conn( const t_long /*ext_index1*/, const t_long ext_index2, t_double tran)
  {
    t_long index1 = mesh->convert_ext_to_int (ext_index);
    t_long index2 = mesh->convert_ext_to_int (ext_index2);
    
    cols_ind[rows_ptr_tmp[index1]] = index2;
    cols_ind[rows_ptr_tmp[index2]] = index1;

    matrix_block_idx_minus[loop->con_num] = rows_ptr[index1];
    matrix_block_idx_minus[loop->con_num + 1] = rows_ptr_tmp[index1];
    matrix_block_idx_plus[loop->con_num + 1] = rows_ptr[index2];
    matrix_block_idx_plus[loop->con_num] = rows_ptr_tmp[index2];
    
    ++rows_ptr_tmp[index1];
    ++rows_ptr_tmp[index2];
    
    cols_ind_transmis[loop->con_num] = index1;
    cols_ind_transmis[loop->con_num + 1] = index2;
    
    values_transmis[loop->con_num] = tran;
    values_transmis[loop->con_num + 1] = -tran;
    loop->con_num += 2;
  }

  bool
  check_column_adjacency (t_long i, t_long j)
  {
    return loop->is_column_adjacent[i + j * nx];
  }
  
  void
  add_boundary (t_long)
  {
  }

  void
  change_by_x (t_long i, t_long j, t_long k, t_long ext_index2, bool is_adjacent)
  {
    t_double tran;
    plane_t plane1;
    element_t element2;
    
    mesh->calc_element(i, j, k, element2);
    fpoint3d center2  = element2.get_center ();
    
    element.get_plane (x_axis_plus, plane1);
    
     
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
    
    // FIXME: uncomment
    //if (tran != 0)
      {
        change_jac_and_flux_conn (ext_index, ext_index2, tran);
      }  
    //else
    //  {
    //    BS_ASSERT (false) (tran) (i) (j) (k) (ext_index) (ext_index2);
    //  }
   
    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index < 1000)
          fprintf (fp, " %d [%d;%d;%d] (%lf)", ext_index2, i, j, k, tran);
    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX 
  }

  void
  change_by_y (t_long i, t_long j, t_long k, t_long ext_index2, bool is_adjacent)
  {
    t_double tran;
    
    plane_t plane1;
    element_t element2;
    
    mesh->calc_element(i, j, k, element2);
    fpoint3d center2  = element2.get_center ();
    element.get_plane (y_axis_plus, plane1);
     
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
     
    // FIXME: uncomment
    //if (tran != 0)
      {
        change_jac_and_flux_conn (ext_index, ext_index2, tran);
      }
    //else
    //  {
    //    BS_ASSERT (false) (tran) (i) (j) (k) (ext_index) (ext_index2);
    //  }

    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index < 1000)
          fprintf (fp, " %d [%d;%d;%d] (%lf)", ext_index2, i, j, k, tran);
    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX 
  }

  void
  change_by_z (t_long i, t_long j, t_long k, t_long ext_index2, bool /*is_adjacent*/)
  {
    t_double tran;
    
    plane_t plane1;
    element_t element2;
    
    mesh->calc_element(i, j, k, element2);
    fpoint3d center2  = element2.get_center ();
    element.get_plane (z_axis_plus, plane1);
    
    tran = mesh->calc_tran (ext_index, ext_index2, plane1, center1, center2, along_dim3);
         
    // FIXME: uncomment
    //if (tran != 0)
      {
        change_jac_and_flux_conn (ext_index, ext_index2, tran);
      }
    //else
    //  {
    //    BS_ASSERT (false) (tran) (i) (j) (k) (ext_index) (ext_index2);
    //  }

    #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index < 1000)
          fprintf (fp, " %d [%d;%d;%d] (%lf)", ext_index2, i, j, k, tran);
    #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
  }

  mesh_t              *mesh;
  loop_t              *loop;
  t_long             *rows_ptr;
  t_long             *cols_ind;
  t_long             *cols_ind_transmis;
  t_float   *values_transmis;
  t_long             *matrix_block_idx_minus;
  t_long             *matrix_block_idx_plus;

  t_long             nx;
  t_long             ny;
  t_long             nz;

  // points to current column position for each row
  // while adding new columns
  stdv_long       rows_ptr_tmp; 

  t_long             ext_index; // index1
  t_long             int_index; // index_ijk

  element_t           element;
  fpoint3d            center1;

};

template <typename L, typename RP, typename CI, typename CT, typename FC>
build_jacobian_cols_class <L>
build_jacobian_cols (mesh_grdecl *m, L *l, RP *rp, CI *ci, CT &conn_trans, FC &flux_conn)
{
  return build_jacobian_cols_class <L> (m, l, rp, ci,
    &(*conn_trans->get_cols_ind ())[0], &(*conn_trans->get_values ())[0],
    &(*flux_conn->get_matrix_block_idx_minus ())[0], &(*flux_conn->get_matrix_block_idx_plus ())[0]);
}


struct build_jacobian_and_flux : boost::noncopyable
{
  typedef mesh_grdecl                 mesh_t;
  typedef  mesh_t::plane_t                plane_t;
  typedef  mesh_t::element_zcorn_t_long  element_zcorn_t_long;

  build_jacobian_and_flux (mesh_grdecl  *mesh)
  : mesh (mesh)
  , nx (mesh->get_nx())
  , ny (mesh->get_ny())
  , nz (mesh->get_nz())
  , con_num (0)
  {
    is_column_adjacent.assign (nx * ny, true);
  }

  template <typename loop_body_t>
  void
  adjacent_columns (loop_body_t &loop_body, t_int const *actnum, t_long i, t_long j, t_long &n_adj_elems)
  {
    for (t_long k = 0; k < nz; ++k)
      {
        t_long ext_index1  = k + (nz * j) + (i * ny * nz);
        
        //skip non-active cells
        if (!actnum[ext_index1])
          continue;
          
        n_adj_elems++;
        
        #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index1 < 1000)
          fprintf (fp, "\n%d [%d;%d;%d]", ext_index1, i, j, k);
        #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
          
        loop_body.prepare (i, j, k);
        
        if (i+1 < nx && actnum[ext_index1 + 1])
          {
            loop_body.change_by_x (i + 1, j, k, ext_index1 + 1, true);
          }
          
        if (j+1 < ny && actnum[ext_index1 + nx])
          {
            loop_body.change_by_y (i, j + 1, k, ext_index1 + nx, true);
          }
         
        if (k+1 < nz && actnum[ext_index1 + nx * ny])
          {
            loop_body.change_by_z (i, j, k + 1, ext_index1 + nx * ny, true);
          }
      }
  }

  template <typename loop_body_t>
  void
  non_adjacent_columns (loop_body_t &loop_body, t_int const *actnum, t_float *zcorn_array, t_long i, t_long j, t_long &n_non_adj_elems)
  {
    t_long last_k_x = 0;
    t_long last_k_y = 0;
    
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
                
    for (t_long k = 0; k < nz; ++k)
      {
        t_long ext_index1  = k + (nz * j) + (i * ny * nz);
        
        //skip non-active cells
        if (!actnum[ext_index1])
          continue;
        
        n_non_adj_elems++;
        
        #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX   
        if (ext_index1 < 1000)
          fprintf (fp, "\n%d [%d;%d;%d]", ext_index1, i, j, k);
        #endif //BS_MESH_WRITE_TRANSMISS_MATRIX  
        
        element_zcorn_t_long zcorn_index1;
        mesh->get_element_zcorn_index (i, j, k, zcorn_index1);
        
        loop_body.prepare (i, j, k);
        
        // if X neighbour exists and current element`s X+ plane is not a line
        if (i + 1 < nx && ((zcorn_array[zcorn_index1[1]] != zcorn_array[zcorn_index1[5]]) || (zcorn_array[zcorn_index1[3]] != zcorn_array[zcorn_index1[7]])))
          {
            t_long k_x = last_k_x - 1;
            
            element_zcorn_t_long zcorn_index2;
            // search first possible neighbour
            do
              {
                k_x++;
                mesh->get_element_zcorn_index (i + 1, j, k_x, zcorn_index2);
              }
            while ((k_x < nz) && ((zcorn_array[zcorn_index1[1]] >= zcorn_array[zcorn_index2[4]]) && (zcorn_array[zcorn_index1[3]] >= zcorn_array[zcorn_index2[6]])));
            
            // calc all neihbours
            while ((k_x < nz) && ((zcorn_array[zcorn_index1[5]] > zcorn_array[zcorn_index2[0]]) || (zcorn_array[zcorn_index1[7]] > zcorn_array[zcorn_index2[2]])))
              {
                if ((zcorn_array[zcorn_index1[5]] >= zcorn_array[zcorn_index2[4]]) && (zcorn_array[zcorn_index1[7]] >= zcorn_array[zcorn_index2[6]]))
                  {
                    // this (i + 1, j, k_x) neigbour won`t touch next (i, j, k + 1) element
                    last_k_x = k_x + 1;
                  }
                  
                t_long ext_index2 = ext_index1 + (k_x - k) + ny * nz;
                
                // if neighbour active and it`s X- plane is not a line
                if (actnum[ext_index2] && ((zcorn_array[zcorn_index2[0]] != zcorn_array[zcorn_index2[4]]) || (zcorn_array[zcorn_index2[2]] != zcorn_array[zcorn_index2[6]])))
                  {
                    loop_body.change_by_x (i + 1, j, k_x, ext_index2, false);
                  }
                k_x++;
                mesh->get_element_zcorn_index (i + 1, j, k_x, zcorn_index2);
              }
          }
          
        // if Y neighbour exists and current element`s Y+ plane is not a line
        if (j + 1 < ny && ((zcorn_array[zcorn_index1[2]] != zcorn_array[zcorn_index1[6]]) || (zcorn_array[zcorn_index1[3]] != zcorn_array[zcorn_index1[7]])))
          {
            t_long k_y = last_k_y - 1;
            
            element_zcorn_t_long zcorn_index2;
            // search first possible neighbour
            do
              {
                k_y++;
                mesh->get_element_zcorn_index (i, j + 1, k_y, zcorn_index2);
              }
            while ((k_y < nz) && ((zcorn_array[zcorn_index1[2]] >= zcorn_array[zcorn_index2[4]]) && (zcorn_array[zcorn_index1[3]] >= zcorn_array[zcorn_index2[5]])));
            
            // calc all neighbours
            while ((k_y < nz) && ((zcorn_array[zcorn_index1[6]] > zcorn_array[zcorn_index2[0]]) || (zcorn_array[zcorn_index1[7]] > zcorn_array[zcorn_index2[1]])))
              {
                if ((zcorn_array[zcorn_index1[6]] >= zcorn_array[zcorn_index2[4]]) && (zcorn_array[zcorn_index1[7]] >= zcorn_array[zcorn_index2[5]]))
                  {
                    // this (i, j + 1, k_y) neigbour won`t touch next (i, j, k + 1) element
                    last_k_y = k_y + 1;
                  }
                  
                t_long ext_index2 = ext_index1 + (k_y - k) + nz;
                
                // if neighbour active and it`s Y- plane is not a line
                if (actnum[ext_index2] && ((zcorn_array[zcorn_index2[0]] != zcorn_array[zcorn_index2[4]]) || (zcorn_array[zcorn_index2[1]] != zcorn_array[zcorn_index2[5]])))
                  {
                    loop_body.change_by_y (i, j + 1, k_y, ext_index2, false);
                  }
                k_y++;
                mesh->get_element_zcorn_index (i, j + 1, k_y, zcorn_index2);
              }
          }
          
        if (k + 1 < nz && actnum[ext_index1 + nx * ny])
          {
            loop_body.change_by_z (i, j, k + 1, ext_index1 + 1, true);
          }    
      }
  }

  template <typename loop_body_t>
  void
  cell_loop (loop_body_t loop_body)
  {
    t_long n_adj_elems = 0, n_non_adj_elems = 0;
    
    t_float *zcorn_array = mesh->zcorn_array;
    t_int const *actnum = mesh->actnum_array->data ();

    for (t_long i = 0; i < nx; ++i)
      {
        for (t_long j = 0; j < ny; ++j)
          {
            bool is_adjacent = loop_body.check_column_adjacency (i, j);
            
            if (is_adjacent)
              {   
                adjacent_columns (loop_body, actnum, i, j, n_adj_elems);
              }
            else
              {
                non_adjacent_columns (loop_body, actnum, zcorn_array, i, j, n_non_adj_elems);
              }
          }
      }
    BOSWARN (section::mesh, level::warning)<< boost::format ("MESH_GRDECL: elements adj %d, non-adj %d, total %d") \
            % n_adj_elems % n_non_adj_elems % (n_adj_elems + n_non_adj_elems) << bs_end;  
    BOSWARN (section::mesh, level::warning)<< boost::format ("MESH_GRDECL: number of tran calcs is %d") % n_tran_calc << bs_end;  
  }

  mesh_grdecl     *mesh;
  t_long                     nx;
  t_long                     ny;
  t_long                     nz;
  shared_vector <bool>        is_column_adjacent;
  t_long                     con_num;
};



int mesh_grdecl::build_jacobian_and_flux_connections_add_boundary (const sp_bcsr_t jacobian,
                                                                   const sp_flux_conn_iface_t flux_conn,
                                                                   spv_long /*boundary_array*/)
{
  write_time_to_log init_time ("Mesh transmissibility calculation", ""); 
  
  
  t_long i, n_non_zeros;
  sp_bcsr_t conn_trans;
  
  n_connections = 0;
  jacobian->get_cols_ind()->clear();
  jacobian->get_rows_ptr()->clear();
  jacobian->init_struct(n_active_elements, n_active_elements, n_active_elements);
  
  t_long *rows_ptr = jacobian->get_rows_ptr()->data ();
  // FIXME: check size of rows_ptr
  rows_ptr[0] = 0;

  std::vector<bool> is_butting(nx*ny,false);

  std::set<t_long, std::less<t_long> > boundary_set;

  build_jacobian_and_flux  build_jacobian (this);
  
 #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX     
  fp = fopen ("transmiss.out", "w");
 #endif //BS_MESH_WRITE_TRANSMISS_MATRIX   

  //first step - define and fill - rows_ptr (jacobian)
  build_jacobian.cell_loop (build_jacobian_rows (this, &build_jacobian, boundary_set, rows_ptr));

  //////jacobian//////////////////////
  //sum rows_ptr
  for (i = 0; i < n_active_elements; ++i)
    {
      rows_ptr[i + 1] += rows_ptr[i] + 1;
    }
    
  n_non_zeros = rows_ptr[n_active_elements];
  n_connections = (n_non_zeros - n_active_elements) / 2;
    
  //create cols_ind
  jacobian->get_cols_ind()->resize(n_non_zeros);
  t_long *cols_ind = jacobian->get_cols_ind()->data ();


  ////////transmis/////////////////////////
  conn_trans = flux_conn->get_conn_trans();
  conn_trans->init (n_connections, 2 * n_connections, 1, 2 * n_connections);

  flux_conn->get_matrix_block_idx_minus ()->resize(n_connections * 2);
  flux_conn->get_matrix_block_idx_plus ()->resize(n_connections * 2);
  
  if (!n_connections)
    {
      for (i = 0; i < n_non_zeros; ++i)
        cols_ind[i] = i;
      return 0;
    }

  t_long *rows_ptr_transmis = conn_trans->get_rows_ptr()->data ();
  for (i = 0; i < n_connections + 1; ++i)
    rows_ptr_transmis[i] = i * 2;

  //second step - fill and define cols_ind
  build_jacobian.con_num = 0;
  build_jacobian.cell_loop (build_jacobian_cols (this, &build_jacobian, rows_ptr, cols_ind, conn_trans, flux_conn));


  //boundary_array.assign(boundary_set.begin(), boundary_set.end());
  
  #ifdef BS_MESH_WRITE_TRANSMISS_MATRIX  
    fflush (fp);
    fclose (fp);
  #endif //BS_MESH_WRITE_TRANSMISS_MATRIX 

  //return (int) boundary_array->size();
  return 0;
}


int mesh_grdecl::intersect_trajectories ()
{
  int n_range_cells = 1000;
  int wrong_cells = 0;

  int n_x_ranges = std::ceil (nx / std::pow(float(n_range_cells), float(1/3)));
  int n_y_ranges = std::ceil (ny / std::pow(float(n_range_cells), float(1/3)));
  int n_z_ranges = std::ceil (nz / std::pow(float(n_range_cells), float(1/3)));

  element_t element;
  t_int *actnum = actnum_array->data ();

  for (t_long i = 0; i < nx; ++i)
    for (t_long j = 0; j < ny; ++j)
      {
        for (t_long k = 0; k < nz; ++k)
          {
            t_long index = k + j * nz + i * ny * nz;

            // miss inactive blocks
            if (actnum[index])
              {
                calc_element (i, j, k, element);
                mesh_element3d::corners_t corns = element.get_corners();
                
                // check X
                if (((corns[1].x - corns[0].x) * (corns[3].x - corns[2].x) < 0) ||
                    ((corns[1].x - corns[0].x) * (corns[5].x - corns[4].x) < 0) ||
                    ((corns[1].x - corns[0].x) * (corns[7].x - corns[6].x) < 0))
                  {
                    actnum[index] = 0;
                    wrong_cells ++;
                    continue;
                  }

                // check Y
                if (((corns[2].y - corns[0].y) * (corns[3].y - corns[1].y) < 0) ||
                    ((corns[2].y - corns[0].y) * (corns[7].y - corns[5].y) < 0) ||
                    ((corns[2].y - corns[0].y) * (corns[6].y - corns[4].y) < 0))
                  {
                    actnum[index] = 0;
                    wrong_cells ++;
                    continue;
                  }
                
                
                // check Z
                if (((corns[4].z - corns[0].z) * (corns[5].z - corns[1].z) < 0) ||
                    ((corns[4].z - corns[0].z) * (corns[6].z - corns[2].z) < 0) ||
                    ((corns[4].z - corns[0].z) * (corns[7].z - corns[3].z) < 0))
                  {
                    actnum[index] = 0;
                    wrong_cells ++;
                    continue;
                  }
              }
        }
    }

  




  return 0;
}

#ifdef _HDF5_MY

int mesh_grdecl::create_array_hdf5(const char *dataset_name, H5::H5File &file_hdf5, H5::DataSet **dataset)
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


bool mesh_grdecl::file_open_activs_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt)
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


bool mesh_grdecl::file_open_cube_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt)
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
      vector<t_double> buf_array(6 * (it - is + 2) * (jt - js  + 2));
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
      hsize_t dims_memory_zcorn_array[] = {8 * (it - is + 1) * (jt - js + 1)};
      H5::DataSpace dataspace_memory_zcorn(1, dims_memory_zcorn);
      H5::DataSpace dataspace_file_zcorn = dataset_zcorn.getSpace();
      zcorn_array.resize(8 * (it - is + 1) * (jt - js  + 1) * (kt - ks + 1));
      // hyperslab settings (common for each plane)
      hsize_t count_zcorn_array[] = {4 * (jt - js + 1)};
      hsize_t stride_zcorn_array[] = {2 * nx};
      hsize_t block_zcorn_array[] = {2 * (it - is + 1)};
      // we read "zcorn" array by planes - kt-ks+1 read operations
      for (int k = 0; k < kt - ks + 1; k++)
        {
          // determine start array for each plane individually
          hsize_t start_zcorn_array[] = {2 * is + 4 * js * nx + 8 * k * nx * ny};
          dataspace_file_zcorn.selectHyperslab(H5S_SELECT_SET, count_zcorn, start_zcorn, stride_zcorn, block_zcorn);
          dataset_zcorn.read(&zcorn_array[k * count_zcorn_array[0] * block_zcorn_array[0]], H5::PredType::NATIVE_FLOAT, dataspace_memory_zcorn, dataspace_file_zcorn);
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


int mesh_grdecl::append_array_hdf5(const t_double *arr, size_t arr_length, H5::DataSet *dataset)
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


bool mesh_grdecl::file_open_cube_with_hdf5_swap(const char* file_name)
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


bool mesh_grdecl::file_open_actnum(const char* /*file_name*/)
{
#if 0
  fstream file(file_name,  ios::in);
  char buf[255];
  char *start_ptr,*end_ptr;
  if (!file.is_open())
    {
      n_active_elements =  accumulate(actnum_array.begin(),actnum_array.end(),0);
      return false;
    }
  actnum_array.clear();

  while (!file.eof())
    {
      file >> buf;
      start_ptr = buf;
      if (strcmp (buf, "/") != 0)
        actnum_array.push_back((int)strtod (start_ptr, &end_ptr));
      else
        break;
    }
  n_active_elements =  accumulate(actnum_array.begin(),actnum_array.end(),0);
#endif
  return true;
}



bool mesh_grdecl::file_open_cube(const char* /*file_name*/)
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

  t_long i_count = 0;
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


//BS_INST_STRAT(mesh_grdecl);
