/*!
	\file mesh_element3d.cpp
	\brief This file defines mesh element in 3 dimensions
	\author Mark Khait
	\date 2009-09-23
 */

#include "bs_mesh_stdafx.h"
#include "mesh_element3d.h"


element_plane_orientation_t get_reverse_orientation (element_plane_orientation_t orient)
  {
    switch (orient)
      {
        case x_axis_plus:   return x_axis_minus;
        case x_axis_minus:  return x_axis_plus;
        case y_axis_plus:   return y_axis_minus;
        case y_axis_minus:  return y_axis_plus;
        case z_axis_plus:   return z_axis_minus;
        case z_axis_minus:  return z_axis_plus;
        default: bs_throw_exception ("Invalid orientation!");
      }
  }

template<class strategy_t>
mesh_element3d<strategy_t>::mesh_element3d()
{
  n_corners = N_ELEMENT_CORNERS;
  n_plane_corners = N_PLANE_CORNERS;
}

template<class strategy_t>
void mesh_element3d<strategy_t>::init(simple_corners_t new_corners)
{
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

  i_type_t i;
  for (i = 0; i < n_corners; ++i)
    {
      corners[i].x = new_corners[i][0];
      corners[i].y = new_corners[i][1];
      corners[i].z = new_corners[i][2];
    }
}



template<class strategy_t>
inline void
mesh_element3d<strategy_t>::get_plane (element_plane_orientation_t orientation, plane_t &plane) const
{
  switch (orientation)
    {
      case x_axis_minus:  //left_cross
        plane[0] = corners[0];
        plane[1] = corners[2];
        plane[2] = corners[4];
        plane[3] = corners[6];
        break;
      case x_axis_plus:   //right_cross
        plane[0] = corners[1];
        plane[1] = corners[3];
        plane[2] = corners[5];
        plane[3] = corners[7];
        break;
      case y_axis_minus:  //top_cross
        plane[0] = corners[0];
        plane[1] = corners[1];
        plane[2] = corners[4];
        plane[3] = corners[5];
        break;
      case y_axis_plus:   //bottom_cross
        plane[0] = corners[2];
        plane[1] = corners[3];
        plane[2] = corners[6];
        plane[3] = corners[7];
        break;
      case z_axis_minus:  //lower_cross
        plane[0] = corners[0];
        plane[1] = corners[2];
        plane[2] = corners[1];
        plane[3] = corners[3];
        break;
      case z_axis_plus:   //upper_cross
        plane[0] = corners[4];
        plane[1] = corners[6];
        plane[2] = corners[5];
        plane[3] = corners[7];
        break;
      
      
      default:
        bs_throw_exception ("Invalid orientation!");;
    }
}

template<class strategy_t>
typename mesh_element3d<strategy_t>::fpoint3d_t
mesh_element3d<strategy_t>::get_center () const
{
  fpoint3d_t center;
  i_type_t i;
  for (i = 0; i < n_corners; i++)
    {
      center += corners[i];
    }

  return center / corners.size();
}




template<class strategy_t>
typename mesh_element3d<strategy_t>::fp_type_t
mesh_element3d<strategy_t>::calc_volume()
{
    fpoint3d_t center = get_center ();
    fp_type_t volume = 0.0;

    //share for 12 tetraidr (6 side -> 2 tetraidr for each)
    volume += calc_tetra_volume(corners[0],corners[1],corners[2], center);
    volume += calc_tetra_volume(corners[1],corners[2],corners[3], center);

    volume += calc_tetra_volume(corners[1],corners[3],corners[7], center);
    volume += calc_tetra_volume(corners[1],corners[5],corners[7], center);

    volume += calc_tetra_volume(corners[0],corners[2],corners[6], center);
    volume += calc_tetra_volume(corners[0],corners[4],corners[6], center);

    volume += calc_tetra_volume(corners[4],corners[5],corners[6], center);
    volume += calc_tetra_volume(corners[5],corners[6],corners[7], center);

    volume += calc_tetra_volume(corners[0],corners[1],corners[4], center);
    volume += calc_tetra_volume(corners[1],corners[4],corners[5], center);

    volume += calc_tetra_volume(corners[2],corners[3],corners[6], center);
    volume += calc_tetra_volume(corners[3],corners[6],corners[7], center);

    return volume;
}

template<class strategy_t>
typename mesh_element3d<strategy_t>::fp_type_t
mesh_element3d<strategy_t>::get_dx ()
{
  plane_t plane1, plane2;
  
  get_plane(x_axis_minus, plane1);
  get_plane(x_axis_plus, plane2);

  fp_type_t dx = 0.0;

  for (i_type_t i = 0; i < n_plane_corners; ++i)
    dx += plane1[i].x - plane2[i].x;

  return fabs (dx / n_plane_corners);
}

template<class strategy_t>
typename mesh_element3d<strategy_t>::fp_type_t
mesh_element3d<strategy_t>::get_dy ()
{
  plane_t plane1, plane2;
  
  get_plane(y_axis_minus, plane1);
  get_plane(y_axis_plus, plane2);
  
  fp_type_t dy = 0.0;

  for (i_type_t i = 0; i < n_plane_corners; ++i)
    dy += plane1[i].y - plane2[i].y;

  return fabs (dy / n_plane_corners);
}

template<class strategy_t>
typename mesh_element3d<strategy_t>::fp_type_t
mesh_element3d<strategy_t>::get_dz ()
{
  plane_t plane1, plane2;
  
  get_plane(z_axis_minus, plane1);
  get_plane(z_axis_plus, plane2);

  fp_type_t dz = 0.0;

  for (i_type_t i = 0; i < n_plane_corners; ++i)
    dz += plane1[i].x - plane2[i].x;

  return fabs (dz / n_plane_corners);
}

template<class strategy_t>
typename mesh_element3d<strategy_t>::point3d_t
mesh_element3d<strategy_t>::get_dx_dy_dz ()
{
  plane_t plane1, plane2, plane3, plane4, plane5, plane6;
  
  get_plane(x_axis_minus, plane1);
  get_plane(x_axis_plus, plane2);
  get_plane(y_axis_minus, plane3);
  get_plane(y_axis_plus, plane4);
  get_plane(z_axis_minus, plane5);
  get_plane(z_axis_plus, plane6);
  
  point3d_t element_size;
  element_size[0] = 0;
  element_size[1] = 0;
  element_size[2] = 0;

  for (i_type_t i = 0; i < n_plane_corners; ++i)
    {
      element_size[0] += plane1[i].x - plane2[i].x;
      element_size[1] += plane3[i].y - plane4[i].y;
      element_size[2] += plane6[i].z - plane6[i].z;
    }

  return element_size;
}


BS_INST_STRAT(mesh_element3d);
