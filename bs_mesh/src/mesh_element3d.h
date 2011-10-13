#ifndef MESH_ELEMENT3D_H
#define MESH_ELEMENT3D_H
/*!
	\file mesh_element3d.h
	\brief This file declares mesh element in 3 dimensions
	\author Mark Khait
	\date 2009-09-23
 */

#include "fpoint3d.h"

using namespace blue_sky;

  //! element plane orientation
  enum element_plane_orientation_t
    {
      x_axis_plus = 0,
      x_axis_minus,
      y_axis_plus,
      y_axis_minus,
      z_axis_plus,
      z_axis_minus
    };
    
  //!
 element_plane_orientation_t get_reverse_orientation (element_plane_orientation_t orient);

#define N_ELEMENT_CORNERS 8
#define N_PLANE_CORNERS 4
    

  
  class BS_API_PLUGIN mesh_element3d
    {

    //-----------------------------------------
    // TYPES
    //-----------------------------------------
    
    public:
    
      ///////////////////////
      // OWN TYPES
      ///////////////////////

      typedef grd_ecl::fpoint3d                             fpoint3d_t;
      
      typedef boost::array <fpoint3d_t, N_ELEMENT_CORNERS>  corners_t;
      typedef boost::array <fpoint3d_t, N_PLANE_CORNERS>    plane_t;
      typedef boost::array <t_double, 3>                     point3d_t;
      typedef boost::array <point3d_t, N_ELEMENT_CORNERS>   simple_corners_t;

    //-----------------------------------------
    //  METHODS
    //-----------------------------------------

    public:
      
      ///////////////////////
      // INIT
      ///////////////////////

      //! default constructor
      mesh_element3d ();

      //! default destructor
      ~mesh_element3d ()	{};

      //! init element with simple corners
      void init (simple_corners_t new_corners);
      
      //! init element with corners
      void init (corners_t new_corners) {corners = new_corners;};
      
      ///////////////////////
      // ACCESS
      ///////////////////////
      
      //! return element corners
      corners_t 
      get_corners ()const
      {
        return corners;
      }
      
      //! return element simple corners
      simple_corners_t 
      get_simple_corners ()const
      {
        simple_corners_t simple_corners;
        return simple_corners;
      }

      //! return element simple corners
      void get_plane (element_plane_orientation_t orientation, plane_t &plane) const;
      
      
      ///////////////////////
      // WORK
      ///////////////////////

      //! get center of element
      fpoint3d_t get_center () const;

      //! find volume of the element
      t_double calc_volume();

      t_double get_dx ();
      t_double get_dy ();
      t_double get_dz ();

      point3d_t get_dx_dy_dz ();

    //-----------------------------------------
    //  VARIABLES
    //-----------------------------------------

    protected:
      corners_t corners;

      t_int n_corners, n_plane_corners;
    };

//! get center of plane
inline void get_plane_center (const boost::array <grd_ecl::fpoint3d, N_PLANE_CORNERS> &plane, grd_ecl::fpoint3d& center)
  {
    t_int i;
    for (i = 0; i < N_PLANE_CORNERS; i++)
      {
        center += plane[i];
      }
    center /= N_PLANE_CORNERS;
  }

#endif //MESH_ELEMENT3D_H
