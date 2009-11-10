/*! 
  \file fpoint3d.h
  \brief addition structure for mesh_grdecl - point3d
  \author Iskhakov Ruslan
  \date 2008-05-01 
 * */

#ifndef FPOINT3D_H
#define FPOINT3D_H

#include "fpoint2d.h"

namespace grd_ecl
{
  struct fpoint3d
  {
    float x,y,z;

    fpoint3d()
    {
      x=0.0f,y=0.0f,z=0.0f;
    }
    fpoint3d(float ax, float ay, float az)
    {
      x = ax;
      y = ay;
      z= az;
    }
    fpoint3d& operator += (const fpoint3d &b)
    {
      x += b.x;
      y += b.y;
      z += b.z;
      return *this;
    }
    fpoint3d& operator *= (const float b)
    {
      x *= b;
      y *= b;
      z *= b;
      return *this;
    }
    fpoint3d& operator *= (const fpoint3d &b)
    {
      x *= b.x;
      y *= b.y;
      z *= b.z;
      return *this;
    }
    fpoint3d& operator /= (const float b)
    {
      x /= b;
      y /= b;
      z /= b;
      return *this;
    }
    bool operator < (const fpoint3d &a)
    {
      if (a.z > z)
        return true;
      return false;
    }

    inline float get_length();

    static fpoint3d abs(const fpoint3d &a)
    {
      return fpoint3d (fabs(a.x),fabs(a.y),fabs(a.z));
    }
  };

  inline fpoint3d operator + (const fpoint3d &a, const fpoint3d  &b)
  {
    return fpoint3d(a.x+b.x,a.y+b.y,a.z+b.z);
  }

  inline fpoint3d operator - (const fpoint3d &a, const fpoint3d &b)
  {
    return fpoint3d(a.x-b.x,a.y-b.y,a.z-b.z);
  }

  inline fpoint3d operator / (const fpoint3d &a, const float d)
  {
    return fpoint3d(a.x/d,a.y/d,a.z/d);
  }

  inline float fpoint3d_mixed_multiplicate(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c)
  {
    return a.x*b.y*c.z+a.y*b.z*c.x+a.z*b.x*c.y-a.z*b.y*c.x-a.y*b.x*c.z-a.x*b.z*c.y;
  }

  inline fpoint3d operator ^ (const fpoint3d &u, const fpoint3d &v) //vector multiplication
  {
    return fpoint3d(u.y*v.z - u.z*v.y,
                    u.z*v.x - u.x*v.z,
                    u.x*v.y - u.y*v.x);
  }

  inline fpoint3d calc_normal (const fpoint3d &a, const fpoint3d &b, const fpoint3d &c)
  {
    fpoint3d n = (b-a)^(c-a);
    float len = n.get_length();
    n /= len;
    return n;
  }

  inline float fpoint3d::get_length() //get length of radius-vector
  {
    return sqrt(x*x+y*y+z*z);
  }

  inline float operator *(const fpoint3d &a, const fpoint3d &b)
  {
    return a.x*b.x+a.y*b.y+a.z*b.z;
  }

  //get cosinus between plane(a,b,c) and d - normal of other plane
  inline float get_cosinus(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c, const fpoint3d &d)
  {
    return calc_normal(a,b,c)*d;
  }

  inline float get_len(const fpoint3d &a, const fpoint3d &b)
  {
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
  }

  inline float formula_gerona3d(const fpoint3d &p1, const fpoint3d &p2, const fpoint3d &p3)
  {
    float a = (p2-p1).get_length();
    float b = (p3-p1).get_length();
    float c = (p3-p2).get_length();

    float p = (float)(a+b+c)/2;
    return sqrtf(p*(p-a)*(p-b)*(p-c));
  }

  inline float find_area_of_side(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c, const fpoint3d &d)
  //get square of side (we know the side numeration)
  {
    float area = formula_gerona3d(a,b,c);
    area += formula_gerona3d(c,d,b);
    return area;
  }
  //calc volume of tetraidr
  inline float calc_tetra_volume(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c, const fpoint3d &d)
  {
    float volume = fpoint3d_mixed_multiplicate(d-a,d-b,d-c);
    return fabs(volume/6.0f);
  }

  inline fpoint3d get_projection_on_all_axis_for_one_side(const boost::array <fpoint3d, 4> &vec)
  {
    fpoint3d A;

    boost::array <fpoint2d, 4> vec2;
    //YoZ
    for (size_t i = 0; i < vec.size(); i++)
      vec2[i] = fpoint2d(vec[i].y,vec[i].z);

    A.x = formula_gerona2d(vec2[0],vec2[1],vec2[2]) + formula_gerona2d(vec2[2],vec2[3],vec2[1]);

    //XoZ
    for (size_t i = 0; i < vec.size(); i++)
      vec2[i] = fpoint2d(vec[i].x,vec[i].z);

    A.y = formula_gerona2d(vec2[0],vec2[1],vec2[2]) + formula_gerona2d(vec2[2],vec2[3],vec2[1]);

    //XoY
    for (size_t i = 0; i < vec.size(); i++)
      vec2[i] = fpoint2d(vec[i].x,vec[i].y);

    A.z = formula_gerona2d(vec2[0],vec2[1],vec2[2]) + formula_gerona2d(vec2[2],vec2[3],vec2[1]);

    return A;
  }
  /*! \brief calculate intersection point for plane and segment
      \param side - vector of sides points
      \param center_side - center of side (for making decision about intersection_point inside side or not)
      \param start_point - ends of segment
      \param delta - difference = end_point - start_point
      \param dif_epsilon - difference between almost equal numbers
      \param intersection_point - return value
      \return false if no intersection*/

  inline bool 
  calc_intersection_of_plane_and_segment(const fpoint3d &p1, const fpoint3d &p2, const fpoint3d &p3, 
                                         const fpoint3d &center_side,
                                         const fpoint3d &start_point, 
                                         const fpoint3d &delta, 
                                         fpoint3d &intersection_point, 
                                         double eps_diff=1.0e-7)
  {
    //get plane a,b,c,d params
    double a,b,c,d;
    a = (p1.y-p2.y)*(p1.z+p2.z)+(p2.y-p3.y)*(p2.z+p3.z)+(p3.y-p1.y)*(p3.z+p1.z);
    b = (p1.z-p2.z)*(p1.x+p2.x)+(p2.z-p3.z)*(p2.x+p3.x)+(p3.z-p1.z)*(p3.x+p1.x);
    c = (p1.x-p2.x)*(p1.y+p2.y)+(p2.x-p3.x)*(p2.y+p3.y)+(p3.x-p1.x)*(p3.y+p1.y);
    d = -(a*p1.x + b*p1.y + c*p1.z);

    //find t params
    double n1 = a*start_point.x + b*start_point.y + c*start_point.z + d;
    double n2 = a*delta.x + b*delta.y + c*delta.z;
    if (!n2)
      return false;
    double t = -n1/n2;

    if (t > 1)
      t = 1;
    if (t < 0)
      t = 0;
    intersection_point.x = start_point.x + t*delta.x;
    intersection_point.y = start_point.y + t*delta.y;
    intersection_point.z = start_point.z + t*delta.z;

    //is it between side[0],side[1],side[2]
    double len = get_len(center_side,p1);
    double len_c = get_len(intersection_point,center_side);

    if (len_c < len || (fabs(len_c-len) < eps_diff))
      return true;
    else
      return false;
  }

};//namespace grd_ecl
#endif
