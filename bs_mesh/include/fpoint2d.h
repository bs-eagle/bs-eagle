#ifndef FPOINT2D_H
#define FPOINT2D_H

/*! \file fpoint2d.h
	\brief addition structure for mesh_grdecl - point2d
	\author Iskhakov Ruslan
	\date 2008-05-01 */

namespace grd_ecl
  {
  struct fpoint2d
    {
      float x,y;

      fpoint2d()
      {
        x=0.0f,y=0.0f;
      }
      fpoint2d(float ax, float ay)
      {
        x = ax;
        y = ay;
      }
      fpoint2d& operator += (const fpoint2d &b)
      {
        x += b.x;
        y += b.y;
        return *this;
      }
      fpoint2d& operator /= (float b)
      {
        x /= b;
        y /= b;
        return *this;
      }
    };
  inline float get_sq_distance(const fpoint2d &a, const fpoint2d &b)
  {
    return (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y);
  }
  inline float get_len(const fpoint2d &a, const fpoint2d &b)
  {
    return sqrtf(get_sq_distance(a,b));
  }

  inline fpoint2d operator + (const fpoint2d &a, const fpoint2d  &b)
  {
    return fpoint2d(a.x+b.x,a.y+b.y);
  }
  inline fpoint2d operator - (fpoint2d &a, fpoint2d &b)
  {
    return fpoint2d(a.x-b.x,a.y-b.y);
  }
  inline fpoint2d operator / (const fpoint2d &a, float d)
  {
    return fpoint2d(a.x/d,a.y/d);
  }

  //additional functions
  inline int crossing_points_2d( //crossing of 2 cut
    const fpoint2d &p11, const fpoint2d &p12,   // coordinates of first cut
    const fpoint2d &p21, const fpoint2d &p22,   // coordinates of second cut
    fpoint2d &resPoint, double eps)
  {

    double Z  = (p12.y-p11.y)*(p21.x-p22.x)-(p21.y-p22.y)*(p12.x-p11.x); //denominator
    //if denominator = 0 => lines are parallel
    if ( fabs(Z) < eps)
      return 0;

    double Ua = ( (p12.y-p11.y)*(p21.x-p11.x)-(p21.y-p11.y)*(p12.x-p11.x) ) / Z; //numerator1
    double Ub = ( (p21.y-p11.y)*(p21.x-p22.x)-(p21.y-p22.y)*(p21.x-p11.x) ) / Z; //numerator2

    resPoint.x = (float)(p11.x + (p12.x - p11.x) * Ub);
    resPoint.y = (float)(p11.y + (p12.y - p11.y) * Ub);

    // if 0<=Ua<=1 && 0<=Ub<=1 -> crossing point inside cut
    if ( (0 <= Ua)&&(Ua <= 1)&&(0 <= Ub)&&(Ub <= 1) )
      return 1;
    return 0;
  }
  inline float formula_gerona2d(const fpoint2d &p1, const fpoint2d &p2, const fpoint2d &p3)
  {
    float a = get_len(p2,p1);
    float b = get_len(p1,p3);
    float c = get_len(p3,p2);

    float p = (a+b+c)/2;
    float x = p*(p-a)*(p-b)*(p-c);
    if (x < 1.0e-16)
      {
        return 0;
      }
    return sqrtf(x);
  }

  inline float get_triangle_crossing_area(const std::vector<fpoint2d> &tri1, const  std::vector<fpoint2d> &tri2, double eps)
  {
    std::vector<fpoint2d> f;
    fpoint2d p;
    //side 0-1
    if (crossing_points_2d(tri1[0],tri1[1],tri2[0],tri2[1],p,eps))
      f.push_back(p);
    if (crossing_points_2d(tri1[0],tri1[1],tri2[1],tri2[2],p,eps))
      f.push_back(p);
    if (crossing_points_2d(tri1[0],tri1[1],tri2[0],tri2[2],p,eps))
      f.push_back(p);
    //side 1-2
    if (crossing_points_2d(tri1[2],tri1[1],tri2[0],tri2[1],p,eps))
      f.push_back(p);
    if (crossing_points_2d(tri1[2],tri1[1],tri2[1],tri2[2],p,eps))
      f.push_back(p);
    if (crossing_points_2d(tri1[2],tri1[1],tri2[0],tri2[2],p,eps))
      f.push_back(p);

    //side 0-2
    if (crossing_points_2d(tri1[0],tri1[2],tri2[0],tri2[1],p,eps))
      f.push_back(p);

    if (!f.size()) //2 point are not triangle
      return 0.0;

    if (crossing_points_2d(tri1[0],tri1[2],tri2[1],tri2[2],p,eps))
      f.push_back(p);
    if (crossing_points_2d(tri1[0],tri1[2],tri2[0],tri2[2],p,eps))
      f.push_back(p);

    //we can have 1 extra point
    size_t ik = 0;
    for (size_t i = 0; i < f.size()-1 && !ik; i++)
      for (size_t j = i+1; j < f.size(); j++)
        if (get_sq_distance(f[i],f[j]) < eps)
          {
            ik = j;
            break;
          }

    if (ik) //got out extra point
      {
        std::vector<fpoint2d>::iterator i;
        i = f.begin()+ik;
        f.erase(i);
      }
    if (f.size() < 3)
      return 0.0;

    float area = formula_gerona2d(f[0],f[1],f[2]);

    if (f.size() == 4) //add second triangle
      area += formula_gerona2d(f[1],f[2],f[3]);

    return area;
  }

  inline float get_polygon_crossing_area(const std::vector<fpoint2d> &polyg1, const std::vector<fpoint2d> &polyg2, double eps)
  {
    //we have only 4-polygon -> we can share it on 2 triangles and after
    //find triangle crossing area -> for our triangles they can
    //have 0, 3 or 4 crossing point->using geron formula we can get the area -> sum = sum(area-triangle)
    //всего может получитс€ 4 больших пересечени€ - искома€ площадь - их сумма (подходит дл€ всех случаев)

    //share first polygon
    std::vector<fpoint2d> tri1, tri2;
    tri1.push_back(polyg1[0]);
    tri1.push_back(polyg1[1]);
    tri1.push_back(polyg1[2]);

    tri2.push_back(polyg1[2]);
    tri2.push_back(polyg1[3]);
    tri2.push_back(polyg1[1]);

    //share second polygon
    std::vector<fpoint2d> tri3, tri4;
    tri3.push_back(polyg2[0]);
    tri3.push_back(polyg2[1]);
    tri3.push_back(polyg2[2]);

    tri4.push_back(polyg2[2]);
    tri4.push_back(polyg2[3]);
    tri4.push_back(polyg2[1]);

    //find their area
    float area = 0.0;
    area += get_triangle_crossing_area(tri1, tri3,eps);
    area += get_triangle_crossing_area(tri1, tri4,eps);
    area += get_triangle_crossing_area(tri2, tri3,eps);
    area += get_triangle_crossing_area(tri2, tri4,eps);

    return area;
  }
};//namespace grd_ecl
#endif
