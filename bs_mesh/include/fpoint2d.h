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
  
  typedef boost::array <fpoint2d, 4>       quadrangle_t;
  typedef boost::array <fpoint2d, 3>       triangle_t;
  
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
  // find cross point of two cuts, return 1 if point is INSIDE both of cuts
  inline int crossing_points_2d( //crossing of 2 cut
    const fpoint2d &p1, const fpoint2d &p2,   // coordinates of first cut
    const fpoint2d &p3, const fpoint2d &p4,   // coordinates of second cut
    fpoint2d &resPoint, double eps)
  {

    double Z  = (p4.y - p3.y) * (p2.x - p1.x) - (p2.y - p1.y) * (p4.x - p3.x); //denominator
    //if denominator = 0 => lines are parallel
    if ( fabs(Z) < eps)
      return 0;

    double Ua = ((p4.x - p3.x) * (p1.y - p3.y) - (p4.y - p3.y) * (p1.x - p3.x) ) / Z; //numerator1
    double Ub = ((p2.x - p1.x) * (p1.y - p3.y) - (p2.y - p1.y) * (p1.x - p3.x) ) / Z; //numerator2

    resPoint.x = (float)(p1.x + (p2.x - p1.x) * Ua);
    resPoint.y = (float)(p1.y + (p2.y - p1.y) * Ua);

    if ( (0 < Ua) && (Ua < 1)&&(0 < Ub) && (Ub < 1) )
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

  inline float get_triangle_crossing_area(const triangle_t &tri1, const  triangle_t &tri2, double eps, int are_opposed)
  {
    boost::array <fpoint2d, 6> f;
    fpoint2d p, p1;
    int i, n = 0;
    int vertice_match1 = 0, vertice_match2 = 0;
    float area = 0.0;
    
    
    if (!are_opposed && (tri1[2].x == tri2[2].x) && (tri1[2].y == tri2[2].y))
      {
        // edges do not cross each other
        f[n++] = tri1[2];
      }
    else
      {
        // search crossing points of triangles non-coord edges (edge between 0 and 1 point is always on COORD)
        
        if (are_opposed && (((tri1[2].x == tri2[0].x) && (tri1[2].y == tri2[0].y)) || ((tri1[2].x == tri2[1].x) && (tri1[2].y == tri2[1].y))))
          {
            // edge 0-2 of tri1 do not cross other edges
            f[n++] = tri1[2];
            vertice_match1 = 1;
          }
        else
          {
            if (crossing_points_2d(tri1[0], tri1[2], tri2[1], tri2[2], p, eps))
              f[n++] = p;
            if (crossing_points_2d(tri1[0], tri1[2], tri2[0], tri2[2], p, eps))
              f[n++] = p;
          }
          
        if (are_opposed && (((tri2[2].x == tri1[0].x) && (tri2[2].y == tri1[0].y)) || ((tri2[2].x == tri1[1].x) && (tri2[2].y == tri1[1].y))))
          {
            // edge 0-2 of tri2 do not cross other edges
            f[n++] = tri2[2];
            vertice_match2 = 1;
          }
        else
          {
            if (crossing_points_2d(tri1[2], tri1[1], tri2[1], tri2[2], p, eps))
              f[n++] = p;
            if (crossing_points_2d(tri1[2], tri1[1], tri2[0], tri2[2], p, eps))
              f[n++] = p;
          }
      }
      
    // here we check points on same COORD, using y coordinate, which is z-coordinate in 3D
    if (are_opposed)
      {
        if (!vertice_match1 && (tri1[2].y > tri2[0].y) && (tri1[2].y < tri2[1].y))
           f[n++] = tri1[2];
           
        if (!vertice_match2 && (tri2[2].y > tri1[0].y) && (tri2[2].y < tri1[1].y))
           f[n++] = tri2[2];  
      }
    else
      {
        if (tri1[0].y > tri2[0].y)
          p = tri1[0];
        else
          p = tri2[0];
          
        if (tri1[1].y < tri2[1].y)
          p1 = tri1[1];
        else
          p1 = tri2[1];
        
        if (p.y < p1.y)
          {
            f[n++] = p;
            f[n++] = p1;  
          }
      }
    
    n -= 2;
    
    for (i = 0; i < n; ++i)
      area += formula_gerona2d(f[0],f[i + 1],f[i + 2]);  
    
    return area;
  }

  inline float get_quadrangle_crossing_area(const quadrangle_t &polyg1, const quadrangle_t &polyg2, double eps)
  {
    // devide each quadrangle into 2 triangles and find 4 areas of triangles crossing
    // crossing area calculation is greatly simplified using MESH_GRDECL properties

    //share first polygon
    triangle_t tri1, tri2;
    tri1[0] = polyg1[0];  // COORD1
    tri1[1] = polyg1[2];  // COORD1
    tri1[2] = polyg1[1];  // COORD2 
    
    tri2[0] = polyg1[1];  // COORD2
    tri2[1] = polyg1[3];  // COORD2
    tri2[2] = polyg1[2];  // COORD1


    //share second polygon
    triangle_t tri3, tri4;
    tri3[0] = polyg2[0];  // COORD1
    tri3[1] = polyg2[2];  // COORD1
    tri3[2] = polyg2[1];  // COORD2 
    
    tri4[0] = polyg2[1];  // COORD2
    tri4[1] = polyg2[3];  // COORD2
    tri4[2] = polyg2[2];  // COORD1

    //find their area
    float area = 0.0;
    
    // tri1 and tri4, tri2 and tri3 are opposed (the have edge on different COORD)
    area += get_triangle_crossing_area(tri1, tri4, eps, true);
    area += get_triangle_crossing_area(tri2, tri3, eps, true);
    
    // tri1 and tri3, tri2 and tri4 are not opposed (the have edge on same COORD)
    area += get_triangle_crossing_area(tri1, tri3, eps, false);
    area += get_triangle_crossing_area(tri2, tri4, eps, false);

    return area;
  }
};//namespace grd_ecl
#endif
