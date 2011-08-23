#include "bs_mesh_stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "conf.h"
#include "export_python_wrapper.h"

namespace bp = boost::python;

namespace {
using namespace std;


struct point3d
{
    t_float x,y,z;

    // FIXME write set_xyz function

};

struct point2d
{
    t_float x,y;
};

class index2d
{
    public:
        t_int x, y;

        index2d()
        {
            x = 0;
            y = 0;
        }
        index2d(t_int index, t_int nx)
        {
            x = index % nx;
            y = index / nx;
        }
};

// finds (x,y,z) coordinates of line & plane intersection point
// plane is specified by points N1, N2 and vector (0,0,1)
// line is specified by points M1, M2
point3d find_edge_point(point3d N1, point3d N2, point3d M1, point3d M2)
{
   point3d X;
   t_float t, A, B, C, D, E, F, G;

   A = N2.y - N1.y;
   B = N2.x - N1.x;

   C = M2.x - M1.x;
   D = M2.y - M1.y;
   E = M2.z - M1.z; 

   F = N2.x - M1.x;
   G = N2.y - M1.y;

   t = (A*F - B*G)/(A*C - B*D);
   
   X.x = M1.x + C*t; 
   X.y = M1.y + D*t;
   X.z = M1.z + E*t;
   
   return X;
}


// finds z-coordinate of line & plane intersection point
// plane is specified by points N1, N2 and vector (0,0,1)
// line is specified by points M1, M2
t_float find_edge_point_z(point2d N1, point2d N2, point3d M1, point3d M2)
{
   t_float t, A, B, C, D, E, F, G, z;

   A = N2.y - N1.y;
   B = N2.x - N1.x;

   C = M2.x - M1.x;
   D = M2.y - M1.y;
   E = M2.z - M1.z; 

   F = N2.x - M1.x;
   G = N2.y - M1.y;

   t = (A*F - B*G)/(A*C - B*D);
   
   z = M1.z + E*t;
   
   return z;
}

// finds z-coordinate of line & plane intersection point
// plane is specified by 3 points N1, N2, N3
// line is specified by point M and vector (0,0,1)
t_float find_point_z(point3d N1, point3d N2, point3d N3, point2d M)
{
   t_float A, B, C, z;

   A = (N1.y - N2.y)*(N1.z - N3.z) - (N1.y - N3.y)*(N1.z - N2.z);
   B = (N1.x - N3.x)*(N1.z - N2.z) - (N1.x - N2.x)*(N1.z - N3.z);
   C = (N1.x - N2.x)*(N1.y - N3.y) - (N1.x - N3.x)*(N1.y - N2.y);
   
   z = N1.z + (A*(N1.x - M.x) + B*(N1.y - M.y))/C;
   
   return z;
}



// finds well-mesh intersection
// returns mesh_points, well_points
bp::tuple make_projection(t_int nx, t_int ny, t_int nz, 
                          spv_int indices_, spv_int faces_, 
                          spv_int internal_,
                          spv_float x_, spv_float y_,   
                          spv_float tops_, spv_float values_, 
                          t_int z_start, t_int z_end) 
{
    v_float& tops = *tops_;
    v_float& values = *values_;
    v_float& cross_x = *x_;
    v_float& cross_y = *y_;
    v_int& indices = *indices_;
    v_int& faces = *faces_;
    v_int& internal = *internal_;


    t_int i, ii, j, n, z; 
    n = indices.size();

    // FIXME
    if (n<2)
        return bp::make_tuple(0,0,0,0,0);
    
    //////////////////////////////////////////////

    t_int n_points, n_z_points, n_l_points;
    n_z_points = z_end-z_start+2;

    list <t_float> points, wpoints, vals;

    spv_float proj_mesh, well_points, scalars;

    proj_mesh = BS_KERNEL.create_object(v_float::bs_type());
    well_points = BS_KERNEL.create_object(v_float::bs_type());
    scalars = BS_KERNEL.create_object(v_float::bs_type());

    t_int my_ind, ind_old;
    index2d my;
    t_int face_in, face_out;

    point3d P[8];
    point2d M[2];

    t_float z1, z2;

    int faces2corners[5][4] = { {0,1,4,5},
                                {1,3,5,7},
                                {2,3,6,7},
                                {0,2,4,6},
                                {0,0,0,0}
                              }; 

    int corners[4];

    t_int index;

    t_float L = 0, l = 0;

    //! for the first column of intersected cells
    i = 0;
    ii = i+1;
    // FIXME change i -> 0

    my_ind = indices[i];
    my = index2d(my_ind, nx);
    ind_old = my_ind;

    //well-cell intersection points
    M[0].x = cross_x[i]; 
    M[0].y = cross_y[i];
    M[1].x = cross_x[ii]; 
    M[1].y = cross_y[ii];

    // well-cell intersection faces and face corners
    face_in = faces[i];

    if (face_in == 4)
    {
        //for the first cell in column
        //find intersection of the cutting plane and left upper edge
        // index = (my_ind.x; my_ind.y; z_start) ZYX
        index = z_start + my.y*nz + my.x*nz*ny;

        for (j=0;j<3;j++)
        {
          P[j].x = tops[8*3*index + 3*j];
          P[j].y = tops[8*3*index + 3*j + 1];
          P[j].z = tops[8*3*index + 3*j + 2];
        }
       
        z1 = find_point_z(P[0],P[1], P[2], M[0]);
           
        points.push_back(z1);
        points.push_back(0);

        //for all cells in column
        //find z-coordinates of intersection points of the cutting plane and left lower edges
        //index = (my_ind.x; my_ind.y; z)

        for (z=z_start; z<z_end+1; z++)
        {
           index = z + my.y*nz + my.x*nz*ny;
           vals.push_back(values[index]);

           for (j=4;j<7;j++)
           {
              P[j].x = tops[8*3*index + 3*j];
              P[j].y = tops[8*3*index + 3*j + 1];
              P[j].z = tops[8*3*index + 3*j + 2];
           }
           
           z2 = find_point_z(P[4],P[5], P[6], M[0]);

           points.push_back(z2);
           points.push_back(0);
        }
        
    }

    //! for all columns of intersected cells
    for (i=0;i<n-1;i++)
    {
        my_ind = indices[i];
        my = index2d(my_ind, nx);

        // well-cell intersection faces and face corners
        face_in = faces[i];

        if (internal[i])
            wpoints.push_back(L);
        
        ii = i + 1;

        // well-cell intersection points
        M[0].x = cross_x[i];
        M[0].y = cross_y[i];
        M[1].x = cross_x[ii]; 
        M[1].y = cross_y[ii];
        
        if (face_in == 4)
        {
            l = sqrt( (M[1].x-M[0].x)*(M[1].x-M[0].x) +  (M[1].y-M[0].y)*(M[1].y-M[0].y));
            L += l;
            continue;
        }

        // every index appears in indices twice: for face_in and face_out
        if (my_ind == ind_old)
           continue;
        ind_old = my_ind;

        corners[0] = faces2corners[face_in][0];
        corners[1] = faces2corners[face_in][1];
        corners[2] = faces2corners[face_in][2];
        corners[3] = faces2corners[face_in][3];


        //for the first cell in column
        //find intersection point of the cutting plane and upper right edge
        //index = (my_ind.x; my_ind.y; z_start) ZYX

        index = z_start + my.y*nz + my.x*nz*ny;

        j = corners[0];
        P[j].x = tops[8*3*index + 3*j];
        P[j].y = tops[8*3*index + 3*j + 1];
        P[j].z = tops[8*3*index + 3*j + 2];

        j = corners[1];
        P[j].x = tops[8*3*index + 3*j];
        P[j].y = tops[8*3*index + 3*j + 1];
        P[j].z = tops[8*3*index + 3*j + 2];

        z1 = find_edge_point_z(M[0], M[1], P[corners[0]], P[corners[1]]);

        points.push_back(z1);
        points.push_back(L);

        //for the other cells in column
        //find z-coordinates of intersection points of the cutting plane and lower right edge
        //index = (my_ind.x; my_ind.y; z) ZYX

        for (z=z_start; z<z_end+1; z++)
        {
           index = z + my.y*nz + my.x*nz*ny;
           
           vals.push_back(values[index]);

           j = corners[2];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];
           
           j = corners[3];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];
           
           z2 = find_edge_point_z(M[0], M[1], P[corners[2]], P[corners[3]]);

           points.push_back(z2);
           points.push_back(L);
        }
        
        l = sqrt( (M[1].x-M[0].x)*(M[1].x-M[0].x) +  (M[1].y-M[0].y)*(M[1].y-M[0].y));
        L += l;
        
    }


    //! for the last column of intersected cells
    i = n-1;
    ii = n-2;

    wpoints.push_back(L);

    my_ind = indices[i];
    my = index2d(my_ind, nx);

    //well-cell intersection points
    M[0].x = cross_x[i]; 
    M[0].y = cross_y[i];
    M[1].x = cross_x[ii]; 
    M[1].y = cross_y[ii];

    // well-cell intersection faces and face corners
    face_out = faces[i];
    corners[0] = faces2corners[face_out][0];
    corners[1] = faces2corners[face_out][1];
    corners[2] = faces2corners[face_out][2];
    corners[3] = faces2corners[face_out][3];

    //for the first cell in column
    //find intersection of the cutting plane and left upper edge
    // index = (my_ind.x; my_ind.y; z_start) ZYX
    index = z_start + my.y*nz + my.x*nz*ny;

    // for face_out
    if (face_out != 4)
    {
       j = corners[0];
       P[j].x = tops[8*3*index + 3*j];
       P[j].y = tops[8*3*index + 3*j + 1];
       P[j].z = tops[8*3*index + 3*j + 2];
       
       j = corners[1];
       P[j].x = tops[8*3*index + 3*j];
       P[j].y = tops[8*3*index + 3*j + 1];
       P[j].z = tops[8*3*index + 3*j + 2];
       
       z1 = find_edge_point_z(M[0], M[1], P[corners[0]], P[corners[1]]);

    }
    else
    {
       for (j=0;j<3;j++)
       {
          P[j].x = tops[8*3*index + 3*j];
          P[j].y = tops[8*3*index + 3*j + 1];
          P[j].z = tops[8*3*index + 3*j + 2];
       }
       
       z1 = find_point_z(P[0],P[1], P[2], M[0]);
    }    
       
    points.push_back(z1);
    points.push_back(L);

    //for all cells in column
    //find z-coordinates of intersection points of the cutting plane and left lower edges
    //index = (my_ind.x; my_ind.y; z)

    for (z=z_start; z<z_end+1; z++)
    {
       index = z + my.y*nz + my.x*nz*ny;

       if (face_out != 4)
       {
           j = corners[2];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];
           
           j = corners[3];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];

           z2 = find_edge_point_z(M[0], M[1], P[corners[2]], P[corners[3]]);
       }
       else
       {
           for (j=4;j<7;j++)
           {
              P[j].x = tops[8*3*index + 3*j];
              P[j].y = tops[8*3*index + 3*j + 1];
              P[j].z = tops[8*3*index + 3*j + 2];
           }
           z2 = find_point_z(P[4],P[5], P[6], M[0]);
       }
       points.push_back(z2);
       points.push_back(L);
    }

    n_points = points.size();
    proj_mesh->resize (n_points);
    well_points->resize (wpoints.size());
    scalars->resize (vals.size());

    copy(points.begin(), points.end(), &(*proj_mesh)[0]);
    copy(wpoints.begin(), wpoints.end(), &(*well_points)[0]);
    copy(vals.begin(), vals.end(), &(*scalars)[0]);

    n_l_points = points.size()/n_z_points/2;

    return bp::make_tuple(proj_mesh, well_points, scalars, n_z_points, n_l_points); 

}

}
namespace blue_sky { namespace python {

 
void py_export_well_edit() {
    def("make_projection", &make_projection);
}

}} 	// eof blue_sky::python

