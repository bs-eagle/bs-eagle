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

};

class index3d
{
    public:
    t_int x, y, z;

    index3d()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    index3d(t_int index, t_int ny, t_int nz)
    {
        z = index % nz;
        y = (index / nz) % ny;
        x = (index / nz) / ny;
    }

};

// mesh2d points from coord
//spv_float mesh2d_from_coord(spv_float coord_, t_int Nx, t_int Ny)
//{
//    spv_float mesh2d;
//    v_float& coord = *coord_;
//    t_int i, j, n;
//
//    n = (Nx+1)*(Ny+1);
//
//    mesh2d = BS_KERNEL.create_object(v_float::bs_type());
//    mesh2d->resize(3*n);
//
//    t_float *mesh;
//    mesh = &(*mesh2d)[0];
//
//    for (i=0;i<n;i++)
//    {
//        for (j=0;j<2;j++)
//            mesh[3*i+j] = coord[6*i+j];
//        mesh[3*i+2] = 0;
//    }
//
//    return mesh2d;
//}

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
t_float find_edge_point_z(point3d N1, point3d N2, point3d M1, point3d M2)
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
t_float find_point_z(point3d N1, point3d N2, point3d N3, point3d M)
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
bp::tuple make_projection(t_int ny, t_int nz, 
                          spv_int indices_, spv_int faces_, 
                          spv_float x_, spv_float y_, spv_float z_,  spv_float tops_ ) 
{
    v_float& tops = *tops_;
    v_float& cross_x = *x_;
    v_float& cross_y = *y_;
    v_float& cross_z = *z_;
    v_int& indices = *indices_;
    v_int& faces = *faces_;


    t_int i, ii, j, n; 
    n = indices.size();


    t_int z, z_start = 1000000, z_end = 0; 

    // calculate z_start, z_end from indices[]
    for (i=0;i<n;i++)
    {
        z = index3d(indices[i], ny, nz).z;
        if (z < z_start)
            z_start = z;
        if (z > z_end)
            z_end = z;
    }

    printf("\n z_start = %d z_end = %d", z_start, z_end);

    //////////////////////////////////////////////

    t_int n_points, n_z_points, n_l_points;
    n_z_points = z_end-z_start+2;

    list <t_float> points, wpoints;

    spv_float proj_mesh, well_points;

    proj_mesh = BS_KERNEL.create_object(v_float::bs_type());
    well_points = BS_KERNEL.create_object(v_float::bs_type());

    t_int my_ind, ind_old;
    index3d my;
    t_int face_in, face_out;

    point3d P[8], M[2], WP;
    t_float z1, z2;

    int faces2corners[7][4] = { {0,1,2,3},
                                {0,1,4,5},
                                {4,5,6,7},
                                {2,3,6,7},
                                {0,2,4,6},
                                {1,3,5,7},
                                {0,0,0,0}
                              }; 

    int corners[8];

    t_int index;

    t_float L = 0, l = 0, wL = 0, wl = 0;

    // for the first column of intersected cells

    i = 0;
    ii = i+1;

    my_ind = indices[i];
    ind_old = my_ind;

    my = index3d(my_ind, ny ,nz);

    printf("\n i=%d my_ind=%d x=%d y=%d z=%d", i, my_ind, my.x, my.y, my.z);

    //FIXME 
    while (faces[ii] == 6 && ii<n-1)
       ii ++;

    //well-cell intersection points
    M[0].x = cross_x[i]; 
    M[0].y = cross_y[i];
    M[0].z = cross_z[i];
    M[1].x = cross_x[ii]; 
    M[1].y = cross_y[ii];
    M[1].z = cross_z[ii];

    // projection length
    l = sqrt( (M[1].x-M[0].x)*(M[1].x-M[0].x) +  (M[1].y-M[0].y)*(M[1].y-M[0].y));
    printf("\n l = %f", l);

    //FIXME always add first well point
    //{
       printf("\n well point");
       //WP.x = cross_x[i]; 
       //WP.y = cross_y[i];
       WP.z = cross_z[i];
       //wl = sqrt( (WP.x-M[0].x)*(WP.x-M[0].x) +  (WP.y-M[0].y)*(WP.y-M[0].y));
       //wL = L - l + wl;
       //wpoints.push_back(WP.z);
       wpoints.push_back(wL);
       printf("\n L=%f l=%f wl=%f wL=%f \n", L, l, wl, wL);
    //}

    face_in = faces[i];
    face_out = faces[ii];

    corners[0] = faces2corners[face_in][0];
    corners[1] = faces2corners[face_in][1];
    corners[2] = faces2corners[face_in][2];
    corners[3] = faces2corners[face_in][3];
    corners[4] = faces2corners[face_out][0];
    corners[5] = faces2corners[face_out][1];
    corners[6] = faces2corners[face_out][2];
    corners[7] = faces2corners[face_out][3];
    //printf("\n %f %f %f %f %f %f", M[0].x, M[0].y, M[0].z, M[1].x, M[1].y, M[1].z);

    //for the first cell in column
    //find intersection of the cutting plane and left upper edge
    // index = (my_ind.x; my_ind.y; z_start) ZYX

    index = z_start + my.y*nz + my.x*nz*ny;
    printf("\n index=%d", index);

    // for face_in
    if (face_in !=0 && face_in != 2 && face_in != 6)
    {
       
       //FIXME dont calculate all 8 points, calculate as much as needed
       //for (j=0;j<8;j++)
       //{
          //P[j].x = tops[8*3*index + 3*j];
          //P[j].y = tops[8*3*index + 3*j + 1];
          //P[j].z = tops[8*3*index + 3*j + 2];
       //}
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

    printf("\n z1 = %f L = %f", z1, L);


    //for all cells in column
    //find z-coordinates of intersection points of the cutting plane and left lower edges
    //index = (my_ind.x; my_ind.y; z)

    for (z=z_start; z<z_end+1; z++)
    {
       index = z + my.y*nz + my.x*nz*ny;

       if (face_in !=0 && face_in != 2 && face_in != 6)
       {
           //FIXME dont calculate all 8 points, calculate as much as needed
           //for (j=0;j<8;j++)
           //{
           //   P[j].x = tops[8*3*index + 3*j];
           //   P[j].y = tops[8*3*index + 3*j + 1];
           //   P[j].z = tops[8*3*index + 3*j + 2];
           //}
           
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
            //FIXME write function for z2 calc
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
       printf("\n z2 = %f L = %f", z2, L);
    }

    //for the first cell in column
    //find intersection of the cutting plane and right upper edge
    //index = (my_ind.x; my_ind.y; z_start) ZYX

    L += l;

    index = z_start + my.y*nz + my.x*nz*ny;
    printf("\n index=%d", index);

    // for face_out
    if (face_out !=0 && face_out != 2)
    {
       
       //FIXME dont calculate all 8 points, calculate as much as needed
       //for (j=0;j<8;j++)
       //{
       //   P[j].x = tops[8*3*index + 3*j];
       //   P[j].y = tops[8*3*index + 3*j + 1];
       //   P[j].z = tops[8*3*index + 3*j + 2];
       //}
       //
       j = corners[4];
       P[j].x = tops[8*3*index + 3*j];
       P[j].y = tops[8*3*index + 3*j + 1];
       P[j].z = tops[8*3*index + 3*j + 2];
       
       j = corners[5];
       P[j].x = tops[8*3*index + 3*j];
       P[j].y = tops[8*3*index + 3*j + 1];
       P[j].z = tops[8*3*index + 3*j + 2];
       
       z1 = find_edge_point_z(M[0], M[1], P[corners[4]], P[corners[5]]);

    points.push_back(z1);
    points.push_back(L);

    printf("\n z1 = %f L = %f", z1, L);
    }

    //for all cells in column
    //find z-coordinates of intersection points of the cutting plane and right lower edges
    //index = (my_ind.x; my_ind.y; z)

    for (z=z_start; z<z_end+1; z++)
    {
       index = z + my.y*nz + my.x*nz*ny;
       if (face_out !=0 && face_out != 2)
       {
           
           //FIXME dont calculate all 8 points, calculate as much as needed
          // for (j=0;j<8;j++)
          // {
          //    P[j].x = tops[8*3*index + 3*j];
          //    P[j].y = tops[8*3*index + 3*j + 1];
          //    P[j].z = tops[8*3*index + 3*j + 2];
          // }
           
           j = corners[6];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];
           
           j = corners[7];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];
           
           z2 = find_edge_point_z(M[0], M[1], P[corners[6]], P[corners[7]]);

       points.push_back(z2);
       points.push_back(L);

       printf("\n z2 = %f L = %f", z2, L);
       }
    }

    // for the other columns of intersected cells

    for (i=1;i<n-1;i++)
    {
       my_ind = indices[i];
       
       my = index3d(my_ind, ny ,nz);
       printf("\n i=%d my_ind=%d x=%d y=%d z=%d", i, my_ind, my.x, my.y, my.z);
       
       //FIXME add well point handler
       if (faces[i] == 6)
       {
           printf("\n\n well point i=%d my_ind=%d",i,my_ind);
           WP.x = cross_x[i]; 
           WP.y = cross_y[i];
           WP.z = cross_z[i];
           wl = sqrt( (WP.x-M[0].x)*(WP.x-M[0].x) +  (WP.y-M[0].y)*(WP.y-M[0].y));
           wL = L - l + wl;
           printf("\n L=%f l=%f wl=%f wL=%f \n", L, l, wl, wL);
           //wpoints.push_back(WP.z);
           wpoints.push_back(wL);
           continue;
       }

       if (my_ind == ind_old)
           continue;

       ind_old = my_ind;

       
       // find on-cell-face well point index
       ii = i + 1;
       while (faces[ii] == 6 && ii<n-1)
           ii ++;
       printf("\n i=%d ii=%d",i,ii);
       //well-cell intersection points
       M[0].x = cross_x[i]; 
       M[0].y = cross_y[i];
       M[0].z = cross_z[i];
       M[1].x = cross_x[ii]; 
       M[1].y = cross_y[ii];
       M[1].z = cross_z[ii];

       printf("\n M[0] %f %f %f M[1] %f %f %f", M[0].x, M[0].y, M[0].z, M[1].x, M[1].y, M[1].z);

       // projection length
       l = sqrt( (M[1].x-M[0].x)*(M[1].x-M[0].x) +  (M[1].y-M[0].y)*(M[1].y-M[0].y));

       //FIXME
       if (l == 0)
           continue;

       printf("\n l = %f", l);
       L += l;

       //face_in = faces[i];
       face_out = faces[ii];

       //corners[0] = faces2corners[face_in][0];
       //corners[1] = faces2corners[face_in][1];
       //corners[2] = faces2corners[face_in][2];
       //corners[3] = faces2corners[face_in][3];
       corners[4] = faces2corners[face_out][0];
       corners[5] = faces2corners[face_out][1];
       corners[6] = faces2corners[face_out][2];
       corners[7] = faces2corners[face_out][3];
       
       //for the first cell in column
       //find intersection point of the cutting plane and upper right edge
       //index = (my_ind.x; my_ind.y; z_start) ZYX
       
       index = z_start + my.y*nz + my.x*nz*ny;

       //printf("\n index=%d", index);
       
       // for face_out
       if (face_out !=0 && face_out != 2 && face_out != 6)
       {
           
           //FIXME dont calculate all 8 points, calculate as much as needed
           //for (j=0;j<8;j++)
           //{
           //   P[j].x = tops[8*3*index + 3*j];
           //   P[j].y = tops[8*3*index + 3*j + 1];
           //   P[j].z = tops[8*3*index + 3*j + 2];
           //}
           
           j = corners[4];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];
           
           j = corners[5];
           P[j].x = tops[8*3*index + 3*j];
           P[j].y = tops[8*3*index + 3*j + 1];
           P[j].z = tops[8*3*index + 3*j + 2];
           
           z1 = find_edge_point_z(M[0], M[1], P[corners[4]], P[corners[5]]);

           points.push_back(z1);
           points.push_back(L);
           printf("\n z1 = %f L = %f", z1, L);
        }
       
       //for the other cells in column
       //find z-coordinates of intersection points of the cutting plane and lower right edge
       //index = (my_ind.x; my_ind.y; z) ZYX

       for (z=z_start; z<z_end+1; z++)
       {
           index = z + my.y*nz + my.x*nz*ny;
           if (face_out !=0 && face_out != 2 && face_out != 6)
           {
               
               //FIXME dont calculate all 8 points, calculate as much as needed
               //for (j=0;j<8;j++)
               //{
               //   P[j].x = tops[8*3*index + 3*j];
               //   P[j].y = tops[8*3*index + 3*j + 1];
               //   P[j].z = tops[8*3*index + 3*j + 2];
               //}
               
               j = corners[6];
               P[j].x = tops[8*3*index + 3*j];
               P[j].y = tops[8*3*index + 3*j + 1];
               P[j].z = tops[8*3*index + 3*j + 2];
               
               j = corners[7];
               P[j].x = tops[8*3*index + 3*j];
               P[j].y = tops[8*3*index + 3*j + 1];
               P[j].z = tops[8*3*index + 3*j + 2];
               
               z2 = find_edge_point_z(M[0], M[1], P[corners[6]], P[corners[7]]);

               points.push_back(z2);
               points.push_back(L);
               printf("\n z2 = %f L = %f", z2, L);
             }
       }

    }

    // for the  last column of intersected cells

    i = n-1;
       //wpoints.push_back(M[1].z);
    wpoints.push_back(L);

    if (faces[i] == 0 || faces[i] == 2 || faces[i] == 6)
    {
       my_ind = indices[i];
       
       my = index3d(my_ind, ny ,nz);
       printf("\n i=%d my_ind=%d x=%d y=%d z=%d", i, my_ind, my.x, my.y, my.z);
       

       
       
       //for the first cell in column
       //find intersection point of the cutting plane and upper right edge
       //index = (my_ind.x; my_ind.y; z_start) ZYX
       
       index = z_start + my.y*nz + my.x*nz*ny;

       //printf("\n index=%d", index);
       
       // for face_out
       
           for (j=0;j<3;j++)
           {
              P[j].x = tops[8*3*index + 3*j];
              P[j].y = tops[8*3*index + 3*j + 1];
              P[j].z = tops[8*3*index + 3*j + 2];
           }
           z1 = find_point_z(P[0],P[1], P[2], M[1]);
       
           points.push_back(z1);
           points.push_back(L);
           printf("\n z1 = %f L = %f", z1, L);
        
       
       //for the other cells in column
       //find z-coordinates of intersection points of the cutting plane and lower right edge
       //index = (my_ind.x; my_ind.y; z) ZYX

       for (z=z_start; z<z_end+1; z++)
       {
           index = z + my.y*nz + my.x*nz*ny;
           
               //FIXME function to calculate z
               for (j=4;j<7;j++)
               {
                  P[j].x = tops[8*3*index + 3*j];
                  P[j].y = tops[8*3*index + 3*j + 1];
                  P[j].z = tops[8*3*index + 3*j + 2];
               }
               z2 = find_point_z(P[4],P[5], P[6], M[1]);
           
               points.push_back(z2);
               points.push_back(L);
               printf("\n z2 = %f L = %f", z2, L);
             
       }

    }


    n_points = points.size();
    proj_mesh->resize (n_points);
    well_points->resize (wpoints.size());
    copy(points.begin(), points.end(), &(*proj_mesh)[0]);
    copy(wpoints.begin(), wpoints.end(), &(*well_points)[0]);

    n_l_points = points.size()/n_z_points/2;

    printf("\n nl=%d nz=%d", n_l_points, n_z_points);

    return bp::make_tuple(proj_mesh, well_points, n_z_points, n_l_points); 
}

}

namespace blue_sky { namespace python {

 
void py_export_well_edit() {
    def("make_projection", &make_projection);
    //def("mesh2d_from_coord", &mesh2d_from_coord);
}

}} 	// eof blue_sky::python

