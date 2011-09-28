#!/usr/bin/python
import bs
import numpy as np
import bs.bs_mesh as bm

# generate initial mesh
cell_x = 25; cell_y = 25; cell_z = 1;
dx = np.array([cell_x], dtype='d') # Lx
dy = np.array([cell_y], dtype='d') # Ly
dz = np.array([cell_z], dtype='d') # Lz
Nx = 50
Ny = 50
Nz = 100
[c, z] = bm.mesh_grdecl.gen_coord_zcorn(Nx, Ny, Nz, dx, dy, dz)
print 'Source mesh generated!'

# generate random points
pnum = 500;
px = np.random.rand(pnum) * cell_x * Nx;
py = np.random.rand(pnum) * cell_y * Ny;
pz = np.random.rand(pnum) * cell_z * Nz;

# build array of points ccords
P = np.c_[px, py, pz];
print P

# find where are points
cell_id = bm.where_is_point(Nx, Ny, c, z, P.reshape([1, -1]));
print 'first point is in cell ', cell_id;

from time import clock
start = clock()
X = bm.where_is_points(Nx, Ny, c, z, P.reshape([1, -1]));
elapsed = (clock() - start)
print '==============================================='
print 'where_is_point time %.3f seconds' % elapsed
print '==============================================='

# write result to file
np.savetxt('points.txt', X.reshape([-1, 2]), fmt = "%10d")

