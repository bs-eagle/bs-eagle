#!/usr/bin/python
import bs
import numpy as np
import bs.bs_mesh as bm

# generate initial mesh
cell_x = 25; cell_y = 25; cell_z = 1;
dx = np.array([cell_x], dtype='d') # Lx
dy = np.array([cell_y], dtype='d') # Ly
dz = np.array([cell_z], dtype='d') # Lz
Nx = 500
Ny = 500
Nz = 100
[c, z] = bm.mesh_grdecl.gen_coord_zcorn(Nx, Ny, Nz, dx, dy, dz)
print 'Source mesh generated!'

# generate randomly mesh trajectory
well_nodes_num = 5000;
wx = np.random.rand(well_nodes_num) * cell_x * Nx;
wy = np.random.rand(well_nodes_num) * cell_y * Ny;
wz = np.random.rand(well_nodes_num) * cell_z * Nz;
print wz

# calc MD
l = 0;
md = np.zeros(well_nodes_num)
for i in range(1, well_nodes_num) :
	l += np.sqrt((wx[i] - wx[i - 1])*(wx[i] - wx[i - 1]) + (wy[i] - wy[i - 1])*(wy[i] - wy[i - 1]) + (wz[i] - wz[i - 1])*(wz[i] - wz[i - 1]));
	md[i] = l;

# build full well specification
W = np.c_[wx, wy, wz, md];
print W

# calc mesh and well intersection
from time import clock
start = clock()
X = bm.well_path_ident(Nx, Ny, c, z, W.reshape([1, -1]), True)
elapsed = (clock() - start)
print '==============================================='
print 'Well_path_ident time %.3f seconds' % elapsed
print 'Result contain %d intersections' % (X.size / 7)
print '==============================================='

# write result to file
np.savetxt('x.txt', X.reshape([-1, 7]), fmt = "%8.4f")

