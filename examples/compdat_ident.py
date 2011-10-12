#!/usr/bin/python
import bs
import numpy as np
import bs.bs_mesh as bm
import bs.sql_well as sw

# generate initial mesh
cell_x = 50; cell_y = 50; cell_z = 50;
dx = np.array([cell_x], dtype='d') # Lx
dy = np.array([cell_y], dtype='d') # Ly
dz = np.array([cell_z], dtype='d') # Lz
Nx = 100
Ny = 100
Nz = 50
[c, z] = bm.mesh_grdecl.gen_coord_zcorn(Nx, Ny, Nz, dx, dy, dz, 95000, 110000, 0);
print 'Source mesh generated!'

# open sql_well db
s = sw.sql_well();
s.open_db('db.db');
sw.completions_ident(s, 33160, Nx, Ny, c, z);

