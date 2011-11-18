#!/usr/bin/python
import bs
import numpy as np
import bs.bs_mesh as bm
import bs.sql_well as sw

# generate initial mesh
cell_x = 50; cell_y = 50; cell_z = 2;
dx = np.array([cell_x], dtype='d') # Lx
dy = np.array([cell_y], dtype='d') # Ly
dz = np.array([cell_z], dtype='d') # Lz
Nx = 100
Ny = 100
Nz = 100
[c, z] = bm.mesh_grdecl.gen_coord_zcorn(Nx, Ny, Nz, dx, dy, dz, 95000, 110000, 2200);
print 'Source mesh generated!'
pool = bs.bs_hdf5_storage.h5_pool()
pool.open_file("h5pool.h5", "/pool")
pool.py_set_pool_dims([Nx, Ny, Nz])
pool.set_fp_data ("COORD", c, 0)
pool.set_fp_data ("ZCORN", z, 0)

# open sql_well db
s = sw.sql_well();
s.open_db('frac.db');
s.save_to_bos_ascii_file ("4.out", pool)
s.close_db ()
pool.close_file ()

