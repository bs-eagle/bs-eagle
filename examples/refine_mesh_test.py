#!/usr/bin/python
import bs
import numpy
import bs.bs_mesh as bm

# generate initial mesh
dx = numpy.array([5], dtype='d') # Lx
dy = numpy.array([5], dtype='d') # Ly
dz = numpy.array([1], dtype='d') # Lz
Nx = 10
Ny = 10
Nz = 1
[c, z] = bm.mesh_grdecl.gen_coord_zcorn(Nx, Ny, Nz, dx, dy, dz)
print 'Source mesh generated!'
print 'Starting refine...'

# generate begin, middle & end coordinates of fault
Wx = 0.05
Wy = 0.05
ax = 2
ay = 2
vx = numpy.array([10,10]) #numpy.random.rand(1) * 50;
vy = numpy.array([10,10]) #numpy.random.rand(1) * 50;
vax = numpy.zeros(vx.size);
vax.fill(ax);
vay = numpy.zeros(vy.size);
vay.fill(ay);
vWx = numpy.array([0.05, 0.5]) #numpy.random.rand(vx.size) * 10;
vWy = numpy.array([0.05, 0.5]) #numpy.random.rand(vy.size) * 10;
p = numpy.reshape(numpy.column_stack([vx, vy, vWx, vWy, vax, vay]), [1, -1]);
#p = numpy.array([
#	1000, 5000, Wx, Wy, ax, ay,
#	5000, 5000, Wx, Wy, ax, ay,
#	9000, 5000, Wx, Wy, ax, ay,
#	], dtype='d')
print p
# generate refined mesh
[c1, z1, Nx, Ny, hit_idx] = bm.mesh_grdecl.refine_mesh(Nx, Ny, c, z, p)
print 'refine finisfed!'
print 'Nx = ', Nx
print 'Ny = ', Ny
print 'model size = ', Nx * Ny * Nz
print 'coord.size = ', c1.size
print 'zcorn.size = ', z1.size

# save to files
print 'Write COORD & ZCORN to text file'
f = open('model_grid.inc', 'w');
f.write('COORD\n')
numpy.savetxt(f, c1.reshape([-1, 6]), fmt="%8.4f");
f.write('/\n\nZCORN\n')
numpy.savetxt(f, z1.reshape([-1, 2]), fmt="%8.4f");
f.write('/')
f.close();

