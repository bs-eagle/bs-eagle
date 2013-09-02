/**
	\file mpi_mesh_grdecl.cpp
 	\brief This file implement class for decompose main mesh among processes
	\author Iskhakov Ruslan
	\date 2008-10-14*/
#include "bs_mesh_stdafx.h"

#include "mpi_mesh_grdecl.h"
#ifdef _MPI_MY
namespace blue_sky
  {
  //! default constructor
  
  mpi_mesh_grdecl::mpi_mesh_grdecl(bs_type_ctor_param)
  {
  }
  
  mpi_mesh_grdecl::mpi_mesh_grdecl(const mpi_mesh_grdecl& src)
  {
    *this = src;
  }

  
  int mpi_mesh_grdecl::par_distribution()
  {
    return 0;
  }

  
  int mpi_mesh_grdecl::print_debug_info()
  {
    printf("numprocs=%d,myid = %d\n",numprocs,myid);
    fflush(0);
    printf("lx=%d, ly=%d, lz=%d\n",lx,ly,lz);
    fflush(0);
    printf("start_x=%d, start_y=%d, start_z=%d\n",i_start,j_start,k_start);
    fflush(0);
    printf("end_x=%d, end_y=%d, end_z=%d\n",i_end,j_end,k_end);
    fflush(0);
    return 0;
  }

  
  int mpi_mesh_grdecl::par_fileOpen_cube_hdf5(const char* file_name)
  {
    hid_t  file_id, dset_id = 0; //file and dataset identifiers
    int cut_dir, n_layers;
    hid_t  filespace = 0, memspace=0; //file and memory dataspace identifiers
    hsize_t  dims_memory;       // dataset dimensions

    //hyperslab selection parameters
    hsize_t stride[1];
    hsize_t block[1];
    hsize_t count[1];
    hsize_t start[1] = {0};
    //bounds for hyperslab selection
    i_start=0, j_start=0, k_start=0, i_end=0, j_end=0, k_end=0;
    file_id = open_file_mpi_hdf5(file_name); // open file
    if (file_id < 0 )
      return -1;

    //! read initial_data
    int initial_data[3];
    dims_memory = 3;
    count[0] = 3;
    stride[0] = 1;
    block[0] = 1;

    mpi_hdf5_read_dataset ("/grid/initial_data", file_id, dims_memory, stride, block, count, TYPE_INT, initial_data, 0,
                           nx, ny, i_start, j_start, k_start, i_end, j_end, k_end, 1);

    nx = initial_data[0];
    ny = initial_data[1];
    nz = initial_data[2];

    if (!myid)
      printf("nx %d, ny %d, nz %d\n", nx, ny, nz);

    //! make 1D-decomposition
    if (nx <= ny)
      cut_dir = CUT_DIR_X;
    else
      cut_dir = CUT_DIR_Y;
    //cut_dir = CUT_DIR_Z;
    if (!myid)
      printf("cut_dir %d\n", cut_dir);

    // at first, init as full
    i_start = 0;
    j_start = 0;
    k_start = 0;
    i_end = nx - 1;
    j_end = ny - 1;
    k_end = nz - 1;

    if (cut_dir == CUT_DIR_X)
      {
        i_start = (int) (myid) * ((double) nx / numprocs);
        i_end = (int) (myid + 1) * ((double) nx / numprocs) - 1;
        n_layers = nx;
      }
    else if (cut_dir == CUT_DIR_Y)
      {
        j_start = (int) myid * ((double) ny / numprocs);
        j_end = (int) (myid + 1) * ((double) ny / numprocs) - 1;
        n_layers = ny;
      }
    else if (cut_dir == CUT_DIR_Z)
      {
        k_start = (int) myid * ((double) nz / numprocs);
        k_end = (int) (myid + 1) * ((double) nz / numprocs) - 1;
        n_layers = nz;
      }

    if (n_layers < numprocs)
      {
        printf("Error: n_procs %d > nx %d!\n", numprocs, nx);
        return -1;
      }
    printf ("BOUND_INDEXES : [%d, %d] , [%d, %d] , [%d, %d]\n", i_start, i_end, j_start, j_end, k_start, k_end);

    if ((i_end < i_start) || (j_end < j_start) || (k_end < k_start))
      {
        printf("Wrong decomp!\n");
        return -1;
      }

    /*
      //! make 1D-decomposition
      if (nx <= ny)
        cut_dir = CUT_DIR_X;
      else
        cut_dir = CUT_DIR_Y;
      if (!myid)
        {printf("cut_dir %d\n", cut_dir); fflush(0);}

      //  at first, init as full
      i_start = 0;
      j_start = 0;
      k_start = 0;
      i_end = nx - 1;
      j_end = ny - 1;
      k_end = nz - 1;

      if (cut_dir == CUT_DIR_X)
        {
          i_start = (int) myid * ((double) nx / numprocs);
          i_end = (int) (myid + 1) * ((double) nx / numprocs) - 1;
          n_layers = nx;
        }
      else if (cut_dir == CUT_DIR_Y)
        {
          j_start = (int) myid * ((double) ny / numprocs);
          j_end = (int) (myid + 1) * ((double) ny / numprocs) - 1;
          n_layers = ny;
        }
      
        //
      else if (cut_dir == CUT_DIR_Z)
        {
          k_start = (int) mpi_rank * ((double) nz / mpi_size);
          k_end = (int) (mpi_rank + 1) * ((double) nz / mpi_size) - 1;
          n_layers = nz;
        }
      if (n_layers < numprocs)
        {
          printf("Error: n_procs %d > nx %d!\n", numprocs, nx);
          return -1;
        }
      printf ("BOUND_INDEXES : [%d, %d] , [%d, %d] , [%d, %d]\n", i_start, i_end, j_start, j_end, k_start, k_end);

      //share between process by equal part of active cells (begin)
      #ifdef _WIN32
        mesh_grdecl::fileOpen_cube_hdf5(file_name,i_start, j_start, k_start, i_end, j_end,k_end);
      //#elif UNIX //hdf5-parallel only in UNIX-like system
        //put misha direction he

      #endif
      printf("mesh_zcorn_number = %d\n",zcorn_array.size());
      return true;
    */
    return 0;
  }


  
  void mpi_mesh_grdecl::set_mpi_comm(MPI_Comm acomm)
  {
    mpi_comm = acomm;
    MPI_Comm_size(mpi_comm, &numprocs);
    MPI_Comm_rank(mpi_comm, &myid);
  }

  
  int mpi_mesh_grdecl::open_file_mpi_hdf5(const char *file_name)
  {
    //Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS); // Creates a new property as an instance of a property li_startt class.

    //Each process of the MPI communicator creates an access template and sets i_end up wi_endh MPI parallel access information. Thi_start i_start2 done wi_endh the H5Pcreate call to obtain the file access property li_startt and the H5Pset_fapl_mpio call to set up parallel I/O access.
    H5Pset_fapl_mpio(plist_id, mpi_comm, mpi_info);

    // Open a new file and release property list identifier
    hid_t file_id = H5Fopen (file_name, H5F_ACC_RDONLY, H5P_DEFAULT); // H5F_ACC_RDONLY - Allow read-only access to file.
    if (file_id < 0)
      {
        printf("Error: Can't open file %s!\n", file_name);
        return -1;
      }
    H5Pclose (plist_id);
    return file_id;
  }
  
  int mpi_mesh_grdecl::mpi_hdf5_read_dataset (const char* dset_title, hid_t file_id, hsize_t dims_memory,
      hsize_t *stride,
      hsize_t *block,
      hsize_t *count,
      const int datatype, int *data_int, float *data_float,
      const int nx, const int ny,
      const int i_start2, const int j_start, const int k_start,
      const int i_end, const int j_end, const int k_end, int mode)
  {
    herr_t status;
    int rank = 1; // read 1-dimensional array
    hid_t  filespace, memspace; // file and memory dataspace identifiers

    int k;
    float *data_float_ptr = data_float;
    int *data_int_ptr = data_int;
    hid_t dset_id;
    hsize_t start[1];

    if (datatype != TYPE_FLOAT && datatype != TYPE_INT)
      {
        printf("Error: wrong type!\n");
        return -1;
      }

    dset_id = open_dataset(rank, &dims_memory, file_id, dset_title);
    if (dset_id < 0)
      return -1; //cant' open dataset!

    filespace = H5Dget_space(dset_id);
    memspace = H5Screate_simple(rank, &dims_memory, NULL);
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    if (mode) // for 1 pass read
      {
        if (mode == 1)// initial_data
          start[0] = 0;
        else if (mode == 2)// coord
          start[0] = (i_start + j_start * (nx + 1)) * 6;

        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, stride, count, block);
        if (datatype == TYPE_FLOAT)
          status = H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data_float);
        else if (datatype == TYPE_INT)
          status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data_int);
      }
    else // multi pass read
      {
        for (k = 0; k < k_end - k_start + 1; k++)
          {
            start[0] = i_start + j_start * nx + (k + k_start) * nx * ny;
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, stride, count, block);

            if (datatype == TYPE_FLOAT)
              status = H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data_float_ptr);
            else if (datatype == TYPE_INT)
              status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data_int_ptr);

            if (status < 0)
              {
                printf("Error: Dataset read error!\n");
                return -1;
              }
            data_int_ptr = data_int_ptr + dims_memory;
            data_float_ptr = data_float_ptr + dims_memory;
          }
      }
    //release resources
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);

    return 0;
  }
  
  int mpi_mesh_grdecl::open_dataset(int rank, hsize_t dims_memory[], hid_t file_id, const char* dset_title)
  {
    hid_t filespace = H5Screate_simple(rank, dims_memory, NULL);
    //Open the dataset and close filespace
    hid_t dset_id = H5Dopen(file_id, dset_title);
    if (dset_id < 0)// check if file has been opened
      {
        printf("Error: Can't open dataset %s!\n", dset_title);
        return -1;
      }
    H5Sclose(filespace);
    return dset_id;
  }

  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (mpi_mesh_grdecl<base_strategy_fif>) , 1, (objbase), "mpi_mesh_grdecl<float, int,float>", "MPI mesh_grdecl class", "MPI mesh ecllipse class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (mpi_mesh_grdecl<base_strategy_did>) , 1, (objbase), "mpi_mesh_grdecl<double, int,double>", "MPI mesh ecllipse class", "MPI mesh ecllipse class", false);

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(mpi_mesh_grdecl, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(mpi_mesh_grdecl, (class));


}; //namespace blue_sky
#endif //_MPI_MY 
