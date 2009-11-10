#ifndef MPI_MESH_H
#define MPI_MESH_H
#ifdef _MPI_MY
#include "mpi_tools.h"
#include "mpi_vector.h"
#include "mpi_csr_comm.h"
#include "mesh_grdecl.h"
#include "mpi.h"
/**
	\file mpi_mesh_grdecl.h
 	\brief This file declare class for decompose main mesh among processes
	\author Iskhakov Ruslan
	\date 2008-10-14
 */

namespace blue_sky
  {

  template <class strategy_t>
  class BS_API_PLUGIN mpi_mesh_grdecl: public  mesh_grdecl<strategy_t> //public mpi::mpi_user
    {
//+++++++++++++++++++++++++++++++++++++++++++
//  INTERNAL TYPE DECLARATION
//===========================================


      enum
      {
        TYPE_INT = 0,
        TYPE_FLOAT
      };

      enum
      {
        CUT_DIR_X = 0,
        CUT_DIR_Y,
        CUT_DIR_Z
      };

    public:
      //typedef smart_ptr<mesh_grdecl<base_strategy_fi>, true> sp_mesh_grdecl_fi;
      typedef mpi_mesh_grdecl<strategy_t> this_t;

//-------------------------------------------
//  METHODS
//===========================================
    public:
      //! blue-sky class declaration
      BLUE_SKY_TYPE_DECL(mpi_mesh_grdecl);
      //! default destructor
      ~mpi_mesh_grdecl () {};
      //! print debug info
      int print_debug_info();
      //! set MPI_Comm
      void set_mpi_comm(MPI_Comm);
      //! set MPI_INFO
      void set_mpi_info(MPI_Info ainfo)
      {
        mpi_info = ainfo;
      }
      //! open hdf5-file and share it between meshes
      int par_fileOpen_cube_hdf5(const char* file_name);

    protected:
      //! share mesh for equal size
      int par_distribution();
    private:
      //! open hdf5-file for mpi-reading
      int open_file_mpi_hdf5(const char *file_name);

      /** \brief read dataset for hdf5-file with mpi approach
          \param dset_title - name of dataset block (for instance "grid","actnum"...)
          \param file_id - id of hdf5-file
          \param dims_memory - size of dataset block which we reading
          \param stride - distance between blocks
          \param block - size of block
          \param count - array of counters, how many blocks necessary to read
          \param datatype - type of data
          \param data_int, data_float - buffer for writing
          \param nx,ny - size of reading area
          \param i_start, j_start, k_start, i_end, j_end, k_end - area of block
          \param mode - single/multi thread ????
      */
      int mpi_hdf5_read_dataset (const char* dset_title, hid_t file_id, hsize_t dims_memory,
                                 hsize_t *stride,
                                 hsize_t *block,
                                 hsize_t *count,
                                 const int datatype, int *data_int, float *data_float,
                                 const int nx, const int ny,
                                 const int i_start2, const int j_start, const int k_start,
                                 const int i_end, const int j_end, const int k_end, int mode);
      //open dataset for reading
      int open_dataset(int rank, hsize_t dims_memory[], hid_t file_id, const char* dset_title);
    protected:
//-------------------------------------------
//  VARIABLES
//===========================================
      //! MPI variables
      int myid, numprocs;
      MPI_Comm mpi_comm;
      MPI_Info mpi_info;

      //! parameters of local mesh
      index_t lx, ly, lz; //!< number of blocks on X&Y&Z in local mesh
      index_t i_start, j_start, k_start; //!< number of start block on X&Y&Z in global mesh for local mesh
      index_t i_end, j_end, k_end; //!< number of end block on X&Y&Z in global mesh for local mesh

      index_t n_layers;
      //! 1-dimension decomposition params
      int cut_dir;

      //! ghost_cells (bool-array)

    }; //class mpi_mesh_grdecl
}; //namespace blue_sky
#endif //_MPI_MY
#endif //MPI_MESH_H
