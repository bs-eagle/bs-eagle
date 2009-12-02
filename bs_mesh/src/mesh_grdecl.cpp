/*! \file mesh_grdecl.cpp
\brief This file implement class for working with grid_ecllipse
\author Iskhakov Ruslan
\date 2008-05-20 */
#include "bs_mesh_stdafx.h"
#include "mesh_grdecl.h"

using namespace grd_ecl;

using namespace rs_mesh_detail;
using namespace blue_sky;

const char filename_hdf5[] = "d:\\grid_swap.h5";

#ifdef _HDF5_MY
template<class strategy_t>
int mesh_grdecl<strategy_t>::create_array_hdf5(const char *dataset_name, H5::H5File &file_hdf5, H5::DataSet **dataset)
{
  // creating dataset_coords
  hsize_t dims[] = {0};
  hsize_t dims_max[] = {H5S_UNLIMITED};
  H5::DataSpace dataspace_coords(1, dims, dims_max);
  // set the dataset to be chunked
  H5::DSetCreatPropList cparms;
  hsize_t chunk_dims[] ={1};
  cparms.setChunk(1, chunk_dims);
  *dataset = new H5::DataSet(file_hdf5.createDataSet(dataset_name, H5::PredType::NATIVE_FLOAT, dataspace_coords, cparms));
  return 0;
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_activs_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt)
{
  return true;
  /*
  // turn off error printing
  H5E_auto2_t old_func;
  void *old_client_data;
  H5::Exception::getAutoPrint(old_func, &old_client_data);
  H5::Exception::dontPrint(); // turn off error printing

  // check if file exists
  hid_t file_hid;
  file_hid = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_hid < 0)
  return false;
  H5Fclose(file_hid);

  H5::H5File file(file_name, H5F_ACC_RDWR);

  try
  {
  actnum_array.resize((it - is + 1) * (jt - js  + 1) * (kt - ks + 1));
  H5::DataSet dataset = file.openDataSet("activ");
  hsize_t dims_memory[] = {((it-is+1))*(jt - js + 1)};
  H5::DataSpace dataspace_memory(1, dims_memory);
  H5::DataSpace dataspace_file = dataset.getSpace();
  hsize_t count[] = {jt - js + 1};
  hsize_t stride[] = {nx};
  hsize_t block[] = {it - is + 1};
  for (int k = 0; k < kt - ks + 1; k++)
  {
  hsize_t start[] = {is + js * nx + (k + ks) * nx * ny};
  dataspace_file.selectHyperslab(H5S_SELECT_SET, count, start, stride, block);
  dataset.read(&actnum_array[k * count[0] * block[0]], H5::PredType::NATIVE_INT, dataspace_memory, dataspace_file);
  }
  }
  catch (H5::FileIException exception)
  {
  exception.printError();
  return false;
  }

  H5::Exception::setAutoPrint(old_func, old_client_data);
  n_active_elements = accumulate(actnum_array.begin(),actnum_array.end(),0);
  return true;
  */
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_cube_hdf5(const char* file_name, int is, int js, int ks, int it, int jt, int kt)
{
  // turn off error printing
  H5E_auto2_t old_func;
  void *old_client_data;
  H5::Exception::getAutoPrint(old_func, &old_client_data);
  H5::Exception::dontPrint();

  // check if file exists
  hid_t file_hid;
  file_hid = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_hid < 0)
    return false;
  H5Fclose(file_hid);

  H5::H5File file(file_name, H5F_ACC_RDONLY);

  try
    {
      // read "coord" array
      H5::DataSet dataset_coord = file.openDataSet("coord");
      hsize_t dims_memory_coord[] = {6 * (it - is + 2) * (jt - js  + 2)};
      H5::DataSpace dataspace_memory_coord(1, dims_memory_coord);
      H5::DataSpace dataspace_file_coord = dataset_coord.getSpace();
      coord_array.clear();
      vector<item_t> buf_array(6 * (it - is + 2) * (jt - js  + 2));
      // hyperslab settings
      hsize_t count_coord[] = {jt - js + 2};
      hsize_t start_coord[] = {(is + js * (nx + 1)) * 6};
      hsize_t stride_coord[] = {6 * (nx + 1)};
      hsize_t block_coord[] = {(it - is + 2) * 6};
      // we read "coord" array by 1 read operation
      dataspace_file_coord.selectHyperslab(H5S_SELECT_SET, count_coord, start_coord, stride_coord, block_coord);
      dataset_coord.read(&buf_array[0], H5::PredType::NATIVE_FLOAT, dataspace_memory_coord, dataspace_file_coord);
      for (int i = 0; i < (it - is + 2) * (jt - js  + 2); i++)
        coord_array.push_back(coordElem (fpoint3d(buf_array[i*6],buf_array[i*6+1],buf_array[i*6+2]),
                                         fpoint3d(buf_array[i*6+3],buf_array[i*6+4],buf_array[i*6+5])));

      // read "zcorn" array
      H5::DataSet dataset_zcorn = file.openDataSet("zcorn");
      hsize_t dims_memory_zcorn[] = {8 * (it - is + 1) * (jt - js + 1)};
      H5::DataSpace dataspace_memory_zcorn(1, dims_memory_zcorn);
      H5::DataSpace dataspace_file_zcorn = dataset_zcorn.getSpace();
      zcorn_array.resize(8 * (it - is + 1) * (jt - js  + 1) * (kt - ks + 1));
      // hyperslab settings (common for each plane)
      hsize_t count_zcorn[] = {4 * (jt - js + 1)};
      hsize_t stride_zcorn[] = {2 * nx};
      hsize_t block_zcorn[] = {2 * (it - is + 1)};
      // we read "zcorn" array by planes - kt-ks+1 read operations
      for (int k = 0; k < kt - ks + 1; k++)
        {
          // determine start array for each plane individually
          hsize_t start_zcorn[] = {2 * is + 4 * js * nx + 8 * k * nx * ny};
          dataspace_file_zcorn.selectHyperslab(H5S_SELECT_SET, count_zcorn, start_zcorn, stride_zcorn, block_zcorn);
          dataset_zcorn.read(&zcorn_array[k * count_zcorn[0] * block_zcorn[0]], H5::PredType::NATIVE_FLOAT, dataspace_memory_zcorn, dataspace_file_zcorn);
        }
    }
  catch (H5::Exception exception)
    {
      exception.printError();
      coord_array.clear();
      zcorn_array.clear();
      return false;
    }

  // enable error printing
  H5::Exception::setAutoPrint(old_func, old_client_data);

  return true;
}

template<class strategy_t>
int mesh_grdecl<strategy_t>::append_array_hdf5(const item_t *arr, size_t arr_length, H5::DataSet *dataset)
{
  // determine new dimensions of dataset
  hsize_t dims_old[1];
  hsize_t dims_memory[1] = {arr_length};
  H5::DataSpace dataspace_file = dataset->getSpace();
  dataspace_file.getSimpleExtentDims(dims_old);
  hsize_t dims_new[1] = {dims_old[0] + dims_memory[0]};
  // extend dataset
  dataset->extend(dims_new);
  dataspace_file = dataset->getSpace();
  H5::DataSpace dataspace_memory(1, dims_memory);
  // select hyperslab to write
  hsize_t count[1] = {dims_memory[0]};
  hsize_t start[1] = {dims_old[0]};
  dataspace_file.selectHyperslab(H5S_SELECT_SET, count, start);
  // write array
  dataset->write(arr, H5::PredType::NATIVE_FLOAT, dataspace_memory, dataspace_file);

  return 0;
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_cube_with_hdf5_swap(const char* file_name)
{
  try
    {
      fstream file(file_name, ios::in);
      if (!file.is_open())
        return false;

      H5::H5File file_hdf5(filename_hdf5, H5F_ACC_TRUNC);
      H5::DataSet *dataset_coords = 0;
      H5::DataSet *dataset_zcorn = 0;

      char buf[100];
      item_array_t buf_array;
      bool array_end;

      /////////////////////////////////
      // read "coord"
      /////////////////////////////////
      while (!file.eof())
        {
          file >> buf;
          if (strcmp(buf, "COORD") == 0)
            break;
        }
      create_array_hdf5("coord", file_hdf5, &dataset_coords);
      array_end = false;
      while (!file.eof())
        {
          // check if it is array end
          file >> buf;
          if (buf[0] == '/')
            array_end = true;

          if (array_end != true)
            {
              buf_array.push_back((float) atof(buf));
            }

          // if maximum space is used
          // or if it is array end
          // write "coords" to the swap file
          if ((buf_array.size() >= HDF5_MAX_SPACE)  ||
              ((array_end == true) && (buf_array.size() > 0))
             )
            {
              // save coord_array to hdf5 swap file
              append_array_hdf5(&buf_array[0], buf_array.size(), dataset_coords);
              buf_array.clear();
            }

          if (array_end == true)
            break;
        }
      delete dataset_coords;

      /////////////////////////////////
      // read "zcorn"
      /////////////////////////////////
      while (!file.eof())
        {
          file >> buf;

          if (strcmp(buf, "ZCORN") == 0)
            break;
        }
      // creating dataset_zcorn
      create_array_hdf5("zcorn", file_hdf5, &dataset_zcorn);
      array_end = false;
      while (!file.eof())
        {
          // check if it is array end
          file >> buf;
          if (buf[0] == '/')
            array_end = true;

          if (array_end != true)
            {
              buf_array.push_back((float) atof(buf));
            }

          // if maximum space is used
          // or if it is array end
          // write "coords" to the swap file
          if ((buf_array.size() >= HDF5_MAX_SPACE)  ||
              ((array_end == true) && (buf_array.size() > 0))
             )
            {
              // save coord_array to hdf5 swap file
              append_array_hdf5(&buf_array[0], buf_array.size(), dataset_zcorn);
              buf_array.clear();
            }

          if (array_end == true)
            break;
        }
      delete dataset_zcorn;
      return true;
    }
  catch (H5::Exception exception)
    {
      exception.printError();
      return false;
    }
}
#endif

template<class strategy_t>
fpoint3d  mesh_grdecl<strategy_t>::cross_coord(const item_t z, const pool_item_t *coord)const
  {
    fpoint3d p;
    p.z = z;
    /*
    float temp = (p.z-m_Coord.pStart.z)/(m_Coord.pEnd.z-m_Coord.pStart.z);
    p.x = temp *(m_Coord.pEnd.x-m_Coord.pStart.x)+m_Coord.pStart.x;
    p.y = temp *(m_Coord.pEnd.y-m_Coord.pStart.y)+m_Coord.pStart.y;
    */
    float temp = (p.z - coord[2]) / (coord[5] - coord[2]);
    p.x = temp *(coord[3] - coord[0]) + coord[0];
    p.y = temp *(coord[4] - coord[1]) + coord[1];
    return p;
  }
template<class strategy_t>
boost::array <fpoint3d, 8> mesh_grdecl<strategy_t>::top_cube(const index_t i, const index_t j, const index_t k) const
  {
    g_fpoint3d_vector cubeVertex;


#ifdef _DEBUG
    BS_ASSERT ( (i >= 0) && (i < nx) && (j >= 0) && (j < ny) && (k >= 0) && (k < nz));
#endif

    //define index
    index_t index1 = i*2+j*4*nx+k*8*nx*ny;//upper side
    index_t index2 = index1+4*nx*ny;//lower side
    index_t iCOORD = i+j*(nx+1);

    /*
    cubeVertex[0] = cross_coord(zcorn_array[index1],coord_array[iCOORD]);
    cubeVertex[1] = cross_coord(zcorn_array[index1+1],coord_array[iCOORD+1]);
    cubeVertex[2] = cross_coord(zcorn_array[index1+2*nx],coord_array[iCOORD+(nx+1)]);
    cubeVertex[3] = cross_coord(zcorn_array[index1+2*nx+1],coord_array[iCOORD+(nx+1)+1]);

    cubeVertex[4] = cross_coord(zcorn_array[index2],coord_array[iCOORD]);
    cubeVertex[5] = cross_coord(zcorn_array[index2+1],coord_array[iCOORD+1]);
    cubeVertex[6] = cross_coord(zcorn_array[index2+2*nx],coord_array[iCOORD+(nx+1)]);
    cubeVertex[7] = cross_coord(zcorn_array[index2+2*nx+1],coord_array[iCOORD+(nx+1)+1]);
    */
    
    cubeVertex[0] = cross_coord(zcorn_array[index1], &coord_array[iCOORD * 6]);
    cubeVertex[1] = cross_coord(zcorn_array[index1+1], &coord_array[(iCOORD+1) * 6]);
    cubeVertex[2] = cross_coord(zcorn_array[index1+2*nx], &coord_array[(iCOORD+(nx+1)) * 6]);
    cubeVertex[3] = cross_coord(zcorn_array[index1+2*nx+1], &coord_array[(iCOORD+(nx+1)+1) * 6]);

    cubeVertex[4] = cross_coord(zcorn_array[index2], &coord_array[(iCOORD) * 6]);
    cubeVertex[5] = cross_coord(zcorn_array[index2+1], &coord_array[(iCOORD+1) * 6]);
    cubeVertex[6] = cross_coord(zcorn_array[index2+2*nx], &coord_array[(iCOORD+(nx+1)) * 6]);
    cubeVertex[7] = cross_coord(zcorn_array[index2+2*nx+1], &coord_array[(iCOORD+(nx+1)+1) * 6]);
    
    return cubeVertex;
  }
template<class strategy_t>
float mesh_grdecl<strategy_t>::get_volume_cube(const g_fpoint3d_vector &cube) const
  {
    fpoint3d center = get_cube_center (cube);
    item_t volume = 0.0;

    //share for 12 tetraidr (6 side -> 2 tetraidr for each)
    volume += grd_ecl::calc_tetra_volume(cube[0],cube[1],cube[2], center);
    volume += grd_ecl::calc_tetra_volume(cube[1],cube[2],cube[3], center);

    volume += grd_ecl::calc_tetra_volume(cube[1],cube[3],cube[7], center);
    volume += grd_ecl::calc_tetra_volume(cube[1],cube[5],cube[7], center);

    volume += grd_ecl::calc_tetra_volume(cube[0],cube[2],cube[6], center);
    volume += grd_ecl::calc_tetra_volume(cube[0],cube[4],cube[6], center);

    volume += grd_ecl::calc_tetra_volume(cube[4],cube[5],cube[6], center);
    volume += grd_ecl::calc_tetra_volume(cube[5],cube[6],cube[7], center);

    volume += grd_ecl::calc_tetra_volume(cube[0],cube[1],cube[4], center);
    volume += grd_ecl::calc_tetra_volume(cube[1],cube[4],cube[5], center);

    volume += grd_ecl::calc_tetra_volume(cube[2],cube[3],cube[6], center);
    volume += grd_ecl::calc_tetra_volume(cube[3],cube[6],cube[7], center);

    return volume;
  }


template<class strategy_t>
int mesh_grdecl<strategy_t>::define_side(const index_t i, const index_t j,const  index_t k, const index_t i1, const index_t j1, const index_t k1) const
  {
    if (j == j1) //one row
      {
        if (i == i1+1)
          return  left;
        if (i == i1-1)
          return  right;
      }
    if (i == i1) //one column
      {
        if (j == j1+1)
          return  top;
        if (j == j1-1)
          return  bottom;
      }
    if (k == k1+1)
      return upper;
    if (k == k1-1)
      return lower;
    return empty;
  }

template<class strategy_t>
bool mesh_grdecl<strategy_t>::is_small(const index_t i, const index_t j, const index_t k, item_t eps)  const
  {
    if (k >= nz)
      return false;

    item_t dz1, dz2, dz3, dz4; //height for each coord
    //define index
    index_t index1 = i*2+j*4*nx+k*8*nx*ny;	//lower side
    index_t index2 = index1+4*nx*ny;			//upper side
    dz1 = zcorn_array[index2]-zcorn_array[index1];
    dz2 = zcorn_array[index2+1]-zcorn_array[index1+1];
    dz3 = zcorn_array[index2+2*nx]-zcorn_array[index1+2*nx];
    dz4 = zcorn_array[index2+2*nx+1]-zcorn_array[index1+2*nx+1];

    if (dz1 <= eps && dz2 <= eps && dz3 <= eps && dz4 <= eps)
      return true;
    return false;
  }

template<class strategy_t>
typename mesh_grdecl<strategy_t>::item_t mesh_grdecl<strategy_t>::calc_block_volume(const index_t i, const index_t j, const index_t k) const
  {
    const g_fpoint3d_vector &cube = top_cube(i, j, k);
    index_t index = i + j * nx + k * nx * ny;
    if (sp_ntg.empty ())
      return (get_volume_cube(cube)); // * (*sp_poro)[index]);
    else
      return (get_volume_cube(cube) * sp_ntg[index]);// * (*sp_poro)[index]);
  }

template<class strategy_t>
fpoint3d mesh_grdecl<strategy_t>::get_side_crossing_projection_on_all_axis(const point_side_t &side1, const point_side_t &side2)const
  {
    fpoint3d sqArea;
    g_fpoint2d_vector polyg1, polyg2;

    //project on YoZ
    for (size_t i = 0; i < side1.size(); i++)
      {
        polyg1.push_back(fpoint2d(side1[i].y, side1[i].z));
        polyg2.push_back(fpoint2d(side2[i].y, side2[i].z));
      }
    sqArea.x = get_polygon_crossing_area(polyg1, polyg2, EPS_SQ);

    //project on XoZ
    polyg1.clear();
    polyg2.clear();
    for (size_t i = 0; i < side1.size(); i++)
      {
        polyg1.push_back(fpoint2d(side1[i].x, side1[i].z));
        polyg2.push_back(fpoint2d(side2[i].x, side2[i].z));
      }
    sqArea.y = get_polygon_crossing_area(polyg1, polyg2, EPS_SQ);

    if (sqArea.x + sqArea.y == 0.0)
      return fpoint3d(0,0,0);

    //project on XoY
    polyg1.clear();
    polyg2.clear();
    for (size_t i = 0; i < side1.size(); i++)
      {
        polyg1.push_back(fpoint2d(side1[i].x, side1[i].y));
        polyg2.push_back(fpoint2d(side2[i].x, side2[i].y));
      }
    sqArea.z = get_polygon_crossing_area(polyg1, polyg2, EPS_SQ);
    return sqArea;
  }

template <typename strategy_t>
typename mesh_grdecl<strategy_t>::center_t
mesh_grdecl<strategy_t>::get_center (index_t i_coord, index_t j_coord, index_t k_coord) const
{
  BS_ASSERT (i_coord != -1 && j_coord != -1 && k_coord != -1) (i_coord) (j_coord) (k_coord);

  grd_ecl::fpoint3d point (get_cube_center (top_cube (i_coord, j_coord, k_coord)));

  center_t res;
  res[0] = point.x;
  res[1] = point.y;
  res[2] = point.z;

  return res;
}

template <typename strategy_t>
typename mesh_grdecl<strategy_t>::center_t
mesh_grdecl<strategy_t>::get_center (index_t n_block) const
{
  BS_ASSERT (n_block != -1) (n_block);

  index_t i_coord = 0, j_coord = 0, k_coord = 0;
  inside_to_XYZ (n_block, i_coord, j_coord, k_coord);
  return get_center (i_coord, j_coord, k_coord);
}

template<class strategy_t>
typename mesh_grdecl<strategy_t>::item_t mesh_grdecl<strategy_t>::get_center_zcorn(const cube_index_t &cube)const
  {
    item_t res = 0.0;
    for (size_t i = 0; i < cube.size(); i++)
      res += zcorn_array[cube[i]];
    return res/cube.size();
  }

template<class strategy_t>
void mesh_grdecl<strategy_t>::change_row(const index_t i, const index_t j,const  index_t k, index_array_t &rows_ptr)
{
  index_t index = XYZ_to_inside(i,j,k);
  rows_ptr[index+1]++;
}
template<class strategy_t>
void mesh_grdecl<strategy_t>::change_col(const index_t index1, const index_t index2, index_array_t &cols_ind, index_array_t& curIndex)
{
  cols_ind[curIndex[index1]] = index2;
  curIndex[index1]++;
}

template<class index_t, typename index_array_t>
void 
change_col_add( const index_t index1, const index_t index2, index_array_t &cols_ind, index_array_t& curIndex,
    index_array_t &m_memory, index_array_t &p_memory, bool is_m_memory,
    bool is_need_to_add, const index_t connection_number, const index_array_t &rows_ptr)
{
  cols_ind[curIndex[index1]] = index2;

  if (is_need_to_add)
    {
      if (is_m_memory)
        {
          m_memory[connection_number] = rows_ptr[index1];
          m_memory[connection_number+1] = curIndex[index1];
        }
      else
        {
          p_memory[connection_number+1] = rows_ptr[index1];
          p_memory[connection_number] = curIndex[index1];
        }
    }
  ++curIndex[index1];
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::set_plus_minus (const index_t index1, const index_t index2, index_array_t &cols_ind,
    index_array_t &m_memory, index_array_t &p_memory, const index_t connection_number_minus,
    const index_t connection_number_plus, const index_array_t &rows_ptr)
{

  int i, i1, i2;

  i1 = rows_ptr[index1] + 1;
  i2 = rows_ptr[index1 + 1];

  //TODO: make search taking into account cols_ind is sorted
  for (i = i1; i < i2; ++i)
    {
      if (cols_ind[i] == index2)
        {
          m_memory[connection_number_minus] = i;
          break;
        }
    }

  i1 = rows_ptr[index2] + 1;
  i2 = rows_ptr[index2 + 1];

  for (i = i1; i < i2; ++i)
    {
      if (cols_ind[i] == index1)
        {
          p_memory[connection_number_plus] = i;
          break;
        }
    }
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_actnum(const char* file_name)
{
#if 0
  fstream file(file_name,  ios::in);
  char buf[255];
  char *start_ptr,*end_ptr;
  if (!file.is_open())
    {
      n_active_elements =  accumulate(sp_actnum.begin(),sp_actnum.end(),0);
      return false;
    }
  sp_actnum.clear();

  while (!file.eof())
    {
      file >> buf;
      start_ptr = buf;
      if (strcmp (buf, "/") != 0)
        sp_actnum.push_back((int)strtod (start_ptr, &end_ptr));
      else
        break;
    }
  n_active_elements =  accumulate(sp_actnum.begin(),sp_actnum.end(),0);
#endif
  return true;
}


template<class strategy_t>
bool mesh_grdecl<strategy_t>::file_open_cube(const char* file_name)
{
  /*
  using namespace std;
  fstream file(file_name,  ios::in);
  char buf[100];
  char *start_ptr, *end_ptr;
  if (!file.is_open())
    return false;

  while (!file.eof())
    {
      file >> buf;
      if (  strcmp (buf, "COORD") == 0)
        break;
    }

  max_x = max_y = -100000;
  min_x = min_y = 1000000;

  index_t i_count = 0;
  const int max_count = 6; //buffer_size
  vector<float> buf_vec(max_count);
  while (!file.eof())
    {
      //coordElem curElem;
      file >> buf;

      if (buf[0] != '/')
        {
          buf_vec[i_count++] = (float)atof(buf);
          if (i_count == max_count) //swap_buffer
            {
              i_count = 0;
              coordElem curElem = coordElem(fpoint3d(buf_vec[0],buf_vec[1],buf_vec[2]),
                                            fpoint3d(buf_vec[3],buf_vec[4],buf_vec[5]));
              coord_array.push_back(curElem);

              //looking for max&min coordinates
              if (curElem.pStart.x > max_x)   max_x = curElem.pStart.x;
              if (curElem.pEnd.x > max_x)     max_x = curElem.pEnd.x;

              if (curElem.pStart.y > max_y)   max_y = curElem.pStart.y;
              if (curElem.pEnd.y > max_y)     max_y = curElem.pEnd.y;

              if (curElem.pStart.x < min_x)   min_x = curElem.pStart.x;
              if (curElem.pEnd.x < min_x)     min_x = curElem.pEnd.x;

              if (curElem.pStart.y < min_y)   min_y = curElem.pStart.y;
              if (curElem.pEnd.y < min_y)     min_y = curElem.pEnd.y;
            }
        }
      else
        break;
    }

  while (!file.eof())
    {
      file >> buf;
      if (  strcmp (buf, "ZCORN") == 0)
        break;
    }

  while (!file.eof())
    {
      file >> buf;
      start_ptr = buf;
      if (strcmp (buf, "/") != 0)
        zcorn_array.push_back((float)strtod (start_ptr, &end_ptr));
      else
        break;
    }
  min_z = *(std::min_element(zcorn_array.begin(),zcorn_array.end()));
  max_z = *(std::max_element(zcorn_array.begin(),zcorn_array.end()));
  return true;
  */
  return false;
}

template<class strategy_t>
int mesh_grdecl<strategy_t>::init_ext_to_int()
{
  item_array_t volumes_temp(n_elements);

  //tools::save_seq_vector ("actnum.bs.txt").save sp_actnum;

  int splicing_num = splicing(volumes_temp);

  //tools::save_seq_vector ("active_blocks.bs.txt").save sp_actnum;


  //make proxy array
  ext_to_int.resize(n_elements,0);
  size_t n_count = 0;

  index_t nn_active = 0, i_index; //number of non-active previous cells
  for (index_t i = 0; i < nz; ++i)
    {
      for (index_t j = 0; j < ny; ++j)
        for (index_t k = 0; k < nx; ++k, ++n_count)
          {
            i_index = k + (nx * j) + (i * nx * ny);
            if (!sp_actnum[i_index])
              {
                nn_active++;
                ext_to_int[n_count] = -1;
              }
            else
              ext_to_int[n_count] = i_index-nn_active;
          }
    }

  //tools::save_seq_vector ("ext_to_int.bs.txt").save (ext_to_int);

  init_int_to_ext();

  //fill volume array (except non-active block and using proxy array)
  volumes.resize(n_active_elements);
  for (int i = 0; i < n_active_elements; ++i)
    volumes[i] = volumes_temp[int_to_ext[i]];

  calc_depths();
  return splicing_num;
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::splice_two_blocks (const index_t i, const index_t j, const index_t k, const index_t k1)
{
  index_t index, index1, index2;

  index1 = i * 2 + j * 4 * nx + k1 * 8 * nx * ny; //upper side of [i,j,k1]
  index2 = i * 2 + j * 4 * nx + (k1 * 8 + 4) * nx * ny; //lower side of [i,j,k1]
  if (k > k1)
    {
      index = i * 2 + j * 4 * nx + k * 8 * nx * ny; //upper side of [i,j,k]

      zcorn_array[index] = zcorn_array[index1];
      zcorn_array[index + 1] = zcorn_array[index1 + 1];
      zcorn_array[index + 2 * nx] = zcorn_array[index1 + 2 * nx];
      zcorn_array[index + 2 * nx + 1] = zcorn_array[index1 + 2 * nx + 1];

      zcorn_array[index2] = zcorn_array[index1];
      zcorn_array[index2 + 1] = zcorn_array[index1 + 1];
      zcorn_array[index2 + 2 * nx] = zcorn_array[index1 + 2 * nx];
      zcorn_array[index2 + 2 * nx + 1] = zcorn_array[index1 + 2 * nx + 1];
    }
  else
    {
      index = i * 2 + j * 4 * nx + (k * 8 + 4) * nx * ny; //lower side of [i,j,k]

      zcorn_array[index] = zcorn_array[index2];
      zcorn_array[index + 1] = zcorn_array[index2 + 1];
      zcorn_array[index + 2 * nx] = zcorn_array[index2 + 2 * nx];
      zcorn_array[index + 2 * nx + 1] = zcorn_array[index2 + 2 * nx + 1];

      zcorn_array[index1] = zcorn_array[index2];
      zcorn_array[index1 + 1] = zcorn_array[index2 + 1];
      zcorn_array[index1 + 2 * nx] = zcorn_array[index2 + 2 * nx];
      zcorn_array[index1 + 2 * nx + 1] = zcorn_array[index2 + 2 * nx + 1];
    }
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::are_two_blocks_close (const index_t i, const index_t j, const index_t k, const index_t k1)
{
  index_t index, index1;
  if (k > k1)
    {
      index = i * 2 + j * 4 * nx + k * 8 *nx * ny; //upper side of [i,j,k]
      index1 = i * 2 + j * 4 * nx + (k1 * 8 + 4) * nx * ny; //lower side of [i,j,k1]
    }
  else
    {
      index = i * 2 + j * 4 * nx + (k * 8 + 4) * nx * ny; //lower side of [i,j,k]
      index1 = i * 2 + j * 4 * nx + k1 * 8 * nx * ny; //upper side of [i,j,k1]
    }

  if ((fabs (zcorn_array[index1] - zcorn_array[index]) < max_thickness) &&
      (fabs (zcorn_array[index1 + 1] - zcorn_array[index + 1]) < max_thickness) &&
      (fabs (zcorn_array[index1 + 2 * nx] - zcorn_array[index + 2 * nx]) < max_thickness) &&
      (fabs (zcorn_array[index1 + 2 * nx + 1] - zcorn_array[index + 2 * nx + 1]) < max_thickness))
    {
      return true;
    }
  else
    {
      return false;
    }
}

template<class strategy_t>
int mesh_grdecl<strategy_t>::splicing(item_array_t& volumes_temp)
{
  index_t nCount = 0;
  item_t vol_sum, vol_block;
  index_t i, j, k, k1;
  index_t small_block_top, big_block_top;


#ifdef _DEBUG
    BS_ASSERT (zcorn_array.size() && coord_array.size());
#endif

  array_uint8_t &actnum     = sp_actnum;
  const array_float16_t &poro = sp_poro;
  
  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j)
      {
        small_block_top = -1;
        big_block_top = -1;
        vol_sum = 0.0;
        for (k = 0; k < nz; ++k)
          {
            index_t index = i + j * nx + k * nx * ny;
            // miss inactive blocks
            if (!actnum[index])
              {
                // blocks can`t be spliced across inactive blocks
                big_block_top = -1;
                small_block_top = -1;
                vol_sum = 0.0;
                continue;
              }

            vol_block = calc_block_volume(i,j,k);
            item_t vol_block_poro = vol_block * poro[index];

            if (vol_block_poro < minpv * 0.99)
              {
                // block is too small, set as inactive
                actnum[index] = 0;
                ++nCount;
                // blocks can`t be spliced across inactive blocks
                vol_sum = 0.0;
                big_block_top = -1;
                small_block_top = -1;
                continue;
              }
            else if (vol_block_poro < minsv)
              {
                // this is small block

                volumes_temp[index] = vol_block;

                if (big_block_top != -1)
                  {
                    // this block is absorbed by bigger block above
                    splice_two_blocks (i, j, big_block_top, k);
                    // add volume this small block to the volume of big block
                    volumes_temp[i + j * nx + big_block_top * nx * ny] += vol_block;

                    // make block inactive
                    actnum[index] = 0;
                    ++nCount;
                    small_block_top = -1;
                  }

                // check if this block is close enough to next block
                if ((k < (nz - 1)) && (are_two_blocks_close (i, j, k, k + 1)))
                  {
                    if (big_block_top == -1)
                      {
                        if (small_block_top == -1)
                          {
                            // if already not spliced, can be spliced with lower block
                            small_block_top = k;
                          }

                        // aggregate volume of small blocks in case they`ll be spliced with lower block
                        vol_sum += vol_block;
                      }
                  }
                else
                  {
                    // this block and next block are not close, can`t be coupled
                    small_block_top = -1;
                    big_block_top = -1;
                    vol_sum = 0.0;
                    // TODO: make multperm = 0
                  }

              }
            else
              {
                // this is big block

                volumes_temp[index] = vol_block;

                if (small_block_top != -1)
                  {
                    // this block absorbes all smaller blocks above
                    for (k1 = k - 1; k1 >= small_block_top; --k1)
                      {
                        splice_two_blocks (i, j, k, k1);
                        // make small block inactive
                        actnum[i + j * nx + k1 * nx * ny] = 0;
                        ++nCount;
                      }

                    // add volume all small blocks above to the volume of this big block
                    volumes_temp[index] += vol_sum;
                    vol_sum = 0.0;
                    small_block_top = -1;
                  }

                // check if this block is close enough to next block
                if ((k < (nz - 1)) && (are_two_blocks_close (i, j, k, k + 1)))
                  {
                    // this block can absorb lower small blocks
                    big_block_top = k;
                  }
                else
                  {
                    // this block and next block are not close, can`t be coupled
                    small_block_top = -1;
                    big_block_top = -1;
                    // TODO: make multperm = 0
                  }
              }
           }
        }
  n_active_elements -= nCount;
  return nCount;
}


template<class strategy_t>
int mesh_grdecl<strategy_t >::calc_depths ()
{
  depths.resize(n_active_elements,0);
  index_t index = 0;
  for (index_t k = 0; k < nz; ++k)
    for (index_t j = 0; j < ny; ++j)
      for (index_t i = 0; i < nx; ++i, ++index)
        {
          if (sp_actnum[index])
            depths[ext_to_int[index]] = get_center_zcorn(top_cube_index(i,j,k, nx, ny));
        }
  return 0;
}

// calculating method have been taken from td eclipse (page 896)
template<class strategy_t>
typename mesh_grdecl<strategy_t>::item_t mesh_grdecl<strategy_t>::calculate_tran(const index_t index1, const index_t index2,
    const fpoint3d& A, const fpoint3d& D1,  const point_side_t &side2, const fpoint3d &center2, direction d_dir) const
  {
    item_t tran;
    fpoint3d centerSide2 = get_cube_center(side2);
    fpoint3d D2 = fpoint3d::abs(centerSide2-center2);

    float koef1, koef2 ; //koef = (A,Di)/(Di,Di)
    koef1 = A*D1 / (D1*D1);
    koef2 = A*D2 / (D2*D2);

    if (koef1 < 10e-16 || koef2 < 10e-16)
      {
        BOSWARN (section::mesh, level::warning) 
          << boost::format ("For indexes (%d, %d) transmissibility will be set to 0 because koef1 = 0 (%f) or koef2 = 0 (%f)")
          % index1 % index2 % koef1 % koef2 
          << bs_end;

        return 0;
      }

    item_t Ti, Tj;

    item_t ntg_index1 = 1;
    item_t ntg_index2 = 1;
    if (!sp_ntg.empty ())
      {
        ntg_index1 = sp_ntg[index1];
        ntg_index2 = sp_ntg[index2];
      }

    if (d_dir == along_dim1) //lengthwise OX
      {
        Ti = sp_permx[index1]*ntg_index1*koef1;
        Tj = sp_permx[index2]*ntg_index2*koef2;
        tran = darcy_constant / (1 / Ti + 1 / Tj);
        if (!sp_multx.empty ())
          tran *= sp_multx[index1];
      }
    else if (d_dir == along_dim2) //lengthwise OY
      {
        Ti = sp_permy[index1]*ntg_index1*koef1;
        Tj = sp_permy[index2]*ntg_index2*koef2;
        tran = darcy_constant / (1 / Ti + 1 / Tj);
        if (!sp_multy.empty ())
          tran *= sp_multy[index1];
      }
    else //lengthwise OZ
      {
        Ti = sp_permz[index1]*koef1;
        Tj = sp_permz[index2]*koef2;
        tran = darcy_constant / (1 / Ti + 1 / Tj);
        if (!sp_multz.empty ())
          tran *= sp_multz[index1];
      }

    return tran;
  }


template<class strategy_t>
int mesh_grdecl<strategy_t>::build_jacobian_and_flux_connections (const sp_bcsr_t &jacobian,const sp_flux_conn_iface_t & flux_conn,
                                                                  index_array_t &boundary_array)

{
  return build_jacobian_and_flux_connections_add_boundary (jacobian, flux_conn, boundary_array);

//  n_connections = 0;
//
//  index_t max_size = n_active_elements;
//
//  jacobian->get_cols_ind().clear();
//  jacobian->get_rows_ptr().clear();
//  jacobian->init_struct(max_size,max_size, max_size);
//  index_array_t* rows_ptr = &jacobian->get_rows_ptr();
//  sp_bcsr_t conn_trans;
//
//  (*rows_ptr)[0] = 0;
//  std::vector<bool> is_butting(nx*ny,false);
//  int nButting = 0;
//  boundary_array.clear();
//
//
//#pragma region //first step - define and fill - rows_ptr (jacobian)
//  for (index_t i = 0; i < nx; ++i)
//    {
//      for (index_t j = 0; j < ny; ++j)
//        {
//          index_t bottom_Kx = 0;
//          index_t bottom_Ky = 0;
//
//          index_t top_Kx, top_Ky, top_Kz;
//          for (index_t k = 0; k < nz; ++k)
//            {
//              if (!sp_actnum[i+j*nx+k*nx*ny])//skip non-active cells
//                continue;
//
//              index_array_t cube_IJK_index = top_cube_index(i,j,k, nx, ny);
//#pragma region skip butting cells && define is_ butting
//              if (k == 0)
//                {
//                  //define is_butting ([i,j,0]&&[i,j+1,0]) ([i,j,0]&&[i+1,j,0]
//                  bool flag = true;
//
//                  if (i+1 < nx)
//                    {
//                      const index_array_t &cube_IIJK_index = top_cube_index(i+1,j,k, nx, ny);
//
//                      const side_t &rightSide = get_side(right,cube_IJK_index,0);
//                      const side_t &leftSide = get_side(right,cube_IIJK_index,1);
//                      //check for butting
//
//                      for (size_t ii = 0; ii < rightSide.size(); ii++)
//                        flag = flag&&(zcorn_array[rightSide[ii]] == zcorn_array[leftSide[ii]]);
//                    }
//                  if (j+1 < ny && flag)
//                    {
//                      const index_array_t &cube_IJJK_index = top_cube_index(i,j+1,k, nx, ny);
//
//                      const side_t &bottomSide = get_side(bottom,cube_IJK_index,0);
//                      const side_t &topSide = get_side( bottom,cube_IJJK_index,1);
//
//                      //check for butting
//                      for (size_t ii = 0; ii < bottomSide.size(); ii++)
//                        flag = flag&&(zcorn_array[bottomSide[ii]] == zcorn_array[topSide[ii]]);
//                    }
//                  if (flag) nButting++;
//                  is_butting[j*nx+i] = flag;
//                }
//
//              //look only 3 positive-direction side
//              if (is_butting[j*nx+i])
//                {
//                  //speed-up by calculating just one square
//                  if ((i+1 < nx) && sp_actnum[(i+1)+j*nx + k*nx*ny])
//                    {
//                      change_row(i,j,k,*rows_ptr);
//                      change_row(i+1, j, k,*rows_ptr);
//                    }
//                  if ((j+1 < ny) && sp_actnum[i+(j+1)*nx + k*nx*ny])
//                    {
//                      change_row(i,j,k,*rows_ptr);
//                      change_row(i, j+1, k,*rows_ptr);
//                    }
//                  if ((k+1 < nz) && sp_actnum[i+j*nx + (k+1)*nx*ny])
//                    {
//                      change_row(i,j,k,*rows_ptr);
//                      change_row(i, j, k+1,*rows_ptr);
//                    }
//                  continue;
//                }
//#pragma endregion
//
//              //calculating for non butting cell!
//              top_Kx = bottom_Kx-1;
//#pragma region rightSide
//              if (i+1 < nx)
//                {
//                  //looking for incident block at right side and their zmin(kmin) & zmax(kmax)
//                  const side_t &right_side = get_side( right, cube_IJK_index,0);
//                  item_t maxZ_RS = zcorn_array[*std::max_element(right_side.begin(),right_side.end())];
//
//                  //go to top
//                  //looking for top block for right side - minZ_LS <= maxZ_RS
//                  while (top_Kx < nz-1)
//                    {
//                      top_Kx++;
//
//                      if (!sp_actnum[(i+1)+j*nx + top_Kx*nx*ny])//skip non-active cells
//                        continue;
//
//                      const side_t &temp_left_side = get_side( left, top_cube_index(i+1,j,top_Kx, nx, ny),0);
//                      item_t minZ_LS = zcorn_array[*std::min_element(temp_left_side.begin(),temp_left_side.end())];
//
//                      if (minZ_LS >= maxZ_RS) //we can't go more
//                        break;
//                      else if (is_2side_crossing(right_side, temp_left_side, true))// => [i,j,k] && [i+1,j,top_Kx] have common contact
//                        {
//                          change_row(i,j,k,*rows_ptr);
//                          change_row( i+1, j, top_Kx,*rows_ptr);
//                        }
//                    }
//                }
//
//#pragma endregion
//
//              top_Ky = bottom_Ky-1;
//#pragma region bottomSide
//              if (j+1 < ny)
//                {
//                  const side_t &bottom_side = get_side( bottom,cube_IJK_index,0);
//                  item_t maxZ_BS = zcorn_array[*std::max_element(bottom_side.begin(),bottom_side.end())];
//
//                  //go to top
//                  //looking for top block for bottom side - minZ_TS <= maxT_RS
//                  while (top_Ky < nz-1)
//                    {
//                      top_Ky++;
//                      if (!sp_actnum[i+(j+1)*nx + top_Ky*nx*ny])//skip non-active cells
//                        continue;
//
//                      const side_t &temp_top_side = get_side( left, top_cube_index(i,j+1,top_Ky, nx, ny),0);
//                      item_t minZ_TS = zcorn_array[*std::min_element(temp_top_side.begin(),temp_top_side.end())];
//                      if (minZ_TS >= maxZ_BS) //we can't go more
//                        break;
//                      else if (is_2side_crossing(bottom_side, temp_top_side,!true))// => [i,j,k] && [i,j+1,top_Ky] have common contact
//                        {
//                          change_row(i,j,k,*rows_ptr);
//                          change_row( i, j+1, top_Ky,*rows_ptr);
//                        }
//                    }
//                }
//
//#pragma endregion
//
//#pragma region upperSide
//              if (k+1 < nz)
//                {
//                  top_Kz = k+1;
//                  while (top_Kz < nz)
//                    {
//                      //TODO add null volume case
//                      if (!sp_actnum[i+j*nx+top_Kz*nx*ny])
//                        {
//                          top_Kz++;
//                          continue;
//                        }
//                      break;
//                    }
//                  if (top_Kz == nz) //no connections
//                    continue;
//                  else if (sp_actnum[i+j*nx+top_Kz*nx*ny])
//                    {
//                      // => [i,j,k] && [i,j,k_Topz] have common contact (take square upperSide) -> findAreaOfSide(upperSide[0],upperSide[1],upperSide[2],upperSide[3]);
//                      change_row(i,j,k,*rows_ptr);
//                      change_row(i, j, top_Kz,*rows_ptr);
//                    }
//                }
//#pragma endregion
//              bottom_Kx = top_Kx;
//              bottom_Ky = top_Ky;
//
//            }
//        }
//    }
//#pragma endregion
//
//  //////jacobian//////////////////////
//  //sum rows_ptr
//  for (size_t i = 1; i < rows_ptr->size(); ++i)
//    {
//      (*rows_ptr)[i]++;
//      (*rows_ptr)[i] += (*rows_ptr)[i-1];
//    }
//  //create cols_ind
//  index_array_t* cols_ind = &jacobian->get_cols_ind();
//
//  cols_ind->resize((*rows_ptr)[rows_ptr->size()-1],-1);
//  ////////transmis/////////////////////////
//  index_t cols_ind_n = (index_t)cols_ind->size();
//
//
//  index_t con_num = (cols_ind_n-max_size)/2;//connection number
//  n_connections = con_num;
//  
//  conn_trans = flux_conn->get_conn_trans();
//  conn_trans->init_struct(con_num, 2*con_num, 2*con_num);
//
//  index_array_t *rows_ptr_transmis = &conn_trans->get_rows_ptr();
//  index_array_t *cols_ind_transmis = &conn_trans->get_cols_ind();
//  rhs_item_array_t *values_transmis = &conn_trans->get_values();
//  index_array_t &matrix_block_idx_plus = flux_conn->get_matrix_block_idx_plus();
//  index_array_t &matrix_block_idx_minus = flux_conn->get_matrix_block_idx_minus();
//
//  matrix_block_idx_minus.resize(con_num*2,-1);
//  matrix_block_idx_plus.resize(con_num*2,-1);
//
//  if (!con_num)
//    {
//      for (int i = 0; i < cols_ind_n; ++i)
//        (*cols_ind)[i] = i;
//      return 0;
//    }
//  for (index_t i = 0; i < con_num+1; ++i)
//    (*rows_ptr_transmis)[i] = i*2;
//
//  //additional array for current index in rows_ptr
//  index_array_t curIndex;
//  curIndex.clear();//curIndex.push_back(0);
//  for (size_t i = 0; i < rows_ptr->size(); ++i)
//    curIndex.push_back((*rows_ptr)[i]);
//
//  con_num = 0;
//
//#pragma region //second step - fill and define cols_ind
//  for (index_t i = 0; i < nx; ++i)
//    {
//      for (index_t j = 0; j < ny; ++j)
//        {
//          index_t bottom_Kx = 0;
//          index_t bottom_Ky = 0;
//
//          index_t top_Kx, top_Ky, top_Kz;
//          for (index_t k = 0; k < nz; ++k)
//            {
//              index_t index1_extern = i+j*nx+k*nx*ny;
//
//              if (!sp_actnum[index1_extern])//skip-non-active cells
//                continue;
//
//              index_t index1 = XYZ_to_inside(i,j,k);
//
//              //condition of first place for own (always)
//              index_t tmp = curIndex[index1];
//              curIndex[index1] = (*rows_ptr)[index1];
//              change_col_add(index1,index1,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,false,con_num,*rows_ptr);
//              if (tmp == (*rows_ptr)[index1])
//                tmp++;
//              curIndex[index1] = tmp;
//
//              const ijk_cube_t &cube_IJK = top_cube(i,j,k);
//              const index_array_t &cube_IJK_index = top_cube_index(i,j,k, nx, ny);
//
//              fpoint3d center1 = get_cube_center (cube_IJK);
//
//              //look only 3 positive-direction side
//#pragma region quick calculating
//              if (is_butting[j*nx+i])
//                {
//                  //speed-up by calculating just one square
//                  if ((i+1 < nx) && sp_actnum[(i+1)+j*nx + k*nx*ny])
//                    {
//                      const point_side_t &rightSide = get_side( right,cube_IJK,0);
//
//                      fpoint3d sqRS = get_projection_on_all_axis_for_one_side(rightSide);
//                      //change jacobian
//                      index_t index2 = XYZ_to_inside(i+1,j,k);
//
//                      change_col_add(index1, index2,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,true,con_num,*rows_ptr);
//                      if (curIndex[index2] == (*rows_ptr)[index2])
//                        curIndex[index2]++;
//                      change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, false,true,con_num,*rows_ptr);
//
//                      //change flux_connection
//                      const ijk_cube_t &cube2 = top_cube(i+1,j,k);
//                      fpoint3d D1 = fpoint3d::abs(get_cube_center(rightSide)-center1);
//                      fpoint3d center2 = get_cube_center(cube2);
//                      item_t tran = calculate_tran(index1_extern, ((i+1)+j*nx+k*nx*ny), sqRS,D1, rightSide, center2, along_dim1);
//
//                      (*cols_ind_transmis)[con_num++] = index1;
//                      (*cols_ind_transmis)[con_num++] = index2;
//                      values_transmis->push_back(tran);
//                      values_transmis->push_back(-tran);
//                    }
//                  if ((j+1 < ny) && sp_actnum[i+(j+1)*nx + k*nx*ny])//skip non-active
//                    {
//                      const point_side_t &bottomSide = get_side( bottom,cube_IJK,0);
//                      fpoint3d sqBS = get_projection_on_all_axis_for_one_side(bottomSide);
//
//                      const ijk_cube_t &cube2 = top_cube(i,j+1,k);
//                      index_t index2 = XYZ_to_inside(i,j+1,k);
//
//                      change_col_add(index1, index2, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,true,con_num,*rows_ptr);
//                      if (curIndex[index2] == (*rows_ptr)[index2])
//                        curIndex[index2]++;
//                      change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, false,true,con_num,*rows_ptr);
//
//                      //change flux_connection
//                      fpoint3d center2 = get_cube_center(cube2);
//                      fpoint3d D1 = fpoint3d::abs(get_cube_center(bottomSide)-center1);
//                      item_t tran = calculate_tran(index1_extern, (i+(j+1)*nx+k*nx*ny), sqBS, D1, bottomSide,center2,along_dim2);
//
//                      (*cols_ind_transmis)[con_num++] = index1;
//                      (*cols_ind_transmis)[con_num++] = index2;
//                      values_transmis->push_back(tran);
//                      values_transmis->push_back(-tran);
//                    }
//                  if ((k+1 < nz) && sp_actnum[i+j*nx + (k+1)*nx*ny])//skip non-active
//                    {
//                      const point_side_t &upperSide = get_side( upper,cube_IJK,0);
//
//                      //change jacobian
//                      index_t index2 = XYZ_to_inside(i,j,k+1);
//
//                      change_col_add(index1, index2,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,true,con_num,*rows_ptr);
//                      if (curIndex[index2] == (*rows_ptr)[index2])
//                        curIndex[index2]++;
//                      change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, false,true,con_num,*rows_ptr);
//
//                      fpoint3d sqUS = get_projection_on_all_axis_for_one_side(upperSide);
//                      fpoint3d D1 = fpoint3d::abs(get_cube_center(upperSide)-center1);
//
//                      const ijk_cube_t &cube2 = top_cube(i,j,k+1);
//
//                      //change flux_connection
//                      fpoint3d center2 = get_cube_center(cube2);
//                      item_t tran = calculate_tran(index1_extern, (i+j*nx+(k+1)*nx*ny), sqUS, D1, upperSide,center2,along_dim3);
//
//                      (*cols_ind_transmis)[con_num++] = index1;
//                      (*cols_ind_transmis)[con_num++] = index2;
//                      values_transmis->push_back(tran);
//                      values_transmis->push_back(-tran);
//
//                    }
//                  continue;
//                }
//#pragma endregion
//
//              top_Kx = bottom_Kx-1;
//#pragma region rightSide
//              if (i+1 < nx)
//                {
//                  //looking for incident block at right side and their zmin(kmin) & zmax(kmax)
//                  const side_t &right_side = get_side( right, cube_IJK_index,0);
//                  item_t maxZ_RS = zcorn_array[*std::max_element(right_side.begin(),right_side.end())];
//
//                  const point_side_t &rightSide = get_side( right,cube_IJK,0);
//                  fpoint3d centerSide1 = get_cube_center(rightSide);
//                  fpoint3d D1 = fpoint3d::abs(centerSide1-center1);
//
//                  //go to top
//                  //looking for top block for right side - minZ_LS <= maxZ_RS
//                  while (top_Kx < nz-1)
//                    {
//                      top_Kx++;
//                      if (!sp_actnum[(i+1)+j*nx + top_Kx*nx*ny])//skip non-active
//                        continue;
//
//                      const side_t &temp_left_side = get_side( left, top_cube_index(i+1,j,top_Kx, nx, ny),0);
//                      item_t minZ_LS = zcorn_array[*std::min_element(temp_left_side.begin(),temp_left_side.end())];
//                      if (minZ_LS >= maxZ_RS) //we can't go more
//                        break;
//                      else if (is_2side_crossing(right_side, temp_left_side, true))// => [i,j,k] && [i+1,j,top_Kx] have common contact
//                        {
//                          //looking for square of intersection
//                          const ijk_cube_t &cube2 = top_cube(i+1,j,top_Kx);
//                          const point_side_t &leftSide = get_side( left, cube2,0);
//                          fpoint3d A = get_side_crossing_projection_on_all_axis(leftSide, rightSide);
//
//                          //change jacobian
//                          index_t index2 = XYZ_to_inside(i+1,j,top_Kx);
//
//                          change_col_add(index1, index2, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,true,con_num,*rows_ptr);
//                          if (curIndex[index2] == (*rows_ptr)[index2])
//                            curIndex[index2]++;
//                          change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, false,true,con_num,*rows_ptr);
//
//                          //change flux_connection
//                          fpoint3d center2 = get_cube_center(cube2);
//                          item_t tran = calculate_tran(index1_extern, ((i+1)+j*nx+top_Kx*nx*ny), A,D1, leftSide,center2,along_dim1);
//
//                          (*cols_ind_transmis)[con_num++] = index1;//add [i,j,k]->index1 && [i+1,j,bottom_K]->index2=>index1&index2
//                          (*cols_ind_transmis)[con_num++] = index2;
//                          values_transmis->push_back(tran);
//                          values_transmis->push_back(-tran);
//                        }
//                    }
//                }
//#pragma endregion
//
//              top_Ky = bottom_Ky-1;
//#pragma region bottomSide
//              if (j+1 < ny)
//                {
//                  //looking for incident block at bottom side and their zmin(kmin) & zmax(kmax)
//                  const point_side_t &bottomSide = get_side( bottom,cube_IJK,0);
//                  fpoint3d centerSide1 = get_cube_center(bottomSide);
//                  fpoint3d D1 = fpoint3d::abs(centerSide1-center1);
//
//                  const side_t &bottom_side = get_side( bottom,cube_IJK_index,0);
//                  item_t maxZ_BS = zcorn_array[*std::max_element(bottom_side.begin(),bottom_side.end())];
//
//                  //go to top
//                  //looking for top block for bottom side - minZ_TS <= maxT_RS
//                  while (top_Ky < nz-1)
//                    {
//                      top_Ky++;
//                      if (!sp_actnum[i+(j+1)*nx + top_Ky*nx*ny])//skip non-active cells
//                        continue;
//
//                      const side_t &temp_top_side  = get_side( left, top_cube_index(i,j+1,top_Ky, nx, ny),0);
//                      item_t minZ_TS = zcorn_array[*std::min_element(temp_top_side.begin(),temp_top_side.end())];
//                      if (minZ_TS >= maxZ_BS) //we can't go more
//                        break;
//                      else if (is_2side_crossing(bottom_side, temp_top_side,!true))// => [i,j,k] && [i,j+1,top_Ky] have common contact
//                        {
//                          const ijk_cube_t &cube2 = top_cube(i,j+1,top_Ky);
//                          const point_side_t &topSide = get_side( top, cube2,0);
//                          fpoint3d A = get_side_crossing_projection_on_all_axis(topSide, bottomSide);
//
//                          //change jacobian
//                          index_t index2 = XYZ_to_inside(i,j+1,top_Ky);
//
//                          change_col_add(index1, index2,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,true,con_num,*rows_ptr);
//                          if (curIndex[index2] == (*rows_ptr)[index2])
//                            curIndex[index2]++;
//                          change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, false,true,con_num,*rows_ptr);
//
//                          //change flux_connection
//                          fpoint3d center2 = get_cube_center(cube2);
//                          item_t tran = calculate_tran(index1_extern, (i+(j+1)*nx+top_Ky*nx*ny), A, D1, topSide,center2,along_dim2);
//
//                          (*cols_ind_transmis)[con_num++] = index1;//add [i,j,k]->index1 && [i+1,j,bottom_K]->index2=>index1&index2
//                          (*cols_ind_transmis)[con_num++] = index2;
//                          values_transmis->push_back(tran);
//                          values_transmis->push_back(-tran);
//
//                        }
//                    }
//                }
//#pragma endregion
//
//#pragma region 	upperSide
//              if (k+1 < nz)
//                {
//                  top_Kz = k+1;
//                  while (top_Kz < nz)
//                    {
//                      if (!sp_actnum[i+j*nx+top_Kz*nx*ny])
//                        {
//                          top_Kz++;
//                          continue;
//                        }
//                      break;
//                    }
//                  if (top_Kz == nz) //no connection
//                    continue;
//                  else if (sp_actnum[i+j*nx+top_Kz*nx*ny])
//                    {
//                      //=> [i,j,k] && [i,j,top_K] have common contact  -> 	findAreaOfSide(upperSide[0],upperSide[1],upperSide[2],upperSide[3]);
//                      const point_side_t &upperSide = get_side( upper,cube_IJK,0);
//                      fpoint3d D1 = fpoint3d::abs(get_cube_center(upperSide)-center1);
//
//                      //change jacobian
//                      index_t index2 = XYZ_to_inside(i,j,top_Kz);
//
//                      change_col_add(index1, index2,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,true,con_num,*rows_ptr);
//                      if (curIndex[index2] == (*rows_ptr)[index2])
//                        curIndex[index2]++;
//                      change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, false,true,con_num,*rows_ptr);
//
//                      const ijk_cube_t &cube2 = top_cube(i,j,top_Kz);
//                      fpoint3d A = get_projection_on_all_axis_for_one_side(upperSide);//because it's too expensive calculate as triangles...
//
//                      //change flux_connection
//                      fpoint3d center2 = get_cube_center(cube2);
//                      item_t tran = calculate_tran(index1_extern, (i+j*nx+top_Kz*nx*ny), A, D1, upperSide,center2,along_dim3);
//
//                      (*cols_ind_transmis)[con_num++] = index1;//add [i,j,k]->index1 && [i+1,j,bottom_K]->index2=>index1&index2
//                      (*cols_ind_transmis)[con_num++] = index2;
//                      values_transmis->push_back(tran);
//                      values_transmis->push_back(-tran);
//                    }
//                }
//#pragma endregion
//              bottom_Kx = top_Kx;
//              bottom_Ky = top_Ky;
//            }
//        }
//    }
//#pragma endregion
//  return 0;
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::init_props(const sp_idata_t &idata)
{
  base_t::init_props (idata);
  
  // init ZCORN
  zcorn_array = idata->get_float_non_empty_array("ZCORN");
  min_z = *(std::min_element(zcorn_array.begin(),zcorn_array.end()));
  max_z = *(std::max_element(zcorn_array.begin(),zcorn_array.end()));
  
  // init COORD
  coord_array = idata->get_float_non_empty_array("COORD");
  
  max_x = min_x = coord_array[0];
  max_y = min_y = coord_array[1];

  array_float16_t::iterator it, e = coord_array.end ();
  
  for (it = coord_array.begin (); it != e; it += 6)
    {
      // move matching points apart
      if (it[2] == it[5])
        it[5] += 1.0f;
        
      //looking for max&min coordinates
      if (min_x > it[0]) min_x = it[0];
      if (min_x > it[3]) min_x = it[3];
      
      if (min_y > it[1]) min_y = it[1];
      if (min_y > it[4]) min_y = it[4];
      
      if (max_x < it[0]) max_x = it[0];
      if (max_x < it[3]) max_x = it[3];

      if (max_y < it[1]) max_y = it[1];
      if (max_y < it[4]) max_y = it[4];
      
    }
  
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::check_data() const
{
  base_t::check_data ();
  
  if (min_x < 0)
    bs_throw_exception (boost::format ("min_x = %d is out of range")% min_x);
  if (min_y < 0)
    bs_throw_exception (boost::format ("min_y = %d is out of range")% min_y);
  if (min_z < 0)
    bs_throw_exception (boost::format ("min_z = %d is out of range")% min_z);
    
  if (!coord_array.size ())
    bs_throw_exception ("COORD array is not initialized");
  if (!zcorn_array.size ())
    bs_throw_exception ("ZCORN array is not initialized");
}

template<class strategy_t>
void mesh_grdecl<strategy_t>:: get_block_dx_dy_dz(index_t n_elem, item_t &dx, item_t &dy, item_t &dz) const
  {
    index_t i,j,k;

    inside_to_XYZ(n_elem,i,j,k);
    const ijk_cube_t &cube_IJK = top_cube(i,j,k);

    const point_side_t &tmp_left_side = get_side(left,cube_IJK,false);
    const point_side_t &tmp_right_side = get_side(left,cube_IJK,true);

    const point_side_t &tmp_top_side = get_side(top,cube_IJK,false);
    const point_side_t &tmp_bottom_side = get_side(top,cube_IJK,true);

    const point_side_t &tmp_upper_side = get_side(upper,cube_IJK,false);
    const point_side_t &tmp_lower_side = get_side(upper,cube_IJK,true);

    dx = dy = dz = 0.0;
    for (index_t ii = 0; ii < 4; ++ii)
      {
        dx += tmp_left_side[ii].x - tmp_right_side[ii].x;
        dy += tmp_top_side[ii].y - tmp_bottom_side[ii].y;
        dz += tmp_upper_side[ii].z - tmp_lower_side[ii].z;
      }

    dx = fabs(dx/4);
    dy = fabs(dy/4);
    dz = fabs(dz/4);
  }

template<class strategy_t>
float mesh_grdecl<strategy_t>:: get_block_dx(index_t n_elem) const
  {
    index_t i,j,k;
    inside_to_XYZ(n_elem,i,j,k);
    const ijk_cube_t &cube_IJK = top_cube(i,j,k);

    const point_side_t &tmp_left_side = get_side(left,cube_IJK,false);
    const point_side_t &tmp_right_side = get_side(left,cube_IJK,true);

    double dx = 0.0;
    for (index_t ii = 0; ii < 4; ++ii)
      dx += tmp_left_side[ii].x - tmp_right_side[ii].x;
    return fabs(dx/4);
  }

template<class strategy_t>
float mesh_grdecl<strategy_t>:: get_block_dy(index_t n_elem) const
  {
    index_t i,j,k;
    inside_to_XYZ(n_elem,i,j,k);
    const ijk_cube_t &cube_IJK = top_cube(i,j,k);

    const point_side_t &tmp_top_side = get_side(top,cube_IJK,false);
    const point_side_t &tmp_bottom_side = get_side(top,cube_IJK,true);

    double dy = 0.0;
    for (index_t ii = 0; ii < 4; ++ii)
      dy += tmp_top_side[ii].y - tmp_bottom_side[ii].y;
    return fabs(dy/4);
  }
  
template<class strategy_t>
float mesh_grdecl<strategy_t>:: get_block_dz(index_t n_elem) const
  {
    index_t i,j,k;
    inside_to_XYZ(n_elem,i,j,k);
    const ijk_cube_t &cube_IJK = top_cube (i,j,k);

    const point_side_t &tmp_upper_side = get_side (upper,cube_IJK,false);
    const point_side_t &tmp_lower_side = get_side (upper,cube_IJK,true);

    double dz = 0.0;
    for (index_t ii = 0; ii < 4; ++ii)
      dz += tmp_upper_side[ii].z - tmp_lower_side[ii].z;

    return fabs(dz/4);
  }

template<class strategy_t>
float  mesh_grdecl<strategy_t>:: get_dtop(index_t n_elem) const
{
  BS_ASSERT (n_elem >= 0) (n_elem);

  index_t i,j,k;
  inside_to_XYZ(n_elem,i,j,k);
  const ijk_cube_t &cube_IJK = top_cube(i,j,k);

  const point_side_t &tmp_upper_side = get_side(upper,cube_IJK,false);
  const point_side_t &tmp_lower_side = get_side(upper,cube_IJK,true);

  double dz = 0.0;
  for (index_t ii = 0; ii < 4; ++ii)
    dz += tmp_upper_side[ii].z - tmp_lower_side[ii].z;
  dz = fabs(dz/8);

  return get_cube_center(cube_IJK).z - dz;
}

template<class strategy_t>
bool mesh_grdecl<strategy_t>::is_2side_crossing (const side_t &side1, const side_t &side2, bool x_dir)
{
  //size_t j = 0;
  if (x_dir)
    {
      for (size_t i = 0, j = 0; i < side1.size(); ++i, j = (j+1)%2)
        {
          if ((zcorn_array[side2[i]] > zcorn_array[side1[j]]) && (zcorn_array[side2[i]] < zcorn_array[side1[j+2]]))
            return true;
        }
    }
  else
    {
      for (size_t i = 0, j = 0; i < side1.size(); ++i, j = (j+1)%2)
        {
          if ((zcorn_array[side2[i]] < zcorn_array[side1[j]]) && (zcorn_array[side2[i]] > zcorn_array[side1[j+2]]))
            return true;
        }
    }
  return false;
}

template<class strategy_t>
void mesh_grdecl<strategy_t>::generate_array()
{
#if 0
  index_t n_size = n_elements;
  sp_poro->clear();
  sp_ntg->clear();

  sp_permx->clear();
  sp_permy->clear();
  sp_permz->clear();

  sp_multx->clear();
  sp_multy->clear();
  sp_multz->clear();


  for (index_t i =0; i < n_size; ++i)
    {
      sp_poro->push_back(0.2f);
      sp_ntg->push_back(0.4f);

      sp_permx->push_back(10.0f);
      sp_permy->push_back(10.0f);
      sp_permz->push_back(10.0f);

      sp_multx->push_back(1.0f);
      sp_multy->push_back(1.0f);
      sp_multz->push_back(1.0f);
    }
#endif 
}


template<class strategy_t>
int mesh_grdecl<strategy_t>::find_neighbours(sp_bcsr_t &jacobian)
{
  index_t max_size = n_active_elements;
  jacobian->init(max_size,max_size, 1, max_size);//1 - n_block_size
  index_array_t* rows_ptr = &jacobian->get_rows_ptr();
  (*rows_ptr)[0] = 0;

  std::vector<bool> is_butting(nx*ny,false);
  int nButting = 0;

#pragma region //first step - define and fill - rows_ptr (jacobian)
  for (index_t i = 0; i < nx; i++)
    {
      for (index_t j = 0; j < ny; j++)
        {
          index_t bottom_Kx = 0;
          index_t bottom_Ky = 0;

          index_t top_Kx, top_Ky, top_Kz;
          for (index_t k = 0; k < nz; k++)
            {
              index_t index1_extern = i+j*nx+k*nx*ny;
              if (!sp_actnum[index1_extern])//skip non-active cells
                continue;

              const cube_index_t &cube_IJK_index = top_cube_index(i,j,k, nx, ny);
              index_t index2_extern;
#pragma region skip butting cells && define is Butting
              if (k == 0)
                {
                  //define is_butting ([i,j,0]&&[i,j+1,0])
                  //				   ([i,j,0]&&[i+1,j,0]
                  bool flag = true;

                  if (i+1 < nx)
                    {
                      const cube_index_t &cube_IIJK_index = top_cube_index(i+1,j,k, nx, ny);

                      const side_t &rightSide = get_side( right,cube_IJK_index,0);
                      const side_t &leftSide = get_side( right,cube_IIJK_index,1);
                      //check for butting

                      for (size_t ii = 0; ii < rightSide.size(); ii++)
                        flag = flag&&(zcorn_array[rightSide[ii]] == zcorn_array[leftSide[ii]]);
                    }
                  if (j+1 < ny && flag)
                    {
                      const cube_index_t &cube_IJJK_index = top_cube_index(i,j+1,k, nx, ny);

                      const side_t &bottomSide = get_side( bottom,cube_IJK_index,0);
                      const side_t &topSide = get_side( bottom,cube_IJJK_index,1);
                      //check for butting
                      for (size_t ii = 0; ii < bottomSide.size(); ii++)
                        flag = flag&&(zcorn_array[bottomSide[ii]] == zcorn_array[topSide[ii]]);
                    }
                  if (flag) nButting++;
                  is_butting[j*nx+i] = flag;
                }

              //look only 3 positive-direction side
              if (is_butting[j*nx+i])
                {
                  //speed-up by calculating just one square
                  index2_extern = (i+1)+j*nx + k*nx*ny;
                  if (i+1 < nx)
                    {
                      if (sp_actnum[index2_extern])
                        {
                          ((*rows_ptr)[index1_extern+1])++;
                          ((*rows_ptr)[index2_extern+1])++;
                        }
                    }
                  index2_extern = i+(j+1)*nx + k*nx*ny;
                  if (j+1 < ny)
                    {
                      if  (sp_actnum[index2_extern])
                        {
                          ((*rows_ptr)[index1_extern+1])++;
                          ((*rows_ptr)[index2_extern+1])++;
                        }
                    }
                  index2_extern = i+j*nx + (k+1)*nx*ny;
                  if (k+1 < nz)
                    {
                      if (sp_actnum[index2_extern])
                        {
                          ((*rows_ptr)[index1_extern+1])++;
                          ((*rows_ptr)[index2_extern+1])++;
                        }
                    }
                  continue;
                }
#pragma endregion

              //calculating for non butting cell!
              top_Kx = bottom_Kx-1;
#pragma region rightSide
              if (i+1 < nx)
                {
                  //looking for incident block at right side and their zmin(kmin) & zmax(kmax)
                  const side_t &right_side = get_side( right, cube_IJK_index,0);
                  item_t maxZ_RS = zcorn_array[*std::max_element(right_side.begin(),right_side.end())];

                  //go to top
                  //looking for top block for right side - minZ_LS <= maxZ_RS
                  while (top_Kx < nz-1)
                    {
                      top_Kx++;
                      index2_extern = (i+1)+j*nx + top_Kx*nx*ny;
                      if (!sp_actnum[index2_extern])//skip non-active cells
                        continue;

                      const side_t &temp_left_side = get_side( left, top_cube_index(i+1,j,top_Kx, nx, ny),0);
                      item_t minZ_LS = zcorn_array[*std::min_element(temp_left_side.begin(),temp_left_side.end())];

                      if (minZ_LS >= maxZ_RS) //we can't go more
                        break;
                      else if (is_2side_crossing(right_side, temp_left_side, true))// => [i,j,k] && [i+1,j,top_Kx] have common contact
                        {
                          ((*rows_ptr)[index1_extern+1])++;
                          ((*rows_ptr)[index2_extern+1])++;
                        }
                    }
                }

#pragma endregion

              top_Ky = bottom_Ky-1;
#pragma region bottomSide
              if (j+1 < ny)
                {
                  const side_t &bottom_side = get_side( bottom,cube_IJK_index,0);
                  item_t maxZ_BS = zcorn_array[*std::max_element(bottom_side.begin(),bottom_side.end())];

                  //go to top
                  //looking for top block for bottom side - minZ_TS <= maxT_RS
                  while (top_Ky < nz-1)
                    {
                      top_Ky++;
                      index2_extern = i+(j+1)*nx + top_Ky*nx*ny;

                      if (!sp_actnum[index2_extern])//skip non-active cells
                        continue;

                      const side_t &temp_top_side = get_side( left, top_cube_index(i,j+1,top_Ky, nx, ny),0);
                      item_t minZ_TS = zcorn_array[*std::min_element(temp_top_side.begin(),temp_top_side.end())];
                      if (minZ_TS >= maxZ_BS) //we can't go more
                        break;
                      else if (is_2side_crossing(bottom_side, temp_top_side,!true))// => [i,j,k] && [i,j+1,top_Ky] have common contact
                        {
                          ((*rows_ptr)[index1_extern+1])++;
                          ((*rows_ptr)[index2_extern+1])++;
                        }
                    }
                }

#pragma endregion

#pragma region upperSide
              if (k+1 < nz)
                {
                  top_Kz = k+1;
                  index2_extern = i+j*nx+top_Kz*nx*ny;

                  while (top_Kz < nz)
                    {
                      if (!sp_actnum[index2_extern])
                        {
                          top_Kz++;
                          index2_extern = i+j*nx+top_Kz*nx*ny;
                          continue;
                        }
                      break;
                    }
                  if (top_Kz == nz) //no connections
                    continue;
                  else if (sp_actnum[index2_extern])
                    {
                      // => [i,j,k] && [i,j,k_Topz] have common contact (take square upperSide) -> findAreaOfSide(upperSide[0],upperSide[1],upperSide[2],upperSide[3]);
                      ((*rows_ptr)[index1_extern+1])++;
                      ((*rows_ptr)[index2_extern+1])++;
                    }
                }
#pragma endregion
              bottom_Kx = top_Kx;
              bottom_Ky = top_Ky;
            }
        }
    }
#pragma endregion

  //////jacobian//////////////////////
  //sum rows_ptr
  for (size_t i = 1; i < rows_ptr->size(); i++)
    (*rows_ptr)[i] += (*rows_ptr)[i-1];

  //create cols_ind
  index_array_t* cols_ind = &jacobian->get_cols_ind();

  cols_ind->resize((*rows_ptr)[rows_ptr->size()-1],-1);
  ////////transmis/////////////////////////
  //index_t cols_ind_n = (index_t)cols_ind->size();


  //additional array for current index in rows_ptr
  index_array_t curIndex;
  curIndex.clear();
  for (size_t i = 0; i < rows_ptr->size(); i++)
    curIndex.push_back((*rows_ptr)[i]);

  //we can optimize, if we using our structure, and just define square
#pragma region //second step - fill and define cols_ind
  for (index_t i = 0; i < nx; i++)
    {
      for (index_t j = 0; j < ny; j++)
        {
          index_t bottom_Kx = 0;
          index_t bottom_Ky = 0;

          index_t top_Kx, top_Ky, top_Kz;
          for (index_t k = 0; k < nz; k++)
            {
              index_t index1_extern = i+j*nx+k*nx*ny;

              if (!sp_actnum[index1_extern])//skip-non-active cells
                continue;

              const ijk_cube_t &cube_IJK = top_cube(i,j,k);
              const cube_index_t &cube_IJK_index = top_cube_index(i,j,k, nx, ny);

              fpoint3d center1 = get_cube_center(cube_IJK);
              index_t index2_extern;

              //look only 3 positive-direction side
#pragma region quick calculating
              if (is_butting[j*nx+i])
                {
                  //speed-up by calculating just one square
                  index2_extern = (i+1)+j*nx + k*nx*ny;
                  if (i+1 < nx)
                    {
                      if (sp_actnum[index2_extern])
                        {
                          //change jacobian
                          change_col(index1_extern, index2_extern, *cols_ind,curIndex);
                          change_col(index2_extern, index1_extern, *cols_ind,curIndex);
                        }
                    }
                  index2_extern = i+(j+1)*nx + k*nx*ny;
                  if (j+1 < ny)
                    {
                      if (sp_actnum[index2_extern])
                        {
                          change_col(index1_extern, index2_extern, *cols_ind,curIndex);
                          change_col(index2_extern, index1_extern, *cols_ind,curIndex);
                        }
                    }
                  index2_extern = i+j*nx + (k+1)*nx*ny;
                  if (k+1 < nz)
                    {
                      if (sp_actnum[index2_extern])
                        {
                          //change jacobian
                          change_col(index1_extern, index2_extern,*cols_ind,curIndex);
                          change_col(index2_extern, index1_extern, *cols_ind,curIndex);
                        }
                    }
                  continue;
                }
#pragma endregion

              top_Kx = bottom_Kx-1;
#pragma region rightSide
              if (i+1 < nx)
                {
                  //looking for incident block at right side and their zmin(kmin) & zmax(kmax)
                  const side_t &right_side = get_side( right, cube_IJK_index,0);
                  item_t maxZ_RS = zcorn_array[*std::max_element(right_side.begin(),right_side.end())];

                  //go to top
                  //looking for top block for right side - minZ_LS <= maxZ_RS
                  while (top_Kx < nz-1)
                    {
                      top_Kx++;
                      index2_extern = (i+1)+j*nx + top_Kx*nx*ny;

                      if (!sp_actnum[index2_extern])//skip non-active
                        continue;

                      const side_t &temp_left_side = get_side( left, top_cube_index(i+1,j,top_Kx, nx, ny),0);
                      item_t minZ_LS = zcorn_array[*std::min_element(temp_left_side.begin(),temp_left_side.end())];
                      if (minZ_LS >= maxZ_RS) //we can't go more
                        break;
                      else if (is_2side_crossing(right_side, temp_left_side, true))// => [i,j,k] && [i+1,j,top_Kx] have common contact
                        {
                          //change jacobian
                          change_col(index1_extern, index2_extern, *cols_ind,curIndex);
                          change_col(index2_extern, index1_extern, *cols_ind,curIndex);
                        }
                    }
                }
#pragma endregion

              top_Ky = bottom_Ky-1;
#pragma region bottomSide
              if (j+1 < ny)
                {
                  //looking for incident block at bottom side and their zmin(kmin) & zmax(kmax)
                  const side_t &bottom_side = get_side( bottom,cube_IJK_index,0);
                  item_t maxZ_BS = zcorn_array[*std::max_element(bottom_side.begin(),bottom_side.end())];

                  //go to top
                  //looking for top block for bottom side - minZ_TS <= maxT_RS
                  while (top_Ky < nz-1)
                    {
                      top_Ky++;
                      index2_extern = i+(j+1)*nx + top_Ky*nx*ny;
                      if (!sp_actnum[index2_extern])//skip non-active cells
                        continue;

                      const side_t &temp_top_side  = get_side( left, top_cube_index(i,j+1,top_Ky, nx, ny),0);
                      item_t minZ_TS = zcorn_array[*std::min_element(temp_top_side.begin(),temp_top_side.end())];
                      if (minZ_TS >= maxZ_BS) //we can't go more
                        break;
                      else if (is_2side_crossing(bottom_side, temp_top_side,!true))// => [i,j,k] && [i,j+1,top_Ky] have common contact
                        {
                          //change jacobian
                          change_col(index1_extern, index2_extern,*cols_ind,curIndex);
                          change_col(index2_extern, index1_extern, *cols_ind,curIndex);
                        }
                    }
                }
#pragma endregion

#pragma region 	upperSide
              if (k+1 < nz)
                {
                  top_Kz = k+1;
                  index2_extern = i+j*nx + top_Kz*nx*ny;

                  while (top_Kz < nz)
                    {
                      if (!sp_actnum[index2_extern])
                        {
                          top_Kz++;
                          index2_extern = i+j*nx + top_Kz*nx*ny;
                          continue;
                        }
                      break;
                    }
                  if (top_Kz == nz) //no connection
                    continue;
                  else if (sp_actnum[index2_extern])
                    {
                      //=> [i,j,k] && [i,j,top_K] have common contact  -> 	findAreaOfSide(upperSide[0],upperSide[1],upperSide[2],upperSide[3]);
                      //change jacobian
                      change_col(index1_extern, index2_extern,*cols_ind,curIndex);
                      change_col(index2_extern, index1_extern, *cols_ind,curIndex);
                    }
                }
#pragma endregion
              bottom_Kx = top_Kx;
              bottom_Ky = top_Ky;
            }
        }
    }
#pragma endregion
  return 0;
}

template <typename strategy_t, typename loop_t>
struct build_jacobian_rows_class
{
  typedef typename strategy_t::index_t          index_t;
  typedef typename strategy_t::item_t           item_t;
  typedef typename strategy_t::index_array_t    index_array_t;
  typedef typename strategy_t::item_array_t     item_array_t;

  typedef mesh_grdecl <strategy_t>              mesh_t;
  typedef typename mesh_t::side_t               side_t;
  typedef typename mesh_t::cube_index_t         cube_index_t;

  build_jacobian_rows_class (mesh_grdecl <strategy_t> *mesh, loop_t *loop, std::set <index_t, std::less <index_t> > &boundary_set, index_array_t *rows_ptr)
  : mesh (mesh)
  , loop (loop)
  , boundary_set (boundary_set)
  , rows_ptr (rows_ptr)
  , nx (mesh->nx)
  , ny (mesh->ny)
  , nz (mesh->nz)
  {
  }

  void
  prepare (index_t, index_t, index_t, const cube_index_t *cube_IJK_index_, bool)
  {
    cube_IJK_index = cube_IJK_index_;
  }

  //skip butting cells && define is Butting
  void
  define_butting (index_t i, index_t j, index_t k)
  {
    if (k == 0)
      {
        //define is_butting ([i,j,0]&&[i,j+1,0])
        //				   ([i,j,0]&&[i+1,j,0]
        bool flag = true;

        if (i+1 < nx)
          {
            const cube_index_t &cube_IIJK_index   = top_cube_index(i+1,j,k, nx, ny);
            const side_t &rightSide               = get_side (right, *cube_IJK_index, 0);
            const side_t &leftSide                = get_side (right, cube_IIJK_index, 1);

            //check for butting
            for (size_t ii = 0; ii < rightSide.size(); ii++)
              flag = flag && (mesh->zcorn_array[rightSide[ii]] == mesh->zcorn_array[leftSide[ii]]);
          }
        if (j+1 < ny && flag)
          {
            const cube_index_t &cube_IJJK_index   = top_cube_index(i,j+1,k, nx, ny);
            const side_t &bottomSide              = get_side (bottom, *cube_IJK_index, 0);
            const side_t &topSide                 = get_side (bottom, cube_IJJK_index, 1);

            //check for butting
            for (size_t ii = 0; ii < bottomSide.size(); ii++)
              flag = flag && (mesh->zcorn_array[bottomSide[ii]] == mesh->zcorn_array[topSide[ii]]);
          }
        //if (flag) 
        //  nButting++;

        loop->is_butting[j*nx+i] = flag;
      }
  }

  void
  change_by_x (index_t i, index_t j, index_t km, index_t kp, bool)
  {
    (*rows_ptr)[mesh->XYZ_to_inside (i,     j, km) + 1]++;
    (*rows_ptr)[mesh->XYZ_to_inside (i + 1, j, kp) + 1]++;
  }

  void
  change_by_y (index_t i, index_t j, index_t km, index_t kp, bool)
  {
    (*rows_ptr)[mesh->XYZ_to_inside (i, j,     km) + 1]++;
    (*rows_ptr)[mesh->XYZ_to_inside (i, j + 1, kp) + 1]++;
  }

  void
  change_by_z (index_t i, index_t j, index_t km, index_t kp, bool)
  {
    (*rows_ptr)[mesh->XYZ_to_inside (i, j, km) + 1]++;
    (*rows_ptr)[mesh->XYZ_to_inside (i, j, kp) + 1]++;
  }

  void
  add_boundary (index_t external_cell_index)
  {
    boundary_set.insert (external_cell_index);
  }

  mesh_grdecl <strategy_t>                  *mesh;
  loop_t                                    *loop;
  std::set <index_t, std::less <index_t> >  &boundary_set;
  index_array_t                             *rows_ptr;
  index_t                                   nx;
  index_t                                   ny;
  index_t                                   nz;
  const cube_index_t                        *cube_IJK_index;
};

template <typename T, typename L, typename BS, typename RP>
build_jacobian_rows_class <T, L>
build_jacobian_rows (mesh_grdecl <T> *mesh, L *l, BS &bs, RP *rp)
{
  return build_jacobian_rows_class <T, L> (mesh, l, bs, rp);
}

template <typename strategy_t, typename loop_t>
struct build_jacobian_cols_class
{
  typedef typename strategy_t::index_t          index_t;
  typedef typename strategy_t::item_t           item_t;
  typedef typename strategy_t::index_array_t    index_array_t;
  typedef typename strategy_t::item_array_t     item_array_t;
  typedef typename strategy_t::rhs_item_array_t rhs_item_array_t;

  typedef mesh_grdecl <strategy_t>              mesh_t;
  typedef typename mesh_t::ijk_cube_t           ijk_cube_t;
  typedef typename mesh_t::point_side_t         point_side_t;
  typedef typename mesh_t::cube_index_t         cube_index_t;

  build_jacobian_cols_class (mesh_t *mesh, loop_t *loop, index_array_t *rows_ptr, index_array_t *cols_ind,
    index_array_t &cols_ind_transmis, rhs_item_array_t &values_transmis,
    index_array_t &matrix_block_idx_minus, index_array_t &matrix_block_idx_plus)
  : mesh (mesh)
  , loop (loop)
  , rows_ptr (rows_ptr)
  , cols_ind (cols_ind)
  , cols_ind_transmis (cols_ind_transmis)
  , values_transmis (values_transmis)
  , matrix_block_idx_minus (matrix_block_idx_minus)
  , matrix_block_idx_plus (matrix_block_idx_plus)
  , nx (mesh->nx)
  , ny (mesh->ny)
  , nz (mesh->nz)
  {
    curIndex.assign (rows_ptr->begin (), rows_ptr->end ());
  }

  void
  prepare (index_t i, index_t j, index_t k, const cube_index_t *cube_IJK_index_, bool is_active_cell)
  {
    if (is_active_cell)
      {
        cube_IJK_index  = cube_IJK_index_;
        index1          = mesh->XYZ_to_inside(i,j,k);
        index_ijk       = i + j * nx + k * nx * ny;

        //condition of first place for own (always)
        index_t tmp = curIndex[index1];
        curIndex[index1] = (*rows_ptr)[index1];
        change_col_add(index1,index1,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true, false,0,*rows_ptr);
        if (tmp == (*rows_ptr)[index1])
          tmp++;
        curIndex[index1] = tmp;

        cube_IJK = mesh->top_cube(i,j,k);
        center1  = get_cube_center(cube_IJK);
      }
  }

  void
  define_butting (index_t, index_t, index_t)
  {
  }
  void 
  add_boundary (index_t)
  {
  }

  item_t
  calc_tran_is_butting_by_x (index_t i, index_t j, index_t km, index_t kp)
  {
    const point_side_t &rightSide = get_side( right,cube_IJK,0);
    const ijk_cube_t &cube2       = mesh->top_cube(i+1,j,km);

    fpoint3d sqRS     = get_projection_on_all_axis_for_one_side(rightSide);
    fpoint3d D1       = fpoint3d::abs(get_cube_center(rightSide)-center1);
    fpoint3d center2  = get_cube_center(cube2);

    item_t tran = mesh->calculate_tran(index_ijk, ((i+1)+j*nx+kp*nx*ny), sqRS, D1, rightSide, center2, along_dim1);

    return tran;
  }

  item_t
  calc_tran_is_butting_by_y (index_t i, index_t j, index_t km, index_t kp)
  {
    const point_side_t &bottomSide = get_side( bottom,cube_IJK,0);
    const ijk_cube_t &cube2 = mesh->top_cube(i,j+1,kp);

    fpoint3d sqBS     = get_projection_on_all_axis_for_one_side(bottomSide);
    fpoint3d center2  = get_cube_center(cube2);
    fpoint3d D1       = fpoint3d::abs(get_cube_center(bottomSide)-center1);

    item_t tran = mesh->calculate_tran(index_ijk, (i+(j+1)*nx+kp*nx*ny), sqBS, D1, bottomSide,center2,along_dim2);

    return tran;
  }

  item_t
  calc_tran_by_x (index_t i, index_t j, index_t km, index_t kp)
  {
    const ijk_cube_t &cube2       = mesh->top_cube (i + 1, j, kp);
    const point_side_t &leftSide  = get_side (left, cube2, 0);
    const point_side_t &rightSide = get_side (right, cube_IJK, 0);

    fpoint3d centerSide1  = get_cube_center(rightSide);
    fpoint3d D1           = fpoint3d::abs(centerSide1-center1);
    fpoint3d A            = mesh->get_side_crossing_projection_on_all_axis(leftSide, rightSide);
    fpoint3d center2      = get_cube_center(cube2);

    item_t tran = mesh->calculate_tran(index_ijk, ((i+1)+j*nx+kp*nx*ny), A,D1, leftSide,center2,along_dim1);
    return tran;
  }

  item_t 
  calc_tran_by_y (index_t i, index_t j, index_t km, index_t kp)
  {
    const ijk_cube_t &cube2         = mesh->top_cube(i,j+1,kp);
    const point_side_t &topSide     = get_side( top, cube2,0);
    const point_side_t &bottomSide  = get_side( bottom,cube_IJK,0);

    fpoint3d centerSide1  = get_cube_center(bottomSide);
    fpoint3d D1           = fpoint3d::abs(centerSide1-center1);
    fpoint3d A            = mesh->get_side_crossing_projection_on_all_axis(topSide, bottomSide);
    fpoint3d center2      = get_cube_center(cube2);

    item_t tran = mesh->calculate_tran(index_ijk, (i+(j+1)*nx+kp*nx*ny), A, D1, topSide,center2,along_dim2);

    return tran;
  }

  item_t 
  calc_tran_is_butting_by_z (index_t i, index_t j, index_t km, index_t kp)
  {
    const point_side_t &upperSide = get_side( upper,cube_IJK,0);
    const ijk_cube_t &cube2       = mesh->top_cube(i,j,kp);

    fpoint3d center2  = get_cube_center(cube2);
    fpoint3d sqUS     = get_projection_on_all_axis_for_one_side(upperSide);
    fpoint3d D1       = fpoint3d::abs(get_cube_center(upperSide)-center1);

    item_t tran = mesh->calculate_tran(index_ijk, (i + j * nx + kp * nx * ny), sqUS, D1, upperSide,center2,along_dim3);

    return tran;
  }

  item_t
  calc_tran_by_z (index_t i, index_t j, index_t km, index_t kp)
  {
    const point_side_t &upperSide = get_side(upper,cube_IJK,0);
    const ijk_cube_t &cube2       = mesh->top_cube(i,j,kp);

    fpoint3d D1       = fpoint3d::abs(get_cube_center(upperSide)-center1);
    fpoint3d A        = get_projection_on_all_axis_for_one_side(upperSide);//because it's too expensive calculate as triangles...
    fpoint3d center2  = get_cube_center(cube2);

    item_t tran = mesh->calculate_tran(index_ijk, (i+j*nx+kp*nx*ny), A, D1, upperSide,center2,along_dim3);

    return tran;
  }

  void
  change_by_x (index_t i, index_t j, index_t km, index_t kp, bool is_butting)
  {
    //change jacobian
    index_t index1 = mesh->XYZ_to_inside (i,    j, km);
    index_t index2 = mesh->XYZ_to_inside (i + 1,j, kp);

    change_col_add(index1, index2,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, true,true,loop->con_num,*rows_ptr);
    if (curIndex[index2] == (*rows_ptr)[index2])
      curIndex[index2]++;
    change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus, false,true,loop->con_num,*rows_ptr);

    item_t tran = is_butting ? calc_tran_is_butting_by_x (i, j, km, kp) : calc_tran_by_x (i, j, km, kp);

    //change flux_connection
    cols_ind_transmis[loop->con_num++] = index1;
    cols_ind_transmis[loop->con_num++] = index2;
    values_transmis.push_back(tran);
    values_transmis.push_back(-tran);
  }

  void
  change_by_y (index_t i, index_t j, index_t km, index_t kp, bool is_butting)
  {
    index_t index1 = mesh->XYZ_to_inside (i, j,     km);
    index_t index2 = mesh->XYZ_to_inside (i, j + 1, kp);

    change_col_add(index1, index2, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus,true,true,loop->con_num,*rows_ptr);
    if (curIndex[index2] == (*rows_ptr)[index2])
      curIndex[index2]++;
    change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus,false,true,loop->con_num,*rows_ptr);

    item_t tran = is_butting ? calc_tran_is_butting_by_y (i, j, km, kp) : calc_tran_by_y (i, j, km, kp);

    //change flux_connection
    cols_ind_transmis[loop->con_num++] = index1;
    cols_ind_transmis[loop->con_num++] = index2;
    values_transmis.push_back(tran);
    values_transmis.push_back(-tran);
  }

  void
  change_by_z (index_t i, index_t j, index_t km, index_t kp, bool is_butting)
  {
    //change jacobian
    index_t index1 = mesh->XYZ_to_inside (i, j, km);
    index_t index2 = mesh->XYZ_to_inside (i, j, kp);

    change_col_add(index1, index2,*cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus,true,true,loop->con_num,*rows_ptr);
    if (curIndex[index2] == (*rows_ptr)[index2])
      curIndex[index2]++;
    change_col_add(index2, index1, *cols_ind,curIndex, matrix_block_idx_minus, matrix_block_idx_plus,false,true,loop->con_num,*rows_ptr);

    item_t tran = is_butting ? calc_tran_is_butting_by_z (i, j, km, kp) : calc_tran_by_z (i, j, km, kp);

    //change flux_connection
    cols_ind_transmis[loop->con_num++] = index1;
    cols_ind_transmis[loop->con_num++] = index2;
    values_transmis.push_back(tran);
    values_transmis.push_back(-tran);
  }

  mesh_t              *mesh;
  loop_t              *loop;
  index_array_t       *rows_ptr;
  index_array_t       *cols_ind;
  index_array_t       &cols_ind_transmis;
  rhs_item_array_t    &values_transmis;
  index_array_t       &matrix_block_idx_minus;
  index_array_t       &matrix_block_idx_plus;

  index_t             nx;
  index_t             ny;
  index_t             nz;

  index_array_t       curIndex;

  const cube_index_t  *cube_IJK_index;
  index_t             index1;
  index_t             index_ijk;

  ijk_cube_t          cube_IJK;
  fpoint3d            center1;
};

template <typename M, typename L, typename RP, typename CI, typename CT, typename FC>
build_jacobian_cols_class <M, L>
build_jacobian_cols (mesh_grdecl <M> *m, L *l, RP *rp, CI *ci, CT &conn_trans, FC &flux_conn)
{
  return build_jacobian_cols_class <M, L> (m, l, rp, ci, 
    conn_trans->get_cols_ind (), conn_trans->get_values (),
    flux_conn->get_matrix_block_idx_minus (), flux_conn->get_matrix_block_idx_plus ());
}

template <typename strategy_t>
struct build_jacobian_and_flux : boost::noncopyable
{
  typedef typename strategy_t::index_t          index_t;
  typedef typename strategy_t::item_t           item_t;
  typedef typename strategy_t::index_array_t    index_array_t;
  typedef typename strategy_t::item_array_t     item_array_t;

  typedef mesh_grdecl <strategy_t>              mesh_t;
  typedef typename mesh_t::side_t               side_t;
  typedef typename mesh_t::cube_index_t         cube_index_t;

  build_jacobian_and_flux (mesh_grdecl <strategy_t> *mesh)
  : mesh (mesh)
  , nx (mesh->nx)
  , ny (mesh->ny)
  , nz (mesh->nz)
  , con_num (0)
  {
    is_butting.assign (nx * ny, false);
  }

  template <typename loop_body_t>
  void
  cell_loop (loop_body_t loop_body)
  {
    for (index_t i = 0; i < nx; ++i)
      {
        for (index_t j = 0; j < ny; ++j)
          {
            index_t bottom_Kx = 0;
            index_t bottom_Ky = 0;

            index_t top_Kx, top_Ky, top_Kz;
            for (index_t k = 0; k < nz; ++k)
              {
                index_t current_cell_index  = i+j*nx+k*nx*ny;
                index_t external_cell_index = mesh->XYZ_to_inside(i,j,k);
                //skip non-active cells but use it for adding new boundary cells
                bool is_current_cell_active = mesh->sp_actnum[current_cell_index];

                const cube_index_t &cube_IJK_index = top_cube_index (i, j, k, nx, ny);
                bool flag_is_boundary = false;
                if (is_current_cell_active && (i == 0 || j == 0 || k == 0 || i == nx-1 || j == ny-1 || k == nz-1))
                  {
                    loop_body.add_boundary (external_cell_index);
                    flag_is_boundary = true;
                  }

                loop_body.prepare (i, j, k, &cube_IJK_index, is_current_cell_active);
                loop_body.define_butting (i, j, k);

                //look only 3 positive-direction side
                if (is_butting[j*nx+i])
                  {
                    //speed-up by calculating just one square
                    if (i+1 < nx)
                      {
                        if (is_current_cell_active && mesh->sp_actnum[(i+1)+j*nx + k*nx*ny])
                          {
                            loop_body.change_by_x (i, j, k, k, true);
                          }
                        else if (!flag_is_boundary)
                          {
                            if (is_current_cell_active)
                              {
                                loop_body.add_boundary (external_cell_index);
                                flag_is_boundary = true;
                              }
                            else
                              {
                                loop_body.add_boundary (mesh->XYZ_to_inside(i + 1, j, k));
                              }
                          }
                      }
                    if (j+1 < ny)
                      {
                        if (is_current_cell_active && mesh->sp_actnum[i+(j+1)*nx + k*nx*ny])
                          {
                            loop_body.change_by_y (i, j, k, k, true);
                          }
                        else if (!flag_is_boundary)
                          {
                            if (is_current_cell_active)
                              {
                                loop_body.add_boundary (external_cell_index);
                                flag_is_boundary = true;
                              }
                            else
                              {
                                loop_body.add_boundary (mesh->XYZ_to_inside(i,j+1,k));
                              }
                          }
                      }
                    if (is_current_cell_active && (k+1 < nz))
                      {
                        if (mesh->sp_actnum[i+j*nx + (k+1)*nx*ny])
                          {
                            loop_body.change_by_z (i, j, k, k + 1, true);
                          }
                      }
                    continue;
                  }
#pragma endregion

                //calculating for non butting cell!
                top_Kx = bottom_Kx-1;
#pragma region rightSide
                if (i+1 < nx)
                  {
                    //looking for incident block at right side and their zmin(kmin) & zmax(kmax)
                    const side_t &right_side = get_side( right, cube_IJK_index,0);
                    item_t maxZ_RS = mesh->zcorn_array[*std::max_element(right_side.begin(),right_side.end())];

                    //go to top
                    //looking for top block for right side - minZ_LS <= maxZ_RS
                    while (top_Kx < nz-1)
                      {
                        top_Kx++;
                        if (flag_is_boundary && !mesh->sp_actnum[i+1+j*nx + top_Kx*nx*ny])//skip non-active cells
                          continue;

                        const side_t &temp_left_side = get_side( left, top_cube_index(i+1,j,top_Kx, nx, ny),0);
                        item_t minZ_LS = mesh->zcorn_array[*std::min_element(temp_left_side.begin(),temp_left_side.end())];

                        if (minZ_LS >= maxZ_RS) //we can't go more
                          break;
                        else if (mesh->is_2side_crossing(right_side, temp_left_side, true))// => [i,j,k] && [i+1,j,top_Kx] have common contact
                          {
                            if (is_current_cell_active)
                              {
                                if (mesh->sp_actnum[(i+1)+j*nx + top_Kx*nx*ny])
                                  {
                                    loop_body.change_by_x (i, j, k, top_Kx, false);
                                  }
                                else if (!flag_is_boundary)
                                  {
                                    loop_body.add_boundary (external_cell_index);
                                    flag_is_boundary = true;
                                  }
                              }
                            else if (mesh->sp_actnum[(i+1)+j*nx + top_Kx*nx*ny])
                              {
                                loop_body.add_boundary (mesh->XYZ_to_inside(i+1,j,top_Kx));
                              }
                          }
                      }
                  }

#pragma endregion

                top_Ky = bottom_Ky-1;
#pragma region bottomSide
                if (j+1 < ny)
                  {
                    const side_t &bottom_side = get_side( bottom,cube_IJK_index,0);
                    item_t maxZ_BS = mesh->zcorn_array[*std::max_element(bottom_side.begin(),bottom_side.end())];

                    //go to top
                    //looking for top block for bottom side - minZ_TS <= maxT_RS
                    while (top_Ky < nz-1)
                      {
                        top_Ky++;
                        if (flag_is_boundary && !mesh->sp_actnum[i+(j+1)*nx + top_Ky*nx*ny])//skip non-active cells
                          continue;

                        const side_t &temp_top_side = get_side( left, top_cube_index(i,j+1,top_Ky, nx, ny),0);
                        item_t minZ_TS = mesh->zcorn_array[*std::min_element(temp_top_side.begin(),temp_top_side.end())];
                        if (minZ_TS >= maxZ_BS) //we can't go more
                          break;
                        else if (mesh->is_2side_crossing(bottom_side, temp_top_side,!true))// => [i,j,k] && [i,j+1,top_Ky] have common contact
                          {
                            if (is_current_cell_active) //standart case
                              {
                                if (mesh->sp_actnum[i+(j+1)*nx + top_Ky*nx*ny])
                                  {
                                    loop_body.change_by_y (i, j, k, top_Ky, false);
                                  }
                                else if (!flag_is_boundary)
                                  {
                                    loop_body.add_boundary (external_cell_index);
                                    flag_is_boundary = true;
                                  }
                              }
                            else if (mesh->sp_actnum[i+(j+1)*nx + top_Ky*nx*ny])
                              {
                                loop_body.add_boundary (mesh->XYZ_to_inside(i,j+1,top_Ky));
                              }
                          }
                      }
                  }

#pragma endregion

#pragma region upperSide
                if (k+1 < nz && is_current_cell_active)
                  {
                    top_Kz = k+1;
                    while (top_Kz < nz)
                      {
                        if (!mesh->sp_actnum[i+j*nx+top_Kz*nx*ny])
                          {
                            top_Kz++;
                            continue;
                          }
                        break;
                      }
                    if (top_Kz == nz) //no connections
                      continue;
                    else if (mesh->sp_actnum[i+j*nx+top_Kz*nx*ny])
                      {
                        // => [i,j,k] && [i,j,k_Topz] have common contact (take square upperSide) -> findAreaOfSide(upperSide[0],upperSide[1],upperSide[2],upperSide[3]);
                        loop_body.change_by_z (i, j, k, top_Kz, false);
                      }
                  }
#pragma endregion
                bottom_Kx = top_Kx;
                bottom_Ky = top_Ky;

              }
          }
      }
  }

  mesh_grdecl <strategy_t>    *mesh;
  index_t                     nx;
  index_t                     ny;
  index_t                     nz;
  shared_vector <bool>        is_butting;
  index_t                     con_num;
};


template<class strategy_t>
int mesh_grdecl<strategy_t>::build_jacobian_and_flux_connections_add_boundary (const sp_bcsr_t &jacobian, 
                                                                               const sp_flux_conn_iface_t &flux_conn,
                                                                               index_array_t &boundary_array)
{
  n_connections = 0;

  index_t max_size = n_active_elements;
  jacobian->get_cols_ind().clear();
  jacobian->get_rows_ptr().clear();
  jacobian->init_struct(max_size,max_size, max_size);
  index_array_t* rows_ptr = &jacobian->get_rows_ptr();
  
  sp_bcsr_t conn_trans;

  (*rows_ptr)[0] = 0;

  std::vector<bool> is_butting(nx*ny,false);
  int nButting = 0;

  std::set<index_t, std::less<index_t> > boundary_set;

  build_jacobian_and_flux <strategy_t> build_jacobian (this);

  //first step - define and fill - rows_ptr (jacobian)
  build_jacobian.cell_loop (build_jacobian_rows (this, &build_jacobian, boundary_set, rows_ptr));

  //////jacobian//////////////////////
  //sum rows_ptr
  for (size_t i = 1; i < rows_ptr->size(); ++i)
    {
      (*rows_ptr)[i]++;
      (*rows_ptr)[i] += (*rows_ptr)[i-1];
    }
  //create cols_ind
  index_array_t* cols_ind = &jacobian->get_cols_ind();

  cols_ind->resize((*rows_ptr)[rows_ptr->size()-1],-1);
  ////////transmis/////////////////////////
  index_t cols_ind_n = (index_t)cols_ind->size();


  index_t con_num = (cols_ind_n-max_size)/2;//connection number
  n_connections = con_num;
  conn_trans = flux_conn->get_conn_trans();
  conn_trans->init_struct(con_num,2*con_num,2*con_num);

  index_array_t *rows_ptr_transmis = &conn_trans->get_rows_ptr();
  index_array_t *cols_ind_transmis = &conn_trans->get_cols_ind();
  rhs_item_array_t *values_transmis = &conn_trans->get_values();
  index_array_t &matrix_block_idx_minus = flux_conn->get_matrix_block_idx_minus ();
  index_array_t &matrix_block_idx_plus = flux_conn->get_matrix_block_idx_plus ();

  matrix_block_idx_minus.resize(con_num*2,-1);
  matrix_block_idx_plus.resize(con_num*2,-1);

  if (!con_num)
    {
      for (int i = 0; i < cols_ind_n; ++i)
        (*cols_ind)[i] = i;
      return 0;
    }

  for (index_t i = 0; i < con_num+1; ++i)
    (*rows_ptr_transmis)[i] = i*2;

  //second step - fill and define cols_ind
  build_jacobian.con_num = 0;
  build_jacobian.cell_loop (build_jacobian_cols (this, &build_jacobian, rows_ptr, cols_ind, conn_trans, flux_conn));


  boundary_array.assign(boundary_set.begin(), boundary_set.end());
  return (int) boundary_array.size();
}

BS_INST_STRAT(mesh_grdecl);
