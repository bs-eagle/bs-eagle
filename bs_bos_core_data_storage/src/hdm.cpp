#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm.h"
#include "main_def.h"
#include "arrays.h"
#include "localization.h"
#include "arrays_tables.h"
#include "data_class.h"

#include "write_time_to_log.h"

      #include <boost/multi_index_container.hpp>
    #include <boost/multi_index/ordered_index.hpp>
    #include <boost/multi_index/identity.hpp>
    #include <boost/multi_index/member.hpp>
    #include "boost/date_time/posix_time/posix_time.hpp"
    
    #include <iostream>

    using namespace boost;
    using namespace boost::multi_index;    
    using namespace boost::posix_time;
    using namespace boost::gregorian;
    using namespace std;
 
/////////////////////////////////////////////////////////////////////////////
// You can add new variable just find next Words
// for symple variable                  -- VAR_V
// for array                                    -- ARRAY_V
// for strings                                  -- STRING_V

//! Default minimal DX value allowed for active cells
#define DEFAULT_MINIMAL_DX              1.e-6
//! Default minimal DY value allowed for active cells
#define DEFAULT_MINIMAL_DY              1.e-6
//! Default minimal DZ value allowed for active cells
#define DEFAULT_MINIMAL_DZ              1.e-10
//! Default minimal PERM value allowed for active cells
#define DEFAULT_MINIMAL_PERM            1e-16
#define DEFAULT_MINIMAL_NTG             1.e-8
#define DEFAULT_MINIMAL_TOPS            1.e-8
//! default blank value
#define DEFAULT_BLANK_VALUE             1.e43
//! check blank value
#define CHECK_BLANK_VALUE(X)  (fabs ((X) - DEFAULT_BLANK_VALUE) < 10)

//! Default minimal PORO value allowed for active cells
#define DEFAULT_MINIMAL_PORO                        1.e-10
#define DEFAULT_MINIMAL_PORO_VOLUME_MULT            1.e-12
/////////////////////////////////////////////////////////////////////////////


// Please don't delete commented old style outputs .. like rep->blah blah blah
namespace blue_sky
  {
  
#define DM_ASSERT_EXCEPTION \
    BS_ASSERT(false) (out_s.str());\
    throw bs_exception("Hydrodynamic model class",out_s.str().c_str());

  
  hdm::~hdm ()
  {

  }

  
  hdm::hdm(bs_type_ctor_param /*param*/): lkeeper ("C", LC_ALL)
  {
    this->reader = BS_KERNEL.create_object(FRead::bs_type());
    this->data = BS_KERNEL.create_object(idata::bs_type());
    this->km = BS_KERNEL.create_object(keyword_manager::bs_type());
    
    this->pvt_3p_ = BS_KERNEL.create_object ("pvt_3p");
    this->event_manager_ = BS_KERNEL.create_object ("event_manager");
  }

  hdm::hdm(const hdm& src):bs_refcounter (src), lkeeper ("C", LC_ALL)
  {
    *this = src;
  }
 
  void
  hdm::init_fluids(t_int n_scal_regions, t_int n_pvt_regions)
  {
    int n_phases; 
    
    pvt_3p_->init_pvt_arrays (n_pvt_regions, data->props->get_b("oil_phase"),
                                             data->props->get_b("water_phase"),
                                             data->props->get_b("gas_phase"));
    
    n_phases = data->props->get_b("oil_phase");
    n_phases += data->props->get_b("water_phase");
    n_phases += data->props->get_b("gas_phase");
    if (n_phases > 1)
      {
        this->scal_3p_ = BS_KERNEL.create_object ("scal_3p");
        scal_3p_->init_scal_input_table_arrays (n_scal_regions, data->props->get_b("oil_phase"),
                                                                data->props->get_b("water_phase"),
                                                                data->props->get_b("gas_phase"));
      }
  }
  
   void
  hdm::init(const std::string &model_name)
  {
    smart_ptr <hdm_iface, true> hdm = this;
    keyword_params kp;
    int n_pases;
    
    kp.hdm = this;
    data->h5_pool->open_file (model_name + ".h5", "/pool");
    km->init(hdm);
    
    
    switch (data->props->get_i("mesh"))
    {
      case 0:
        km->handle_keyword_reactor ("MESH_IJK", kp);
        break;
      case 1:
        km->handle_keyword_reactor ("MESH_GRDECL", kp);
        break;   
          
     default:
        bs_throw_exception ("init: wrong mesh choice");  
    }
    
    switch (data->props->get_i("init"))
    {
      case 0:
        km->handle_keyword_reactor ("EXPLICIT_MODEL", kp);
        break;
      case 1:
        km->handle_keyword_reactor ("EQUIL_MODEL", kp);
        break;   
          
     default:
        bs_throw_exception ("init: wrong mesh choice");  
    }
    
  }
  
  
  void 
  hdm::read_keyword_file(const std::string filename)
  {
    char buf[CHAR_BUF_LEN];
    char key[CHAR_BUF_LEN];
    keyword_params kp;
    int flag;
    int len;
    
    write_time_to_log init_time ("Read model", "");

    reader->init (filename, filename);
    kp.hdm = this;

    // start of loop for data file reading
    flag = 1;

    for (; flag;)
      {
        // reading keyword
        len = reader->read_line (buf, CHAR_BUF_LEN);
        if (len < 0)              // break if EOF
          {
            switch (len)
              {
              case YS_NOTHING_TO_READ:
                break;
              default:
                return;
              }
            //rep->print (LOG_READ_SECTION, LOG_DEBUG, "Finish reading\n");
            BOSOUT (section::read_data, level::low) << "Finish reading" << bs_end;
            break;
          }
        if (sscanf (buf, "%s", key) != 1) // look from buffer keyword
          continue;

        std::string keywrd(key);

        if (key[0] == '/')
          {
            continue;
          }

        if (std::string (key) == "END")
          {
            BOSOUT (section::read_data, level::low) << "Finish reading with END keyword" << bs_end;
            break;
          }

        km->handle_keyword (keywrd, kp);
      }
  }


  //! macros for checking
#define CHECK_FOR_TWO_PARAMS(VAL,ACTNUM,MIN_DEF,COUNTER)                                \
  if ((VAL) < (MIN_DEF) || CHECK_BLANK_VALUE(VAL))                                      \
    {                                                                                   \
      if ((ACTNUM))                                                                     \
        {                                                                               \
          t_long ix, iy, iz;                                                               \
          iz = i / (nx * ny);                                 \
          iy = (i - iz * (nx * ny)) / nx;        \
          ix = i - iz * (nx * ny) - iy * nx;     \
          BOSWARN (section::check_data, level::low)                                     \
            << "blocks[" << (ix + 1) << ", " << (iy + 1) << ", " << (iz + 1)            \
            << "] will be set inactive because of " << #VAL << " = " << VAL << bs_end;  \
          (ACTNUM) = 0;                                                                 \
          ++(COUNTER);                                                                  \
          continue;                                                                     \
        }                                                                               \
    }

#define CHECK_FOR_THREE_PARAMS(VAL,ACTNUM,MIN_DEF,MAX_DEF,COUNTER)                      \
  if ((VAL) < (MIN_DEF) || CHECK_BLANK_VALUE(VAL) || (VAL) > (MAX_DEF))                 \
    {                                                                                   \
      if ((ACTNUM))                                                                     \
        {                                                                               \
          int ix, iy, iz;                                                               \
          iz = i / (nx * ny);                                                           \
          iy = (i - iz * (nx * ny)) / nx;                                               \
          ix = i - iz * (nx * ny) - iy * nx;                                            \
          BOSWARN (section::check_data, level::low)                                     \
            << "blocks[" << (ix + 1) << ", " << (iy + 1) << ", " << (iz + 1)            \
            << "] will be set inactive because of " << #VAL << bs_end;                  \
          (ACTNUM) = 0;                                                                 \
          ++(COUNTER);                                                                  \
          continue;                                                                     \
        }                                                                               \
    }

#define CHECK_DIRECTION_FOR_XYZ(case,num)   \
  if ((case) && actnum [(num)])             \
    {                                       \
      n++;                                  \
      valx += dx[(num)];                    \
      valy += dy[(num)];                    \
      valz += dz[(num)];                    \
    }

#define CHECK_DIRECTION_FOR_TOPS(case,num)  \
  if ((case) && actnum[(num)])              \
    {                                       \
      n++;                                  \
      valtops += tops[(num)];               \
    }

  /*
  
  void hdm::update_geometry() const
    {
      int i, ix, iy, iz, n;
      std::ostringstream out_s;
      int nb = data->nx * data->ny * data->nz;

      int n_tops, n_xyz;
      double dx_var, dy_var, dz_var, tops_var;
      double valx, valy, valz, valtops;

      array_float16_t dx = tools::get_non_empty (((*data->d_map)[DX].array));
      array_float16_t dy = tools::get_non_empty (((*data->d_map)[DY].array));
      array_float16_t dz = tools::get_non_empty (((*data->d_map)[DZ].array));
      array_float16_t tops = tools::get_non_empty (((*data->d_map)[TOPS].array));
      array_uint8_t   actnum = tools::get_non_empty (((*data->i_map)[ACTNUM].array));

      // In case of restart we got all geometry from mesh restored from
      // the restart file. So we skip this test.
      if (data->restart)
        return;

      if (!nb)
        {
          bs_throw_exception ("update_geometry: all dimensions must be greater than 0");
        }

      dx_var = dy_var = dz_var = tops_var = 0.;
      valx = valy = valz = valtops = 0.;
      n_tops = n_xyz = 0;
      if (data->geom_def_flag == GEOM_FLAG_DX_DY_DZ)
        {
          for (iz = 0, i = 0; iz < data->nz; ++iz)
            {
              for (iy = 0; iy < data->ny; ++iy)
                {
                  for (ix = 0; ix < data->nx; ++ix, ++i)
                    {
                      if (actnum[i])
                        {
                          dx_var += dx[i];
                          dy_var += dy[i];
                          dz_var += dz[i];
                          n_xyz++;
                        }
                    }
                }
            }

          if (!n_xyz)
            {
              bs_throw_exception ("update_geometry: No active blocks");
            }

          for (iy = 0, i = 0; iy < data->ny; ++iy)
            {
              for (ix = 0; ix < data->nx; ++ix, ++i)
                {
                  if (actnum[i])
                    {
                      tops_var += tops[i];
                      ++n_tops;
                    }
                }
            }
          if (!n_tops)
            {
              bs_throw_exception ("update_geometry: No active blocks");
            }

          dx_var /= n_xyz;
          dy_var /= n_xyz;
          dz_var /= n_xyz;
          tops_var /= n_tops;
          for (iy = 0, i = 0; iy < data->ny; ++iy)
            {
              for (ix = 0; ix < data->nx; ++ix, ++i)
                {
                  if (!actnum[i])
                    {
                      valtops = 0.;
                      n = 0;
                      CHECK_DIRECTION_FOR_TOPS (ix > 0, i - 1);
                      CHECK_DIRECTION_FOR_TOPS (iy > 0, i - data->nx);
                      CHECK_DIRECTION_FOR_TOPS (ix < data->nx - 1 , i + 1);
                      CHECK_DIRECTION_FOR_TOPS (iy < data->ny - 1 , i + data->nx);
                      if (n)
                        {
                          tops[i] = valtops / n;
                        }
                      else
                        {
                          tops[i] = tops_var;

                          BOSWARN (section::check_data, level::low) 
                            << "new tops in node (" << (ix + 1)
                            << ", " << (iy + 1) << ", " << (iz + 1) << ") = " << tops[i]
                            << " (average = " << tops_var << ")" << bs_end;
                        }
                    }
                }
            }

          for (iz = 0, i = 0; iz < data->nz; ++iz)
            {
              for (iy = 0; iy < data->ny; ++iy)
                {
                  for (ix = 0; ix < data->nx; ++ix, ++i)
                    {
                      if (!actnum[i])
                        {
                          valx = valy = valz = 0.;
                          n = 0;
                          CHECK_DIRECTION_FOR_XYZ (ix > 0, i - 1);
                          CHECK_DIRECTION_FOR_XYZ (iy > 0, i - data->nx);
                          CHECK_DIRECTION_FOR_XYZ (iz > 0, i - data->nx * data->ny );
                          CHECK_DIRECTION_FOR_XYZ (ix < data->nx - 1 , i + 1);
                          CHECK_DIRECTION_FOR_XYZ (iy < data->ny - 1 , i + data->nx);
                          CHECK_DIRECTION_FOR_XYZ (iz < data->nz - 1 , i + data->nx * data->ny);
                          if (n)
                            {
                              dx[i] = valx / n;
                              dy[i] = valy / n;
                              dz[i] = valz / n;
                            }
                          else
                            {
                              dx[i] = dx_var;
                              dy[i] = dy_var;
                              dz[i] = dz_var;
                            }

                          BOSOUT (section::check_data, level::low)
                            << "Info: new dx in node (" << (ix + 1) << ", "
                              << (iy + 1) << ", " << (iz + 1) << ") = " << dx[i] << " (aver dx = "
                              << dx_var << ")" 
                              << "\n"
                            << "Info: new dy in node (" << (ix + 1) << ", "
                              << (iy + 1) << ", " << (iz + 1) << ") = " << dy[i] << " (aver dy = "
                              << dy_var << ")" 
                              << "\n"
                            << "Info: new dz in node (" << (ix + 1) << ", "
                              << (iy + 1) << ", " << (iz + 1) << ") = " << dz[i] << " (aver dz = "
                              << dz_var << ")" 
                              << bs_end;
                        }
                    }
                }
            }
        }
      else if (data->geom_def_flag == GEOM_FLAG_ZCORN_COORD)
        {
#if 0
          for (iz = 0, i = 0; iz <= nz; ++iz)
            {
              for (iy = 0; iy <= ny; ++iy)
                {
                  for (ix = 0; ix <= nx; ++ix, ++i)
                    {
                      if ((*depth)[i] < 0)
                        {
                          BOSERR (section::check_data, level::error) << "Depth of node (" << (ix + 1) << ", "
                          << (iy + 1) << ", " << (iz + 1) << ") = " << depth[i]
                          << " is out of range" << bs_end;
                          return YS_BAD_VALUE;
                        }
                    }
                }
            }
#endif
        }
      else
        {
          bs_throw_exception ("Unsupported geometry specification");
        }
    }
*/
  
  void hdm::check_arrays_for_inactive_blocks () const
    {
      spv_int actnum_  = data->get_i_array ("ACTNUM");
      spv_float permx_ = data->get_fp_array ("PERMX");
      spv_float permy_ = data->get_fp_array ("PERMY");
      spv_float permz_ = data->get_fp_array ("PERMZ");
      spv_float poro_  = data->get_fp_array ("PORO");
      spv_float ntg_   = data->get_fp_array ("NTG");
      
      t_int* actnum = actnum_->data ();
      const t_float *permx  = permx_->data ();
      const t_float *permy  = permy_->data ();
      const t_float *permz  = permz_->data ();
      const t_float *poro   = poro_->data ();
      const t_float *ntg    = ntg_->data ();

      t_long nx = data->props->get_i("nx");
      t_long ny = data->props->get_i("ny");
      t_long nz = data->props->get_i("nz");
      t_long nb = nx * ny * nz;

      t_long permx_counter = 0;
      t_long permy_counter = 0;
      t_long permz_counter = 0;
      t_long poro_counter = 0;
      t_long ntg_counter = 0;
      for (t_long i = 0; i < nb; ++i)
        {
          if (!actnum[i])
            {
              continue;
            }

          if (actnum[i] && permx[i] < DEFAULT_MINIMAL_PERM && permy[i] < DEFAULT_MINIMAL_PERM && permz[i] < DEFAULT_MINIMAL_PERM)
            {
              t_long iz = i / (nx * ny);
              t_long iy = (i - iz * (nx * ny)) / nx;
              t_long ix = i - iz * (nx * ny) - iy * nx;

              BOSWARN (section::check_data, level::low)
                << "blocks " << i << " [" << (ix + 1) << ", " << (iy + 1) << ", "
                << (iz + 1) << "] will be set inactive because of PERMEABILITY" << bs_end;

              actnum[i] = 0;
              ++permx_counter;
              continue;
            }

          CHECK_FOR_TWO_PARAMS (poro[i], actnum[i], DEFAULT_MINIMAL_PORO, poro_counter);
          if (ntg_->size () != 0)
            {
              CHECK_FOR_TWO_PARAMS (ntg[i], actnum[i], DEFAULT_MINIMAL_NTG, ntg_counter);
            }
        }

      if (permx_counter)
        BOSWARN (section::check_data, level::warning) << permx_counter << " blocks will be set inactive because of PERMX" << bs_end;
      if (permy_counter)
        BOSWARN (section::check_data, level::warning) << permy_counter << " blocks will be set inactive because of PERMY" << bs_end;
      if (permz_counter)
        BOSWARN (section::check_data, level::warning) << permz_counter << " blocks will be set inactive because of PERMZ" << bs_end;
      if (poro_counter)
        BOSWARN (section::check_data, level::warning) << poro_counter << " blocks will be set inactive because of PORO" << bs_end;
      if (ntg_counter)
        BOSWARN (section::check_data, level::warning) << ntg_counter << " blocks will be set inactive because of NTG" << bs_end;
    }
    
 
     
  void hdm::test_well_storage ()
  {
     struct well_info
      {
        std::string               name;
        double  date;
        
        int i_params[20];
        float f_params[20];
        
      };
    
    struct comp_date
    {
      bool operator()(const well_info& w1,const well_info& w2)const 
      {
        if (w1.name < w2.name)
         return true;
        else if (w1.name == w2.name)
          return (w1.date < w2.date);
        else
          return false;
      } // && w1.date < w2.date
    }; 
    
    typedef multi_index_container<
      well_info,
      indexed_by<
        ordered_unique<identity<well_info>, comp_date>,
        ordered_non_unique<member<well_info, double, &well_info::date>>
        //ordered_non_unique<member<well_info, int, &well_info::id>>
      > 
    > well_set;
    
    struct comp_name
      {
        bool operator()(std::string x,const well_info& w2)const{return x < w2.name;}

        bool operator()(const well_info& w1,std::string x)const{return w1.name < x;}
      };
    
    
    
    
    well_set ws;
    well_info *iw;
    int i, j, k, n_wells, n_dates;
    days dd(1);
    char tmp[10];
    std::vector <double> res_date;
    std::vector <float> res_param;
    
    
    //typedef well_set::index<name>::type employee_set_by_name;
    
    n_wells = 5000;
    n_dates = 500;
    
    srand (555);
    
    printf ("Fill multiindex\n");
    cout << "Fill started!"; 
    std::pair<well_set::nth_index<0>::type::iterator, bool> res;
    
    iw = new well_info;
    iw->name = "1";
    iw->date = 2;
    iw->i_params[0] = 500;
    res = ws.insert (*iw);
    printf ("Insert: %i\n", res.second);
    
    iw = new well_info;
    iw->name = "1";
    iw->date = 2;
    iw->i_params[0] = 1500;
    res = ws.insert (*iw);
    printf ("Insert: %i\n", res.second);
    
    iw = new well_info;
    iw->name = "1";
    iw->date = 3;
    iw->i_params[0] = 2500;
    res = ws.insert (*iw);
    printf ("Insert: %i\n", res.second);
    
    
    /*
    {  
      write_time_to_log init_time("Fill!", "");
      for (i = 0; i < n_wells; i++)
        {
          iw = new well_info;
          sprintf(tmp, "%d", i);
          iw->name = tmp;
          for  (j = 0; j < n_dates; j++)
            {
              iw->date = rand () % n_wells;
              
              for  (k = 0; k < 20; k++)
                {
                  iw->i_params[k] = rand () % n_wells;
                  iw->f_params[k] = rand () % n_wells * 1.5;
                }
              ws.insert (*iw);
           }
        }
     }
     
    */
    
    const well_set::nth_index<0>::type& name_index = ws.get<0>();
    const well_set::nth_index<1>::type& date_index = ws.get<1>();
    well_set::nth_index<0>::type::iterator name_it, name_b, name_e;;
    well_set::nth_index<1>::type::iterator date_it, date_b, date_e;
     
    int pi = 0;
    float pf = 0.;
    long long n = 0;
      
    {
      write_time_to_log init_time("sum for T!", "");  
      printf ("Calculating sum of param for given T\n");
      date_b = date_index.lower_bound (333);
      date_e = date_index.upper_bound (333);
     
      for (date_it = date_b; date_it != date_e; ++date_it)
        {
          //printf ("%s \n", date_it->name.c_str());
          pi += date_it->i_params[15];
          pf += date_it->f_params[10];
          n++;
        }
    }
    
    printf ("Done %d iters: %d, %lf\n", n, pi, pf);
    
    pi = 0;
    pf = 0;
    
    {
      write_time_to_log init_time("Sum for name!", "");
        
      printf ("Done %d iters\n", n);
      printf ("Calculating sum of param for given well name\n");
      name_b = name_index.lower_bound ("1", comp_name());
      name_e = name_index.upper_bound ("1", comp_name());
      
      n = distance(name_b, name_e);
      res_date.reserve (n);
      res_param.reserve (n);
      printf ("Distance %d \n", n);
      n = 0;
      //printf ("%s %lf \n", name_b->name.c_str(), name_b->date);
      //printf ("%s %lf \n", name_e->name.c_str(), name_e->date);
      for (name_it = name_b; name_it != name_e; ++name_it)
        {
          pi += name_it->i_params[15];
          pf += name_it->f_params[10];
          res_date.push_back(name_it->date);
          res_param.push_back(name_it->f_params[10]);
          n++;
          printf ("%s %lf %d\n", name_it->name.c_str(), name_it->date, name_it->i_params[0]);
        }
      printf ("Done %d iters\n", n);
      printf ("Job finished\n");
    }
   
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(hdm)
  BLUE_SKY_TYPE_STD_COPY(hdm)

  BLUE_SKY_TYPE_IMPL(hdm, objbase, "hdm", "BOS_Core hdm class", "BOS_Core hdm class")
}//ns bs
