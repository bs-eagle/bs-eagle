/// @file sql_well_hdm.cpp
/// @brief sql_well HDM-related functions implementation
/// @author uentity
/// @version 1.0
/// @date 25.03.2014
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_common.h"
#include "sql_well.h"
#include "frac_comp_ident.h"
#include "wpi_iface.h"

#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/mpl/bool.hpp>
// DEBUG
//#include <boost/timer/timer.hpp>
//#include <google/profiler.h>

// shorthand macro to pass values via val2str filter
#define V2S(t) val2str::go(t).c_str()

namespace blue_sky {

// hidden details
namespace {

// the purpose of this helper is to print '1*' if argument == -1
// otherwise number converted to string is unchanged
struct val2str {
  // specialization for floating-point numbers
  template< class T >
  static std::string do_fix(const T& v, const boost::true_type) {
    if(v == -1.0)
      return "1*";
    else
      return boost::lexical_cast< std::string >(v);
  }

  // specialization for other types
  template< class T >
  static std::string do_fix(const T& v, const boost::false_type) {
      return boost::lexical_cast< std::string >(v);
  }

  // main operator
  template< class T >
  static std::string go(const T& v) {
    return do_fix(v, boost::is_floating_point< T >());
  }
};

} // eof hidden namespace

  int
  sql_well::save_to_bos_ascii_file (const std::string &fname, sp_pool_t pool, sp_prop_t prop)
    {
      FILE *fp = fopen (fname.c_str (), "w");
      char s_buf[2048];

      // interface to mesh
      sp_iwpi iwpi = BS_KERNEL.create_object("wpi_iface");
      BS_ASSERT(iwpi);

      //const double eps = 1e-10;
      sp_dt_t dt_t = BS_KERNEL.create_object ("dt_tools");

      if (!fp)
        {
          printf ("Error: Can not open destination file %s\n", fname.c_str ());
          return -1;
        }

      // obtain model dimensions
      boost::python::list dims;
      dims = pool->py_get_pool_dims();
      t_long nx = boost::python::extract<int>(dims[0]);
      t_long ny = boost::python::extract<int>(dims[1]);
      //nz = boost::python::extract<int>(dims[2]);
      //unsigned long nx_ny = nx * ny;
      BS_SP (well_pool_iface) sp_wp = this;

      // precalc trimesh backend to share among all compl & frac builders
      sp_obj trim_backend = iwpi->make_trimesh_backend(nx, ny,
        pool->get_fp_data("COORD"), pool->get_fp_data("ZCORN"),
        fci::strat_t::traits_t::name()
      );

      // here we wil store builders on every date
      typedef std::list< fci::compl_n_frac_builder > cfb_storage_t;
      typedef cfb_storage_t::iterator cfb_iterator;
      typedef cfb_storage_t::const_iterator cfb_citerator;
      cfb_storage_t cfb_storage;

      // get dates list
      std::string sql = "SELECT d FROM dates ORDER BY d ASC";
      if (prepare_sql (sql.c_str ()))
        return -1;

      // fill dates array
      std::list<double> dates;
      for (; !step_sql ();)
        {
          double d = get_sql_real (0);
          dates.push_back (d);
        }
      finalize_sql ();

      // calc COMPDATs on every date
      std::list<double>::const_iterator di = dates.begin();
      for(ulong i = 0; i < dates.size(); ++i, ++di) {
          fci::compl_n_frac_builder cfb;
          cfb.init(nx, ny, trim_backend);
          cfb.init(sp_wp);
          if(i == 0) {
            //boost::timer::auto_cpu_timer t;
            //ProfilerStart("/home/uentity/my_projects/blue-sky.git/plugins/bs-eagle/examples/build_cache.prof");
            cfb.build_cache();
            //ProfilerStop();
          }
          else
            cfb.share_cache_with(cfb_storage.front());

          // find completions & fractures
          cfb.compl_build(*di);
          cfb_storage.push_back(cfb);
      }

      // wells in mesh come here
      std::set< std::string > good_wells;

      // check how many wells do we have
      prepare_sql("SELECT COUNT(*) FROM wells");
      if (!step_sql ())
        {
          BSOUT << "Wells count " << get_sql_int(0) << bs_end;
          //point->resize (3 * get_sql_int(0));
        }
      else {
        finalize_sql();
        return -1;
      }
      finalize_sql();


      /*-----------------------------------------------------------------
       * Main cycle - iterate over all dates
       *----------------------------------------------------------------*/
      cfb_iterator p_cfb = cfb_storage.begin();
      di = dates.begin();
      for (std::list<double>::const_iterator de = dates.end(); di != de; ++di, ++p_cfb)
        {
          char d_buf[1024];
          if (*di - int(*di) == 0) // if *di is integer
            {
              dt_t->d2ecl (*di, d_buf);
              printf ("DATE %s\n", d_buf);
              fprintf (fp, "DATES\n%s\n/\n\n", d_buf);
            }
          else
            {
              std::list<double>::const_iterator di_prev = di;
              di_prev--;
              double dt = *di - *di_prev;
              printf ("TSTEP %.10lf \t DATETIME %.10lf\n", dt, *di);
              fprintf (fp, "TSTEP\n%.10lf\n/\n\n", dt);
            }

          // WELSPECS only on first date
          if(di == dates.begin()) {
            fprintf (fp, "WELSPECS\n");

            if (prepare_sql ("SELECT name FROM wells ORDER BY name ASC"))
              return -1;

            for (int i = 0; !step_sql (); i++)
              {
                // obtain well name
                std::string s = get_sql_str (0);
                // we can really just skip wells without any COMPDATs
                // so we need to check compdats on every date
                for(cfb_citerator p_cfb = cfb_storage.begin(), cfb_end = cfb_storage.end(); p_cfb != cfb_end; ++p_cfb) {
                  const fci::cd_storage& cd = p_cfb->storage_compdat();
                  fci::cd_storage::const_iterator p_cd = fci::compdat::find_first_cd(cd, s);
                  if(p_cd != cd.end()) {
                    // well has at least one COMPDAT, so write it to WELLSPEC
                    // using first COMPDAT cell_id as well's position
                    BSOUT << "Exporting well " << s << bs_end;
                    fprintf (fp,
                        (boost::format("\'%s\' \'FIELD\' %u %u /\n") % s %
                        (p_cd->cell_pos[0] + 1) % (p_cd->cell_pos[1] + 1)).str().c_str());

                    // remember well as good -- not really needed now
                    good_wells.insert(s);
                    break;
                  }
                }
                // quick and dirty check if well is good enough
                if(good_wells.find(s) == good_wells.end()) {
                  BSOUT << "Warning! Well " << s << " has no active COMPDATs on any date and won't be exported!" << bs_end;
                }
              }
            fprintf (fp, "/\n\n");
            finalize_sql ();
          }

          // COMPDAT
          // Building compdat to put first completion I, J to WELLSPECS
          const fci::cd_storage& cd = p_cfb->storage_compdat();
          const fci::cd_storage::const_iterator cde = cd.end();
          const double eps = 1.0e-5;
          const std::set< std::string >::const_iterator good_wells_end = good_wells.end();
          if (!cd.empty ())
          {
            fprintf (fp, "COMPDAT\n");
            for (fci::cd_storage::const_iterator cdi = cd.begin(); cdi != cde; ++cdi)
            {
              // skip out of mesh wells
              if (fabs(cdi->kh_mult) > eps && good_wells.find(cdi->well_name) != good_wells_end)
                {
                  if (cdi->status)
                    fprintf (fp,
                      (boost::format("\'%s\' %u %u %u %u \'OPEN\' 2* %lf 1* %lf 1* \'%c\' /\n") %
                      cdi->well_name % (cdi->cell_pos[0] + 1) % (cdi->cell_pos[1] + 1) % (cdi->cell_pos[2] + 1) %
                      (cdi->cell_pos[3] + 1) % cdi->diam % cdi->skin % cdi->dir).str().c_str());
                  else
                    fprintf (fp,
                      (boost::format("\'%s\' %u %u %u %u \'SHUT\' 2* %lf 1* %lf 1* \'%c\' /\n") %
                      cdi->well_name % (cdi->cell_pos[0] + 1) % (cdi->cell_pos[1] + 1) %
                      (cdi->cell_pos[2] + 1) % (cdi->cell_pos[3] + 1) % cdi->diam % cdi->skin % cdi->dir).str().c_str());
                }
            }
            fprintf (fp, "/\n\n");
          }


          // WPIMULT

          if (!cd.empty ())
          {
            int wpimult_exist = 0;
            for (fci::cd_storage::const_iterator cdi = cd.begin(); cdi != cde; ++cdi)
            {
              if (fabs(cdi->kh_mult - 1.0) > 1e-6 && std::abs(cdi->kh_mult) > eps && good_wells.find(cdi->well_name) != good_wells_end)
                {
                  if (!wpimult_exist)
                    fprintf (fp, "WPIMULT\n");
                  wpimult_exist = 1;
                  fprintf (fp,
                      (boost::format("\'%s\' %lf %u %u %u /\n") % cdi->well_name % cdi->kh_mult % (cdi->cell_pos[0] + 1) %
                      (cdi->cell_pos[1] + 1) % (cdi->cell_pos[2] + 1)).str().c_str());
                }
            }
            if (wpimult_exist)
              fprintf (fp, "/\n\n");
          }


          // FRACTURES
          const fci::frac_storage& ft = p_cfb->frac_build(*di);
          const fci::frac_storage::const_iterator fte = ft.end();

          if (!ft.empty ())
          {
            char main_k_str[1024];
            char perm_str[1024];
            fprintf (fp, "FRACTURE\n");
            for (fci::frac_storage::const_iterator fti = ft.begin (); fti != fte; ++fti)
            {
              // skip out-of-mesh wells
              if(good_wells.find(fti->well_name) == good_wells_end)
                continue;

              if (fti->frac_perm > 0)
                sprintf (perm_str, "%lf", fti->frac_perm);
              else
                sprintf (perm_str, " * ");

              int horiz = 0;
              std::string sql = "SELECT name, horiz FROM wells WHERE name = ";
              sql += std::string("'") + fti->well_name + std::string("'");
              if (prepare_sql (sql.c_str ()))
                return -1;
              if (!step_sql ())
                {
                  horiz = get_sql_int (1);
                }
              finalize_sql ();

              if (horiz)
                sprintf (main_k_str, "%d", fti->frac_main_k + 1);
              else
                sprintf (main_k_str, " * ");

              fprintf (fp,
                (boost::format("\'%s\' %u %u %u %u %lf %lf %lf \'%s\' %lf %s %s %s %s /\n") %
                fti->well_name % (fti->cell_pos[0] + 1) % (fti->cell_pos[1] + 1) % (fti->cell_pos[2] + 1) %
                (fti->cell_pos[3] + 1) % fti->frac_half_length_1 % (fti->frac_angle - 90) % fti->frac_skin %
                (fti->frac_status == 0 ? "SHUT" : "OPEN") % fti->frac_half_thin % perm_str % " * " % " * " %
                main_k_str).str().c_str()
              );
            }
            fprintf (fp, "/\n\n");
          }

          // WCONINJH
          // Injector wells with rate control having historical bhp value
          sprintf (s_buf, "SELECT * FROM well_hist WHERE d=%lf AND ctrl < -1 AND i_bhp >= 0 ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;
          t_uint wconinjh_flag = 0;
          for (; !step_sql ();)
            {
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              if (wconinjh_flag == 0)
                {
                  fprintf (fp, "WCONINJH\n");
                  wconinjh_flag++;
                }

              int status = get_sql_int (14);
              int ctrl = get_sql_int (13);
              double i_or = get_sql_real (8);
              double i_wr = get_sql_real (9);
              double i_gr = get_sql_real (10);
              double i_bhp = get_sql_real (11);
              double rate = i_wr;
              std::string s_status;
              std::string s_phase = "WATER";
              std::string s_rate;
              std::string s_bhp;
              
              switch (status)
                {
                  case STATUS_OPEN:
                    s_status = "OPEN";
                    break;
                  case STATUS_SHUT:
                    s_status = "SHUT";
                    break;
                  case STATUS_CLOSE:
                    s_status = "CLOSE";
                    break;
                  default:
                    s_status = "SHUT";
                }

              if (ctrl == CTRL_I_ORATE)
                {
                  s_phase = "OIL";
                  rate = i_or;
                }
              else if (ctrl == CTRL_I_GRATE)
                {
                  s_phase = "GAS";
                  rate = i_gr;
                }
              
              if (rate > 0)
                { 
                  s_rate = V2S(rate);   
                }
              else
                {
                  s_rate = "*";   
                }

              if (i_bhp > 0)
                { 
                  s_bhp = V2S(i_bhp);   
                }
              else
                {
                  s_bhp = "*";   
                }

                
              fprintf (fp, "\'%s\' \'%s\' \'%s\' %s %s ", s.c_str (), s_phase.c_str(), s_status.c_str(), s_rate.c_str(), s_bhp.c_str());
              fprintf (fp, "/\n");
            }
          if (wconinjh_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();


          // WCONINJE
          sprintf (s_buf, "SELECT * FROM well_hist WHERE d=%lf AND (i_bhp < 0 AND ctrl < -1 OR ctrl = -1) ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;
          t_uint wconinje_flag = 0;
          for (; !step_sql ();)
            {
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              if (wconinje_flag == 0)
                {
                  fprintf (fp, "WCONINJE\n");
                  wconinje_flag++;
                }

              int status = get_sql_int (14);
              int ctrl = get_sql_int (13);
              double i_or = get_sql_real (8);
              double i_wr = get_sql_real (9);
              double i_gr = get_sql_real (10);
              double i_bhp = get_sql_real (11);
//              double lim_bhp = get_sql_real (11);
              double rate = i_wr;
              std::string s_status;
              std::string s_ctrl;
              std::string s_phase;
              std::string s_params;
              switch (status)
                {
                  case STATUS_OPEN:
                    s_status = "OPEN";
                    break;
                  case STATUS_SHUT:
                    s_status = "SHUT";
                    break;
                  case STATUS_CLOSE:
                    s_status = "CLOSE";
                    break;
                  default:
                    s_status = "SHUT";
                }
            
              s_ctrl = "BHP";
              // TODO: add injection phase into DB
              s_phase = "WATER";
              s_params = (boost::format("2* %lf") % i_bhp).str();
              if (ctrl == CTRL_I_BHP)
                {
                  s_ctrl = "BHP";
                  // TODO: add injection phase into DB
                  s_phase = "WATER";
                  s_params = (boost::format("2* %lf") % i_bhp).str();
                }
              else
                {
                  s_ctrl = "RATE";
                  if (ctrl == CTRL_I_WRATE)
                    s_phase = "WATER";
                  else if (ctrl == CTRL_I_ORATE)
                    {
                      s_phase = "OIL";
                      rate = i_or;
                    }
                  else if (ctrl == CTRL_I_GRATE)
                    {
                      s_phase = "GAS";
                      rate = i_gr;
                    }
                  else
                    s_phase = "WATER";

                  s_params = (boost::format("%s 2*") % V2S(rate)).str();
                }
              fprintf (fp, "\'%s\' \'%s\' \'%s\' \'%s\' %s ", s.c_str (), s_phase.c_str(), s_status.c_str(),
                  s_ctrl.c_str(), s_params.c_str());
              /*
              if (rate < 0)
                fprintf (fp, "2* ");
              else
                fprintf (fp, "%lf 1* ", rate);
              if (i_bhp < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", i_bhp);
              */
              fprintf (fp, "/\n");
            }
          if (wconinje_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();

          

          // WCONHIST
          sprintf (s_buf, "SELECT * FROM well_hist WHERE d=%lf AND p_bhp >= 0 AND ctrl > 1 ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;

          t_uint wconhist_flag = 0;
          for (; !step_sql ();)
            {
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              if (wconhist_flag == 0)
                {
                  fprintf (fp, "WCONHIST\n");
                  wconhist_flag++;
                }

              int status = get_sql_int (14);
              int ctrl = get_sql_int (13);
              double p_or = get_sql_real (2);
              double p_wr = get_sql_real (3);
              double p_gr = get_sql_real (4);
              double p_bhp = get_sql_real (6);
              double p_lr = get_sql_real (5);
              std::string s_status;
              std::string s_ctrl;
              std::string s_params;
              switch (status)
                {
                  case STATUS_OPEN:
                    s_status = "OPEN";
                    break;
                  case STATUS_SHUT:
                    s_status = "SHUT";
                    break;
                  case STATUS_CLOSE:
                    s_status = "CLOSE";
                    break;
                  default:
                    s_status = "SHUT";
                }
              if (ctrl == CTRL_P_LRATE)
                {
                  s_ctrl = "LRAT";
                  if (p_or + p_wr == 0)
                    p_or = p_lr;
                  if (status == STATUS_OPEN && (p_or + p_wr == 0))
                    continue;
                }
              else if (ctrl == CTRL_P_ORATE)
                {
                  s_ctrl = "ORAT";
                  if (status == STATUS_OPEN && (p_or == 0))
                    continue;
                }
              else if (ctrl == CTRL_P_WRATE)
                {
                  s_ctrl = "WRAT";
                  if (status == STATUS_OPEN && (p_wr == 0))
                    continue;
                }
              else if (ctrl == CTRL_P_GRATE)
                {
                  s_ctrl = "GRAT";
                  if (status == STATUS_OPEN && (p_gr == 0))
                    continue;
                }
             
              if (p_bhp > 0)
                { 
                  s_params = (boost::format("%s %s %s 3* %s") % V2S(p_or) % V2S(p_wr)% V2S(p_gr) % V2S(p_bhp)).str();
                }
              else
                {
                  s_params = (boost::format("%s %s %s 4*") % V2S(p_or) % V2S(p_wr)% V2S(p_gr)).str();
                }
              
              fprintf (fp, "\'%s\' \'%s\' \'%s\' %s ", s.c_str (), s_status.c_str(), s_ctrl.c_str(), s_params.c_str());
              /*
              if (p_or < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_or);
              if (p_wr < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_wr);
              if (p_gr < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_gr);
              if (p_lr < 0)
                fprintf (fp, "2* ");
              else
                fprintf (fp, "%lf * ", p_lr);
              if (p_bhp < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_bhp);
              */
              fprintf (fp, "/\n");
            }


          if (wconhist_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();


          // WCONPROD
          sprintf (s_buf, "SELECT * FROM well_hist WHERE d=%lf AND (p_bhp < 0 AND ctrl > 1 OR ctrl = 1) ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;

          t_uint wconprod_flag = 0;
          for (; !step_sql ();)
            {
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              if (wconprod_flag == 0)
                {
                  fprintf (fp, "WCONPROD\n");
                  wconprod_flag++;
                }

              int status = get_sql_int (14);
              int ctrl = get_sql_int (13);
              double p_or = get_sql_real (2);
              double p_wr = get_sql_real (3);
              double p_gr = get_sql_real (4);
              double p_lr = get_sql_real (5);
              double p_bhp = get_sql_real (6);
              std::string s_status;
              std::string s_ctrl;
              std::string s_params;
              switch (status)
                {
                  case STATUS_OPEN:
                    s_status = "OPEN";
                    break;
                  case STATUS_SHUT:
                    s_status = "SHUT";
                    break;
                  case STATUS_CLOSE:
                    s_status = "CLOSE";
                    break;
                  default:
                    s_status = "SHUT";
                }
              if (ctrl == CTRL_P_LRATE)
                {
                  s_ctrl = "LRAT";
                  s_params = (boost::format("3* %s 2*") % V2S(p_lr)).str();
                }
              else if (ctrl == CTRL_P_BHP)
                {
                  s_ctrl = "BHP";
                  s_params = (boost::format("5* %s") % V2S(p_bhp)).str();
                }
              else if (ctrl == CTRL_P_ORATE)
                {
                  s_ctrl = "ORAT";
                  s_params = (boost::format("%s 5*") % V2S(p_or)).str();
                }
              else if (ctrl == CTRL_P_WRATE)
                {
                  s_ctrl = "WRAT";
                  s_params = (boost::format("1* %s 4*") % V2S(p_wr)).str();
                }
              else if (ctrl == CTRL_P_GRATE)
                {
                  s_ctrl = "GRAT";
                  s_params = (boost::format("2* %s 3*") % V2S(p_gr)).str();
                }

              fprintf (fp, "\'%s\' \'%s\' \'%s\' %s ", s.c_str (), s_status.c_str(), s_ctrl.c_str(), s_params.c_str());
              /*
              if (p_or < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_or);
              if (p_wr < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_wr);
              if (p_gr < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_gr);
              if (p_lr < 0)
                fprintf (fp, "2* ");
              else
                fprintf (fp, "%lf * ", p_lr);
              if (p_bhp < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_bhp);
              */
              fprintf (fp, "/\n");
            }


          if (wconprod_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();

          // WELTARG BHP limit for injectors
          sprintf (s_buf, "SELECT well_name, lim_i_bhp, lim_p_bhp FROM well_hist WHERE d=%lf AND (lim_i_bhp > 0 OR lim_p_bhp > 0) ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;
          t_uint weltarg_flag = 0;
          for (; !step_sql ();)
            {
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              if (weltarg_flag == 0)
                {
                  fprintf (fp, "WELTARG\n");
                  weltarg_flag++;
                }

              double lim_i_bhp = get_sql_real (1);
              double lim_p_bhp = get_sql_real (2);

              if (lim_i_bhp > 0)
                fprintf (fp, "\'%s\' BHP %s ", s.c_str (), V2S(lim_i_bhp));
              else if (lim_p_bhp > 0)
                fprintf (fp, "\'%s\' BHP %s ", s.c_str (), V2S(lim_p_bhp));

             
              fprintf (fp, "/\n");
            }
          if (weltarg_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();

          // WEFAC
          sprintf (s_buf, "SELECT well_name, wefac FROM well_hist WHERE d=%lf ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;

          t_uint wefac_flag = 0;
          for (; !step_sql ();)
            {
              double wefac = get_sql_real (1);
              //if (wefac == 1.0)
              //  continue;

              std::string s = get_sql_str (0);
              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              if (wefac_flag == 0)
                {
                  fprintf (fp, "WEFAC\n");
                  wefac_flag++;
                }

              fprintf (fp, "\'%s\' %lf", s.c_str (), wefac);
              fprintf (fp, "/\n");
            }


          if (wefac_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();
        }   // eof main cycle over dates

      fclose (fp);
      return 0;
    }

  int
  sql_well::read_from_ascii_file (const std::string &fname, double starting_date)
    {
      //if (fr_file)
      //  delete fr_file;
      //fr_file = new FRead (fname.c_str (), fname.c_str ());
      fr_file = BS_KERNEL.create_object ("bos_reader");
      fr_file->open (fname.c_str (), fname.c_str ());

      int rc = 0;
      char buf[4096];
      char *id = 0;
      char *other = 0;
      double d = 0;

      for (;;)
        {
          rc = fr_file->read_line (buf, 4096, FREAD_DONT_CONVERT_CASE);
          //printf ("%s\n", buf);
          if (rc <= 0)
            break;
          d = starting_date;
          // read date and time
          rc = read_date_and_time (buf, &id, &d);
          if (rc)
            return -4;
          // read id
          fr_file->trim_left (&id);
          other = id;
          rc |= fr_file->get_phrase (&other);
          if (rc)
            return -1;
          //trim_right_s (id);
          fr_file->locale_ucase (id);

          if (!strcmp (id, "W_SPEC"))
            {
              rc = read_w_spec (other);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_BRANCH_F"))
            {
              rc = read_w_branch_f (other);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_COMP"))
            {
              rc = read_w_comp (other, d);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_FRAC"))
            {
              rc = read_w_frac (other, d);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_PROD"))
            {
              rc = read_w_prod (other, d);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_INJ"))
            {
              rc = read_w_inj (other, d);
              if (rc)
                return rc;
            }
          //printf ("date %lf\n", d);
          //printf ("next: %s\n", next);
        }
      return 0;
    }

  int
  sql_well::read_w_branch_f (char *buf)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char branch[1024];
      char parent[1024];
      char fname[1024];
      char fname2[1024];

      char sql[1024];
      double md = -1;
      // read well name
      wname[0] = '\0';
      nx = buf;
      // read well name
      rc = fr_file->get_phrase_str (&nx, wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_BRANCH_F can not be set by default\n");
          return -1;
        }
      // read branch
      strcpy (branch, "main");
      rc = fr_file->get_phrase_str (&nx, branch);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read branch name\n");
          return -1;
        }
      // read parent
      parent[0] = '\0';
      rc = fr_file->get_phrase_str (&nx, parent);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read parent name\n");
          return -1;
        }
      //read md
      rc = fr_file->get_phrase_double (&nx, &md);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F\n");
          return -1;
        }
      // read file name
      fname[0] = '\0';
      rc = fr_file->get_phrase_filepath (&nx, fname);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read file name\n");
          return -1;
        }
      // read well_log file name
      fname2[0] = '\0';
      rc = fr_file->get_phrase_filepath (&nx, fname2);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read file name\n");
          return -1;
        }
      // add to data base
      sprintf (sql, "INSERT OR REPLACE INTO branches(well_name, branch_name, md, parent) VALUES ('%s', '%s', %lf, '%s')",
               wname, branch, md, parent);
      printf ("SQL: %s\n", sql);
      if (exec_sql (sql))
        return -1;
      if (fname[0] != '\0')
        {
          sp_traj_t sp_traj = BS_KERNEL.create_object ("traj");
          rc = sp_traj->read_from_dev_file (fr_file->get_incdir() + std::string(fname));
          if (rc)
            {
              return rc;
            }
          if (add_branch_traj (wname, branch, sp_traj))
            return -6;
          printf ("TRAJ\n %s\n", sp_traj->py_str ().c_str ());
        }
      if (fname2[0] != '\0')
        {
          sp_gis_t sp_gis = BS_KERNEL.create_object ("gis");
          rc = sp_gis->read_from_las_file (fr_file->get_incdir() + std::string(fname2));
          if (rc)
            {
              return rc;
            }
          if (add_branch_gis (wname, branch, sp_gis))
            return -6;
        }

      printf ("W_BRANCH_F %s %s %s %lf\n", wname, branch, parent, md);

      return 0;
    }
  int
  sql_well::read_w_spec (char *buf)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char sql[1024];
      double x = -1, y = -1;
      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_SPEC can not be set by default\n");
          return -1;
        }
      rc = fr_file->get_phrase_double (&nx, &x);
      //printf ("rc: %lf\n", x);
      rc |= fr_file->get_phrase_double (&nx, &y);
      //printf ("rc: %lf\n", y);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      printf ("W_SPEC %s %lf %lf\n", wname, x, y);
      // add to data base
      sprintf (sql, "INSERT INTO wells(name, x, y) VALUES ('%s', %lf, %lf)",
               wname, x, y);
      return exec_sql (sql);
    }
  int
  sql_well::read_w_frac (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char branch[1024];
      char status[1024];
      int i_status;
      char sql[1024];
      double md = -1;
      double angle  = 0;
      double half_length_1 = 50.0;
      double half_length_2 = 50.0;
      double half_up = 5.0;
      double half_down = 5.0;
      double perm = -1;
      double half_thin = 0.005;
      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_FRAC can not be set by default\n");
          return -1;
        }
      strcpy (branch, "main");
      rc = fr_file->get_phrase_str (&nx, branch);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well branch name in keyword W_FRAC\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_FRAC\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &md);
      rc |= fr_file->get_phrase_double (&nx, &angle);
      rc |= fr_file->get_phrase_double (&nx, &half_length_1);
      rc |= fr_file->get_phrase_double (&nx, &half_length_2);
      rc |= fr_file->get_phrase_double (&nx, &half_up);
      rc |= fr_file->get_phrase_double (&nx, &half_down);
      rc |= fr_file->get_phrase_double (&nx, &perm);
      rc |= fr_file->get_phrase_double (&nx, &half_thin);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      if (status[0] == 'S')
        i_status = STATUS_CON_SHUT;
      //else if (status[0] == 'C') // close
      //  i_status = 1;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_CON_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_FRAC not aloowed\n", status);
          return -9;
        }
      if (md < 0)
        {
          fprintf (stderr, "Error: you should specify md for W_FRAC \n");
          return -9;
        }
      if (half_length_1 <= 0)
        {
          fprintf (stderr, "Error: HALF_LENGTH_1 should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_length_2 <= 0)
        {
          fprintf (stderr, "Error: HALF_LENGTH_2 should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_up <= 0)
        {
          fprintf (stderr, "Error: HALF_UP should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_down <= 0)
        {
          fprintf (stderr, "Error: HALF_DOWN should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_thin <= 0)
        {
          fprintf (stderr, "Error: HALF_THIN should be > 0 for W_FRAC \n");
          return -9;
        }

      printf ("W_FRAC %s %s %s %lf %lf %lf %lf %lf %lf %lf %lf\n",
              wname, branch, status, md, angle, half_length_1, half_length_2,
              half_up, half_down, perm, half_thin);
      // add to data base
      sprintf (sql, "INSERT INTO fractures(well_name, branch_name, md, d, status, \
half_up, half_down, angle, half_length_1, half_length_2, perm, half_thin) \
VALUES ('%s', '%s', %lf, %lf, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf)",
               wname, branch, md, d, i_status, half_up, half_down, angle,
               half_length_1, half_length_2, perm, half_thin);
      return exec_sql (sql);
      return 0;
    }
  int
  sql_well::read_w_comp (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char branch[1024];
      char status[1024];
      int i_status;
      char sql[1024];
      double md = -1, length  = 1, rw = 0.08, skin = 0, khmult = 1.0;
      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_COMP can not be set by default\n");
          return -1;
        }
      strcpy (branch, "main");
      rc = fr_file->get_phrase_str (&nx, branch);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well branch name in keyword W_COMP\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_COMP\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &md);
      rc |= fr_file->get_phrase_double (&nx, &length);
      rc |= fr_file->get_phrase_double (&nx, &rw);
      rc |= fr_file->get_phrase_double (&nx, &skin);
      rc |= fr_file->get_phrase_double (&nx, &khmult);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      if (status[0] == 'S')
        i_status = STATUS_CON_SHUT;
      //else if (status[0] == 'C') // close
      //  i_status = 1;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_CON_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_COMP not aloowed\n", status);
          return -9;
        }
      if (md < 0)
        {
          fprintf (stderr, "Error: you should specify md for W_COMP \n");
          return -9;
        }
      if (length <= 0)
        {
          fprintf (stderr, "Error: length should be > 0 for W_COMP \n");
          return -9;
        }
      if (rw <= 0)
        {
          fprintf (stderr, "Error: rw should be > 0 for W_COMP \n");
          return -9;
        }
      if (khmult <= 0)
        {
          fprintf (stderr, "Error: khmult should be > 0 for W_COMP \n");
          return -9;
        }

      printf ("W_COMP %s %s %lf %lf %lf %lf %lf\n",
              wname, branch, md, length, rw, skin, khmult);
      // add to data base
      sprintf (sql, "INSERT INTO completions(well_name, branch_name, md, d, length, status, rw, kh_mult) VALUES ('%s', '%s', %lf, %lf, %lf, %d, %lf, %lf)",
               wname, branch, md, d, length, i_status, rw, khmult);
      return exec_sql (sql);
      return 0;
    }

  int
  sql_well::read_w_inj (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char status[1024];
      char ctrl[1024];
      char fluid[1024];
      int i_status = 0;
      int i_ctrl = 0;
      char sql[1024];
      double bhp = -1;
      double rate = -1;
      double orate = -1;
      double wrate = -1;
      double grate = -1;
      double lim_bhp = -1;
      double lim_rate = -1;
      double lim_orate = -1;
      double lim_wrate = -1;
      double lim_grate = -1;
      double wefac = 1.0;

      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_INJ can not be set by default\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_INJ\n");
          return -1;
        }
      strcpy (ctrl, "BHP");
      rc = fr_file->get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_INJ\n");
          return -1;
        }
      strcpy (fluid, "WATER");
      rc = fr_file->get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_INJ\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &bhp);
      rc |= fr_file->get_phrase_double (&nx, &rate);
      rc = fr_file->get_phrase_double  (&nx, &lim_bhp);
      rc |= fr_file->get_phrase_double (&nx, &lim_rate);
      rc |= fr_file->get_phrase_double (&nx, &wefac);
      if (rc)
        {
          fprintf (stderr, "Error: W_PROD\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      fr_file->locale_ucase (ctrl);
      fr_file->locale_ucase (fluid);
      if (status[0] == 'S')
        i_status = STATUS_SHUT;
      else if (status[0] == 'C') // close
        i_status = STATUS_CLOSE;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_COMP not aloowed\n", status);
          return -9;
        }
      if (ctrl[0] == 'B')
        {
          i_ctrl = CTRL_I_BHP;
          if (i_status == 2 && bhp <= 0)
            {
              fprintf (stderr, "Error: BHP = %lf for CONTROL %s in keyword W_INJ",
                       bhp, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'R') // rate
        {
          if (i_status == STATUS_OPEN && rate <= 0)
            {
              fprintf (stderr, "Error: RATE = %lf for CONTROL %s in keyword W_INJ",
                       wrate, ctrl);
              return -1;
            }
          if (fluid[0] == 'W')
            {
              i_ctrl = CTRL_I_WRATE;
              wrate = rate;
              lim_wrate = lim_rate;
            }
          else if (fluid[0] == 'O')
            {
              i_ctrl = -CTRL_I_ORATE;
              orate = rate;
              lim_orate = lim_rate;
            }
          else if (fluid[0] == 'G')
            {
              i_ctrl = CTRL_I_GRATE;
              grate = rate;
              lim_grate = lim_rate;
            }
          else
            {
              fprintf (stderr, "Error: FLUID = %s not allowed in keyword W_INJ",
                       fluid);
              return -1;
            }
        }
      else
        {
          fprintf (stderr, "Error: control %s for W_INJ not aloowed\n", ctrl);
          return -9;
        }
      if (i_status == STATUS_OPEN && wefac <= 0)
        {
          fprintf (stderr, "Error: WEFAC = %lf in keyword W_INJ",
                   wefac);
          return -1;
        }

      printf ("W_INJ %s %s %s %s %lf %lf %lf %lf %lf\n",
              wname, status, ctrl, fluid, bhp, rate, lim_bhp, lim_rate, wefac);
      // add to data base
      sprintf (sql, "INSERT INTO well_hist(well_name, d, i_or, i_wr, i_gr, \
i_bhp, wefac, ctrl, status, lim_i_or, lim_i_wr, lim_i_gr, lim_i_bhp) \
VALUES ('%s', %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf, %lf)",
               wname, d, orate, wrate, grate, bhp, wefac, i_ctrl, i_status,
               lim_orate, lim_wrate, lim_grate, lim_bhp);
      return exec_sql (sql);
      return 0;
    }

  int
  sql_well::read_w_prod (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char status[1024];
      char ctrl[1024];
      int i_status = 0;
      int i_ctrl = 0;
      char sql[1024];
      double bhp = -1;
      double orate = -1;
      double wrate = -1;
      double grate = -1;
      double lrate = -1;
      double lim_bhp = -1;
      double lim_orate = -1;
      double lim_wrate = -1;
      double lim_grate = -1;
      double lim_lrate = -1;
      double wefac = 1.0;

      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_PROD can not be set by default\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_PROD\n");
          return -1;
        }
      strcpy (ctrl, "BHP");
      rc = fr_file->get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_PROD\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &bhp);
      rc |= fr_file->get_phrase_double (&nx, &wrate);
      rc |= fr_file->get_phrase_double (&nx, &orate);
      rc |= fr_file->get_phrase_double (&nx, &grate);
      rc |= fr_file->get_phrase_double (&nx, &lrate);
      rc = fr_file->get_phrase_double  (&nx, &lim_bhp);
      rc |= fr_file->get_phrase_double (&nx, &lim_wrate);
      rc |= fr_file->get_phrase_double (&nx, &lim_orate);
      rc |= fr_file->get_phrase_double (&nx, &lim_grate);
      rc |= fr_file->get_phrase_double (&nx, &lim_lrate);

      rc |= fr_file->get_phrase_double (&nx, &wefac);
      if (rc)
        {
          fprintf (stderr, "Error: W_PROD\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      fr_file->locale_ucase (ctrl);
      if (status[0] == 'S')
        i_status = STATUS_SHUT;
      else if (status[0] == 'C') // close
        i_status = STATUS_CLOSE;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_COMP not aloowed\n", status);
          return -9;
        }
      if (ctrl[0] == 'B')
        {
          i_ctrl = CTRL_P_BHP;
          if (i_status == 2 && bhp <= 0)
            {
              fprintf (stderr, "Error: BHP = %lf for CONTROL %s in keyword W_PROD",
                       bhp, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'W') // close
        {
          i_ctrl = CTRL_P_WRATE;
          if (i_status == 2 && wrate <= 0)
            {
              fprintf (stderr, "Error: WRATE = %lf for CONTROL %s in keyword W_PROD",
                       wrate, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'O') // close
        {
          i_ctrl = CTRL_P_ORATE;
          if (i_status == 2 && orate <= 0)
            {
              fprintf (stderr, "Error: ORATE = %lf for CONTROL %s in keyword W_PROD",
                       orate, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'G') // close
        {
          i_ctrl = CTRL_P_GRATE;
          if (i_status == 2 && grate <= 0)
            {
              fprintf (stderr, "Error: GRATE = %lf for CONTROL %s in keyword W_PROD",
                       grate, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'L') // close
        {
          i_ctrl = CTRL_P_LRATE;
          if (i_status == 2 && lrate <= 0)
            {
              fprintf (stderr, "Error: LRATE = %lf for CONTROL %s in keyword W_PROD",
                       lrate, ctrl);
              return -1;
            }
        }
      else
        {
          fprintf (stderr, "Error: control %s for W_PROD not aloowed\n", status);
          return -9;
        }
      if (i_status == STATUS_OPEN && wefac <= 0)
        {
          fprintf (stderr, "Error: WEFAC = %lf in keyword W_PROD",
                   wefac);
          return -1;
        }

      printf ("W_PROD %s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
              wname, status, ctrl, bhp, wrate, orate, grate, lrate, lim_bhp,
              lim_wrate, lim_orate, lim_grate, lim_lrate, wefac);
      // add to data base
      sprintf (sql, "INSERT INTO well_hist(well_name, d, p_or, p_wr, p_gr, p_lr, \
p_bhp, wefac, ctrl, status, lim_p_or, lim_p_wr, lim_p_gr, lim_p_lr, lim_p_bhp) \
VALUES ('%s', %lf, %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf, %lf, %lf)",
               wname, d, orate, wrate, grate, lrate, bhp, wefac, i_ctrl, i_status,
               lim_orate, lim_wrate, lim_grate, lim_lrate, lim_bhp);
      return exec_sql (sql);
      return 0;
    }

  int
  sql_well::read_date_and_time (char *buf, char **next_start, double *dd)
    {
      int rc = 0;
      double days;
      double t;
      char *nx;

      *next_start = buf;
      fr_file->trim_left (next_start);
      nx = *next_start;
      rc |= fr_file->get_phrase (&nx);
      if (**next_start == '*')
        days = 0;
      else
        {
          rc |= fr_file->get_dt ()->cstr2d (*next_start, days);
          *dd = (double)days;
        }

      *next_start = nx;
      fr_file->trim_left (next_start);
      nx = *next_start;
      rc |= fr_file->get_phrase (&nx);
      if (**next_start == '*')
        t = 0;
      else
        {
          rc |= fr_file->get_dt ()->cstr2t (*next_start, t);
        }
      *next_start = nx;

      *dd += (double)t;
      return rc;
    }

} /* namespace blue_sky */
