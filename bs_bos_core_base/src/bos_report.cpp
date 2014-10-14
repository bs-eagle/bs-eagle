#include "bos_report.h"
#include "bs_assert.h"
#include "bs_log_scribers.h"

using namespace std;

namespace blue_sky
  {
  void init_bos_logs()
  {
    //string log_file;
    //char *c_dir = NULL;
    //if (!(c_dir = getenv("BS_KERNEL_DIR")))
    //  c_dir = (char *)".";

    STLOG.add_channel(sp_channel(new bs_channel("bos_out")));
    STLOG.add_channel(sp_channel(new bs_channel("bos_warn")));
    STLOG.add_channel(sp_channel(new bs_channel("bos_err")));

    BOSOUT.get_channel ()->attach(sp_stream(new log::detail::cout_scriber ("COUT")));
    BOSOUT.get_channel ()->attach(sp_stream(new log::detail::file_scriber ("FILE", "./bos_main_out.log", ios::out|ios::app)));

    // TODO:
    // now (at 07.07.2009) we can add output_time to log,
    // but we don't do it
    //BOSOUT << output_time;
    BOSOUT (section::app_info, level::debug) << "bos_out log has been registered!" << bs_end;



    BOSWARN.get_channel ()->set_prefix("Warning: ");
    BOSWARN.get_channel ()->attach(sp_stream(new log::detail::cout_scriber ("COUT")));
    BOSWARN.get_channel ()->attach(sp_stream(new log::detail::file_scriber ("FILE", "./bos_main_warn.log", ios::out|ios::app)));
    BOSWARN << output_time;
    BOSOUT (section::app_info, level::debug) << "bos_warn log has been registered!" << bs_end;


    BOSERR.get_channel ()->set_prefix("Error: ");
    BOSERR.get_channel ()->attach(sp_stream(new log::detail::cout_scriber ("COUT")));
    BOSERR.get_channel ()->attach(sp_stream(new log::detail::file_scriber ("FILE", "./bos_main_err.log", ios::out|ios::app)));
    BOSERR << output_time;
    BOSOUT (section::app_info, level::debug) << "bos_err log has been registered!" << bs_end;

    log::detail::add_section_to_channel (*BOSOUT.get_channel ());
    log::detail::add_section_to_channel (*BOSWARN.get_channel ());
    log::detail::add_section_to_channel (*BOSERR.get_channel ());

    //! use MTLOG.add_log_channel(sp_channel(new bs_channel(BOSOUT)));
    //! to copy existing channel in new thread
  }
}
