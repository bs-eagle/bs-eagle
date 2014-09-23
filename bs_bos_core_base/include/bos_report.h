#ifndef BOS_REPORT_H
#define BOS_REPORT_H

#include "bs_report.h"
#include "bs_kernel.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#endif

#define STLOG BS_KERNEL.get_log ()
#define MTLOG BS_KERNEL.get_tlog ()

#define BOSOUT STLOG.get_locked ("bos_out", __FILE__, __LINE__)
#define BOSERR STLOG.get_locked ("bos_err", __FILE__, __LINE__)
#define BOSWARN STLOG.get_locked ("bos_warn", __FILE__, __LINE__)

#define MBOSOUT MTLOG["bos_out"]
#define MBOSERR MTLOG["bos_err"]
#define MBOSWARN MTLOG["bos_warn"]

#define MTLOG_INIT \
	MTLOG.add_log_channel(sp_channel(new bs_channel(*BOSOUT))); \
	MTLOG.add_log_channel(sp_channel(new bs_channel(*BOSERR))); \
	MTLOG.add_log_channel(sp_channel(new bs_channel(*BOSWARN))); \

#define MTLOG_FINI \
	MTLOG.rem_log_channel("bos_out"); \
	MTLOG.rem_log_channel("bos_err"); \
	MTLOG.rem_log_channel("bos_warn"); \

namespace blue_sky
{
  //namespace sn
  //  {
  //  typedef enum
  //  {
  //    sn_start,
  //    read,        /*Messages about reading input data*/
  //    write,       /*Messages about writing input data*/
  //    check,       /*Messages about checking input data*/
  //    init,        /*Messages about initialization of simulation process*/
  //    iters,       /*Messages about iteration process*/
  //    fip,         /*Messages about fip*/
  //    solve,       /*Messages about solving*/
  //    gui,         /*Messages about user interface*/
  //    jacobian,    /*Messages about jacobian*/
  //    rsv,         /*Messages about reservoir*/
  //    group,       /*Messages about well groups*/
  //    well,        /*Messages about wells*/
  //    conn,        /*Messages about connections*/
  //    print_logs,  /*Messages about logs*/
  //    sys_info,    /*Messages about current system*/
  //    thread_info, /*Messages about current thread*/
  //    sn_end
  //  } section_names;
  //}

#define DECL_LOG_SECTIONS_ENUM_I(r, data, i, elem)                          \
  BOOST_PP_TUPLE_ELEM (3, 0, elem),

#define DECL_LOG_LEVELS_ENUM_I(r, data, i, elem)                            \
  BOOST_PP_TUPLE_ELEM (3, 0, elem),

#define DECL_LOG_SECTIONS_KEYWORD_HANDLER_I(r, data, i, elem)               \
  if (BOOST_PP_TUPLE_ELEM (3, 1, elem) == name)                             \
    {                                                                       \
      using namespace section;                                              \
      BOSOUT.set_priority (BOOST_PP_TUPLE_ELEM (3, 0, elem), priority);     \
      return true;                                                          \
    }                                                                       \
  else 

#define DECL_LOG_SECTIONS_SET_PRIORITY_STREAM_I(r, data, i, elem)           \
  if (BOOST_PP_TUPLE_ELEM (3, 1, elem) == name)                             \
    {                                                                       \
      using namespace section;                                              \
      stream->set_priority (priority (BOOST_PP_TUPLE_ELEM (3, 0, elem), level));    \
      return true;                                                          \
    }                                                                       \
  else

#define DECL_LOG_LEVELS_FUNCTION_I(r, data, i, elem)                        \
  if (BOOST_PP_TUPLE_ELEM (3, 1, elem) == name)                             \
    {                                                                       \
      using namespace blue_sky::level;                                      \
      return BOOST_PP_TUPLE_ELEM (3, 0, elem);                              \
    }                                                                       \
  else  

#define DECL_LOG_SECTIONS_ADD_SECTION_I(r, data, i, elem)                   \
  .add_section (BOOST_PP_TUPLE_ELEM (3, 0, elem),                           \
    BOOST_PP_TUPLE_ELEM (3, 2, elem))

#define DECLARE_LOG_SECTIONS_ENUM(section_name, seq)                        \
  namespace section_name                                                    \
  {                                                                         \
    enum                                                                    \
    {                                                                       \
      section_begin,                                                        \
      BOOST_PP_SEQ_FOR_EACH_I (DECL_LOG_SECTIONS_ENUM_I, _, seq)            \
      section_end,                                                          \
    };                                                                      \
  }

#define DECLARE_LOG_SECTIONS_KEYWORD_HANDLER(function_name, seq)            \
  namespace log { namespace detail {                                        \
    inline bool                                                             \
    function_name (const std::string &name, int priority)                   \
    {                                                                       \
      BOOST_PP_SEQ_FOR_EACH_I (DECL_LOG_SECTIONS_KEYWORD_HANDLER_I, _, seq) \
        {                                                                   \
          return false;                                                     \
        }                                                                   \
    }                                                                       \
    inline bool                                                                 \
    function_name (const std::string &name, int level, const sp_stream &stream) \
    {                                                                           \
      BOOST_PP_SEQ_FOR_EACH_I (DECL_LOG_SECTIONS_SET_PRIORITY_STREAM_I, _, seq) \
        {                                                                   \
          return false;                                                     \
        }                                                                   \
    }                                                                       \
  } }

#define DECLARE_LOG_LEVELS_ENUM(name, seq)                                  \
  namespace name                                                            \
  {                                                                         \
    enum                                                                    \
    {                                                                       \
      BOOST_PP_SEQ_FOR_EACH_I (DECL_LOG_LEVELS_ENUM_I, _, seq)              \
    };                                                                      \
  }

#define DECLARE_LOG_LEVELS_FUNCTION(function_name, seq)                     \
  namespace log { namespace detail {                                        \
    inline int                                                              \
    function_name (const std::string &name)                                 \
    {                                                                       \
      BOOST_PP_SEQ_FOR_EACH_I (DECL_LOG_LEVELS_FUNCTION_I, _, seq)          \
        {                                                                   \
          return -1;                                                        \
        }                                                                   \
    }                                                                       \
  } }

#define DECLARE_LOG_SECTIONS_ADD_SECTION(seq)                               \
  namespace log { namespace detail {                                        \
    inline void                                                             \
    add_section_to_channel (bs_channel &ch)                                 \
    {                                                                       \
      using namespace section;                                              \
      ch                                                                    \
      BOOST_PP_SEQ_FOR_EACH_I (DECL_LOG_SECTIONS_ADD_SECTION_I, _, seq)     \
      ;                                                                     \
    }                                                                       \
  } }

#define DECLARE_LOG_SECTIONS(section_name, function_name, sections_seq)     \
  DECLARE_LOG_SECTIONS_ENUM (section_name, sections_seq)                    \
  DECLARE_LOG_SECTIONS_KEYWORD_HANDLER (function_name, sections_seq)        \
  DECLARE_LOG_SECTIONS_ADD_SECTION (sections_seq)

#define DECLARE_LOG_LEVELS(level_name, function_name, levels_seq)           \
  DECLARE_LOG_LEVELS_ENUM (level_name, levels_seq)                          \
  DECLARE_LOG_LEVELS_FUNCTION (function_name, levels_seq)

  DECLARE_LOG_LEVELS (level, get_log_level_by_name,
      ((critical,    "CRITICAL",  ""))  /* Critical errors such as "Cannot allocate memory" */ 
      ((error,       "ERROR",     ""))  /* Program errors occurs because of input data */
      ((warning,     "WARNING",   ""))  /* Program warnings */
      ((highest,     "HIGHEST",   ""))  /* Highest message priority level */
      ((high,        "HIGH",      ""))  /* High message priority level */
      ((medium,      "MEDIUM",    ""))  /* Medium message priority level */
      ((low,         "LOW",       ""))  /* Low message priority level */
      ((debug,       "DEBUG",     ""))  /* Developing messages not for end users */
    );

  DECLARE_LOG_SECTIONS (section, set_section_level_by_name,
      ((iters,       "ITERS",       level::low))      /* Messages about iteration process */                           
      ((wells,       "WELLS",       level::low))      /* Messages about wells */
      ((mesh,        "MESH",        level::low))      /* Messages about mesh */
      ((solvers,     "SOLVERS",     level::low))      /* Messages about linear solvers except AMG */
      ((amg,         "AMG",         level::low))      /* Messages about AMG solver */
      ((read_data,   "READ_DATA",   level::low))      /* Messages about reading input data */
      ((init_data,   "INIT_DATA",   level::low))      /* Messages about initialization data */
      ((check_data,  "CHECK_DATA",  level::low))      /* Messages about checking data */
      ((save_data,   "SAVE_DATA",   level::low))      /* Messages about save data */
      ((main_loop,   "MAIN_LOOP",   level::low))      /* Messages about main_loop */
      ((pvt,         "PVT",         level::low))      /* Messages about pvt */
      ((scal,        "SCAL",        level::low))      /* Messages about scal */
      ((schedule,    "SCHEDULE",    level::low))      /* Messages about schedule */
      ((keywords,    "KEYWORDS",    level::low))      /* Messages about keywords */
      ((fip,         "FIP",         level::low))      /* Messages about fip */
      ((arithmetic,  "ARITHMETIC",  level::low))      /* Messages about arithmetic */
      ((print_reports, "PRINT_REPORTS", level::low))  /* Report messages */
      ((app_info,    "APP_INFO",    level::debug))    /* Messages about internal application state (for debug only) */
      ((sys_info,    "SYS_INFO",    level::debug))    /* Messages about external system info */
      ((h5,          "H5",          level::low))      /* H5 specific messages */
    );

  BS_API_PLUGIN 
  void init_bos_logs();

} // namespace blue_sky

#endif // BOS_REPORT_H

