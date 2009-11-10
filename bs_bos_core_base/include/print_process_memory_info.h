/**
 * \file print_process_memory_info.h
 * \brief helper to output process memory info (for debug purposes)
 * \author Sergey Miryanov
 * \date 22.01.2009
 * */
#ifndef BS_TOOLS_PRINT_PROCESS_MEMORY_INFO_H_
#define BS_TOOLS_PRINT_PROCESS_MEMORY_INFO_H_


#ifndef UNIX
#ifdef _DEBUG
#include <psapi.h>
#pragma comment (lib, "psapi.lib")
#endif
#endif //!UNIX

namespace blue_sky {
namespace debug {

  inline void
  print_memory_info (const std::string &
#ifndef UNIX
#ifdef _DEBUG
      context_name
#endif //_DEBUG
#endif //UNIX
      )
  {
#ifndef UNIX
#ifdef _DEBUG

    PROCESS_MEMORY_COUNTERS pcmex = {0};
    HANDLE process = GetCurrentProcess ();
    if (process == NULL)
      {
        BOSERR (section::sys_info, level::error) << context_name << "\n\tget_current_process failed: " << GetLastError () << bs_end;
      }
    else
      {
        BOOL memory_get_result = GetProcessMemoryInfo (process, &pcmex, sizeof (pcmex));
        if (memory_get_result)
          {
            //BOSOUT (section::system_info, level::medium) << context_name << "\n"
            //  << "\tpeak ws size: " << (double)pcmex.PeakWorkingSetSize / (1024 * 1024) << "\n"
            //  << "\tws size: " << (double)pcmex.WorkingSetSize / (1024 * 1024) << "\n"
            //  << "\tpagefile usage: " << (double)pcmex.PeakPagefileUsage / (1024 * 1024) 
            //  << bs_end;
          }
        else
          {
            BOSERR (section::sys_info, level::error) << context_name << "\n\tget_process_memory_info failed: " << GetLastError () << bs_end;
          }
      }
#endif 
#endif
  }


} // namespace tools
} // namespace blue_sky


#endif  // #ifndef BS_TOOLS_PRINT_PROCESS_MEMORY_INFO_H_

