#include <string>
#include <time.h>

void int2str (int a, char *str);

// return date time as string in format dd.MM.yyyy hh:mm:ss
std::string get_date_time_str ();

// properly compare strings formatted by function above
struct str_dt_compare {
	bool operator()(const std::string& dt1, const std::string& dt2) const;
};

