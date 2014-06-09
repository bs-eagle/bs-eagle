#include "bs_common.h"
#include "date_helper.h"

#include <stdio.h>
#include <string>
#include <time.h>
//#include <boost/lexical_cast.hpp>
#include <cstdlib>

void int2str (int a, char *str)
{
	if (a < 10)
		sprintf(str, "0%d", a);
	else
		sprintf(str, "%d", a);
}


// return date time as string in format dd.MM.yyyy hh:mm:ss
std::string get_date_time_str ()
{
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	char str[1024];
	int2str(timeinfo->tm_mday, &str[0]);
	int2str(1 + timeinfo->tm_mon, &str[3]);
	int2str(1900 + timeinfo->tm_year, &str[6]);
	int2str(timeinfo->tm_hour, &str[11]);
	int2str(timeinfo->tm_min, &str[14]);
	int2str(timeinfo->tm_sec, &str[17]);
	str[2] = str[5] = '.';
	str[10] = '_';
	str[13] = str[16] = ':';

	return std::string (str);
}

std::vector< int > extract_date_time(const std::string& dt) {
	using namespace boost;
	// assume that strings were formatted by get_date_time_str()
	// extract data in following order:
	// year, month, day, hour, min, sec
	std::vector< int > res(6, 0);
	if(dt.size() < 19)
		return res;
	// take last 19 symbols
	const std::string dt_ = dt.substr(dt.size() - 19);

	// atoi is much faster than lexical_cast that is an overkill here
	// year
	res[0] = std::atoi(dt_.substr(6, 4).c_str());
	// month
	res[1] = std::atoi(dt_.substr(3, 2).c_str());
	// day
	res[2] = std::atoi(dt_.substr(0, 2).c_str());
	// hour
	res[3] = std::atoi(dt_.substr(11, 2).c_str());
	// min
	res[4] = std::atoi(dt_.substr(14, 2).c_str());
	// sec
	res[5] = std::atoi(dt_.substr(17, 2).c_str());

	//// year
	//res[0] = lexical_cast< int >(dt_.substr(6, 4));
	//// month
	//res[1] = lexical_cast< int >(dt_.substr(3, 2));
	//// day
	//res[2] = lexical_cast< int >(dt_.substr(0, 2));
	//// hour
	//res[3] = lexical_cast< int >(dt_.substr(11, 2));
	//// min
	//res[4] = lexical_cast< int >(dt_.substr(14, 2));
	//// sec
	//res[5] = lexical_cast< int >(dt_.substr(17, 2));
	return res;
}

// implementation of SGI extension algorithm
template< typename _InputIterator1, typename _InputIterator2 >
int lexicographical_compare_3way(_InputIterator1 __first1,
		_InputIterator1 __last1,
		_InputIterator2 __first2,
		_InputIterator2 __last2)
{
	while (__first1 != __last1 && __first2 != __last2)
	{
		if (*__first1 < *__first2)
			return -1;
		if (*__first2 < *__first1)
			return 1;
		++__first1;
		++__first2;
	}
	if (__first2 == __last2)
		return !(__first1 == __last1);
	else
		return -1;
}

bool str_dt_compare::operator()(const std::string& dt1, const std::string& dt2) const {
	const std::vector< int >& dt1_dig = extract_date_time(dt1);
	const std::vector< int >& dt2_dig = extract_date_time(dt2);

	return lexicographical_compare_3way(
		dt1_dig.begin(), dt1_dig.end(), dt2_dig.begin(), dt2_dig.end()
	) < 0;
}

