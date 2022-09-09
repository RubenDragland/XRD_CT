/**
 * @file cxsLib.h
 * @author CXS group
 * @brief 
 * @version 0.1
 * @date 2019-07-19
 * 
 * *-----------------------------------------------------------------------*
 * |                                                                       |
 * |  Except where otherwise noted, this work is licensed under a          |
 * |  Creative Commons Attribution-NonCommercial-ShareAlike 4.0            |
 * |  International (CC BY-NC-SA 4.0) license.                             |
 * |                                                                       |
 * |  Copyright (c) 2017 by Paul Scherrer Institute (http://www.psi.ch)    |
 * |                                                                       |
 * |       Author: CXS group, PSI                                          |
 * *-----------------------------------------------------------------------*
 * You may use this code with the following provisions:
 *                                                                            
 * If the code is fully or partially redistributed, or rewritten in another
 *   computing language this notice should be included in the redistribution.
 * 
 * If this code, or subfunctions or parts of it, is used for research in a 
 *   publication or if it is fully or partially rewritten for another 
 *   computing language the authors and institution should be acknowledged 
 *  in written form in the publication: “Data processing was carried out 
 *   using the 'cSAXS software package' developed by the CXS group,
 *   Paul Scherrer Institut, Switzerland.”
 *   Variations on the latter text can be incorporated upon discussion with
 *   the CXS group if needed to more specifically reflect the use of the package
 *   for the published work.
 * 
 * A publication that focuses on describing features, or parameters, that
 *    are already existing in the code should be first discussed with the
 *    authors.
 *  
 * This code and subroutines are part of a continuous development, they 
 *    are provided “as they are” without guarantees or liability on part
 *    of PSI or the authors. It is the user responsibility to ensure its 
 *    proper use and the correctness of the results.
 * 
 */
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <tuple>
#include "pstream.h"

/**
 * @brief collection of various helper function that ought to be shared between programs
 * 
 */
namespace cxs {
namespace utils {
/**
 * @brief split string by given delimiter and return a vector of strings
 * 
 */
template <typename Out>
void split(const std::string& s, char delim, Out result);
std::vector<std::string> split(const std::string& s, char delim);

/**
 * @brief check if file exists
 * 
 * @param name filename
 * @return true file exists
 * @return false file does not exist
 */
bool exists(const std::string& name);

/**
 * @brief search for given argument in json stream and return the matching value
 * 
 * @param jsonReply json input string 
 * @param searchVal parameter to search for
 * @return std::string output
 */
std::string extractFromJson(const std::string& jsonReply, const std::string& searchVal);

/**
 * @brief Timer class
 * Calculate elapsed time with nanoseconds precision. 
 * E.g.:
 *  Timer tic = Timer();
 *  // do something meaningful //
 *  tic.Stop();
 * Timer's destructor calls Stop() automatically. Like this, a scope based timer can be implemented without explicitly calling <timer_instance>.Stop();
 * E.g.:
 *  {
 *      Timer tic = Timer();
 *      // do something meaningful //
 *  }
 * 
 */
#ifndef TIMER_H
#define TIMER_H
class Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::tuple<std::string, std::string> getOutput();
    bool scopeBased = true;

    double elapsedTime;

    public:
        Timer();
        ~Timer();
        void Stop();

};
#endif
}  // namespace utils

namespace io {

extern int GLOBAL_verbosity;
/**
 * @brief write to std::cout depending on given verbosity level
 * 
 * @param lvl vebosity level
 * @param out string that needs to written to cout
 */
void verbose(const int lvl, const std::string& out);

/**
 * @brief turn bold font on
 * 
 * @return std::string 
 */
std::string bold_on();
/**
 * @brief turn bold font off
 * 
 * @return std::string 
 */
std::string bold_off();

/**
 * @brief system call with a single return value
 * 
 * @param bash_cmd 
 * @return std::string 
 */
std::string system_call(const std::string& bash_cmd);

/**
 * @brief system call with multiple return values, parses as vector of strings
 * 
 * @param bash_cmd 
 * @return std::vector<std::string> 
 */
std::vector<std::string> system_call_multiple(const std::string& bash_cmd);

}  // namespace io
}  // namespace cxs