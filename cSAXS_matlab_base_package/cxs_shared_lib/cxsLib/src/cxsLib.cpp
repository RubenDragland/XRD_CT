/**
 * @file cxsLib.cpp
 * @author CXS group
 * @brief collection of various helper functions that ought to be shared between programs
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
#include "cxsLib.h"

/**
 * @brief collection of various helper functions that ought to be shared between programs
 * 
 */
namespace cxs {

namespace utils {
/**
 * @brief split string by given delimiter and return a vector of strings
 * 
 */
template <typename Out>
void split(const std::string& s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) {
            *(result++) = item;
        }
    }
}
std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

/**
 * @brief search for given argument in json stream and return the matching value
 * 
 * @param jsonReply json input string 
 * @param searchVal parameter to search for
 * @return std::string output
 */
std::string extractFromJson(const std::string& jsonReply, const std::string& searchVal) {
    int searchValStart = jsonReply.find(searchVal) + searchVal.length() + 3;
    int searchValStop = jsonReply.find("\"", searchValStart);
    std::string jsonEntry;

    jsonEntry = jsonReply.substr(searchValStart, searchValStop - searchValStart);

    return jsonEntry;
}

/**
 * @brief check if file exists
 * 
 * @param name filename
 * @return true file exists
 * @return false file does not exist
 */
bool exists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

/**
 * @brief Construct a new cxs::utils::Timer::Timer object
 * 
 */
cxs::utils::Timer::Timer() {
    m_start = std::chrono::high_resolution_clock::now();
    
}
/**
 * @brief Destroy the cxs::utils::Timer::Timer object
 * Calls Stop() function if Timer instance goes out of scope and is still running
 * 
 */
cxs::utils::Timer::~Timer() {
    if (scopeBased)
        Stop();
}

/**
 * @brief Stop() member function to stop the current Timer instance
 * 
 */
void cxs::utils::Timer::Stop(){
    auto m_end = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::time_point_cast<std::chrono::nanoseconds> (m_start).time_since_epoch().count();
    auto stop = std::chrono::time_point_cast<std::chrono::nanoseconds> (m_end).time_since_epoch().count();
    auto duration = stop - start;

    elapsedTime = duration;
    scopeBased = false;
    auto out = getOutput();
    std::cout << "Elapsed time: " << std::get<0>(out) << " " << std::get<1>(out) << std::endl;
}

/**
 * @brief beautify output and adjust the unit
 * 
 * @return std::tuple<std::string, std::string> 
 */
std::tuple<std::string, std::string> cxs::utils::Timer::getOutput(){
    std::string time = std::to_string(elapsedTime);
    std::string unit = "ns";

    double tmpTime = elapsedTime;
    tmpTime/=1000000000;
    if (tmpTime>1){
        time = std::to_string(tmpTime);
        unit = "s";
        return std::make_tuple(time, unit);
    } 
    tmpTime*=1000;
    if (tmpTime>1){
        time = std::to_string(tmpTime);
        unit = "ms";
        return std::make_tuple(time, unit);
    } 
    tmpTime*=1000;
    if (tmpTime>1){
        time = std::to_string(tmpTime);
        unit = "us";
        return std::make_tuple(time, unit);
    } 
    

    return std::make_tuple(time, unit);
}

}  // namespace utils

/**
 * @brief collection of input/output-related helper functions
 * 
 */
namespace io {
/**
 * @brief 
 * 
 * @param lvl 
 * @param out 
 */

/**
 * @brief write to std::cout depending on given verbosity level
 * 
 * @param lvl verbosity level
 * @param out string that needs to written to cout
 */
#ifndef VERBOSE_H
#define VERBOSE_H
int GLOBAL_verbosity = 1;
#endif /* VERBOSE_H */
void verbose(const int lvl, const std::string& out) {
    if (lvl <= cxs::io::GLOBAL_verbosity) {
        std::cout << out << std::endl;
    }
}

/**
 * @brief turn bold font on
 * 
 * @return std::string 
 */
std::string bold_on() {
    return "\e[1m";
}

/**
 * @brief turn bold font off
 * 
 * @return std::string 
 */
std::string bold_off() {
    return "\e[0m";
}

/**
 * @brief system call with a single return value
 * 
 * @param bash_cmd 
 * @return std::string 
 */
std::string system_call(const std::string& bash_cmd) {
   // local variable declaration
    redi::ipstream is(bash_cmd);
    std::string cmd_res;
    std::getline(is, cmd_res);
    is.close();
return cmd_res;
}

/**
 * @brief system call with multiple return values, parses as vector of strings
 * 
 * @param bash_cmd 
 * @return std::vector<std::string> 
 */
std::vector<std::string> system_call_multiple(const std::string& bash_cmd) {
   // local variable declaration
    redi::ipstream is(bash_cmd);
    std::vector<std::string> cmd_res;
    std::string cmd_res_inter;
    while (std::getline(is, cmd_res_inter))
    {
        cmd_res.push_back (cmd_res_inter);
        //std::cout << cmd_res_inter;
    }
    is.close();
return cmd_res;
}

}  // namespace io

}  // namespace cxs
