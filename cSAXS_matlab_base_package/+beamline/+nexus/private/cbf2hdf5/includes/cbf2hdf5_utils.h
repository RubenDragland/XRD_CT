/**
 * @file cbf2hdf5_utils.h
 * @author CXS group
 * @brief 
 * @version 0.1
 * @date 2019-10-31
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
#include<vector>
#include<string>
#include <fstream>
#include <cstring>
#include <regex>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <zlib.h>
#include <cmath>

#include "hdf5.h"
#include "H5DOpublic.h"


void pilatus_read_data(std::string &path, uint fileID, hid_t dataset, hid_t datatype, hid_t dataspace, hsize_t *dim, uint aggression, uint64_t *timeStamp);

void cbf_dims (const std::vector<char> &data, std::vector<long> &dims);

void appendUnits(hid_t gid, const std::string& units);