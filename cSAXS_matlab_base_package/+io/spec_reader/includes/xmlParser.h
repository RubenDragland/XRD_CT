/**
 * @file xmlParser.h
 * @author CXS group
 * @brief 
 * @version 0.1
 * @date 2019-07-22
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
#include "precomp.h"


#ifndef XMLPARSER_H
#define XMLPARSER_H
#include "specRead.h"
#include "spec_reader_utils.h"
// #include "Cpp.h"


/**
 * @brief Helper class to loop through an XML file and write the NEXUS file
 * 
 */
class xmlParser {
   public:
    std::string specDefaultPath;
    std::string collPath;

    xmlParser();
    template <typename h5GroupData>
    void writeAttributes(tinyxml2::XMLElement* element, std::string indent, h5GroupData group);

    void parseXML(tinyxml2::XMLElement* element, std::string indent, H5::Group group, std::map<std::string, std::vector<float>>& motor, std::map<std::string, std::vector<float>>& monitor, std::map<std::string, std::vector<std::string>>& metadataContainer, std::map<std::string, datCont_t> datContainer);
};
#endif