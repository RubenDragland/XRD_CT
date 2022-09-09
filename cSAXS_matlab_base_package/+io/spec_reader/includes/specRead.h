/**
 * @file specRead.h
 * @author CXS group
 * @brief Main class to parse a SPEC file
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
#include "help.h"
#include "spec_reader_utils.h"


#ifndef SPECREAD_H
#define SPECREAD_H
class specRead {
    uint lineidx;
    bool scanexists;
    std::string line;
    std::string scanNr_str;
    std::string startString;  //("#S 5398");
    std::string endString;    //("#X 5398");
    std::string startScanString;
    std::string motorString;
    std::vector<std::vector<std::string>> scanContent;
    std::vector<std::string> scanContentName;
    std::vector<std::string> motorConfig;
    std::vector<std::string> motorUpdate;
    std::vector<std::string> metadata;


   public:
    int scan;
    std::string jsonOut;

    std::vector<std::string> header;
    std::map<std::string, std::vector<float>> motor;
    std::map<std::string, std::vector<float>> monitor;
    std::map<std::string, std::vector<std::string>> metadataContainer;
    std::map<std::string, datCont_t> datContainer;

    specRead(int& scanNr, std::ifstream& datFile);
    void read_from_file(const std::vector<std::string>& inputFile);
    void read_from_orchestra(const std::string& orchestra_path, bool fullPath, std::string groupName);
    void read_datFiles(const std::string& datFiles);
    void output_stream(const bool hdf5, const std::string output, const uint& scanNrIdx, const uint scanNrSz, const std::string& xmlLayoutFile);
    void flush();

};
#endif
