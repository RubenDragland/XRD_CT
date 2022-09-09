/**
 * @file spec_reader.cpp
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
#include "help.h"

using namespace std;

namespace fs = std::experimental::filesystem;

extern int cxs::io::GLOBAL_verbosity;
fs::path binaryName;


int main(int argc, char** argv) {
    binaryName = fs::path(std::string(argv[0]));
    globalArgs_t globalArgs = ProcessArgs(argc, argv);

    if (globalArgs.targetString.size() == 0) {
        std::cout << "The search string has to be specified!" << std::endl;
        display_usage();
        return 1;
    }

    std::string specFile = globalArgs.source;
    std::ifstream datFile(specFile);  //("specES1_started_2016_06_28_0954.dat");
    std::vector<std::string> inputFile;
    std::string line;

    // std::string scanNr_str;
    std::string startString = "#S ";
    std::string startMeta = globalArgs.targetString;//"#C meta " + targetMetaEntry + " int 1 " + targetMetaVal;

    std::string currentScan;
    std::vector<std::string> scans;
    std::string jsonOut;
    jsonOut = "[";
    datFile.seekg(0);
    while (datFile.peek() != EOF) {
        getline(datFile, line);
        if (line.compare(0,startString.size(),startString)==0){
            istringstream iss(line);
            std::vector<std::string> scanContent;
            scanContent = cxs::utils::split(line, ' ');
            currentScan = scanContent[1];
        }
        if ((line.compare(0,startMeta.size(),startMeta)==0) && (line.length()==startMeta.length()+1)){
            scans.push_back(currentScan);
            jsonOut += currentScan + ", ";
            // std::cout << currentScan << std::endl;
        }
    }
    if (jsonOut.length()>1){
        jsonOut.replace(jsonOut.length() - 2, jsonOut.length(), "]\n");
    } else {
        jsonOut += "]";
    }

    std::cout << jsonOut << std::endl;

    datFile.close();
}
