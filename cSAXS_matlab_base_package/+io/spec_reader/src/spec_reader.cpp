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
#include "specRead.h"

using namespace std;
using namespace tinyxml2;
using namespace H5;

namespace fs = std::experimental::filesystem;

extern int cxs::io::GLOBAL_verbosity;
fs::path binaryName;


int main(int argc, char** argv) {
    binaryName = fs::path(std::string(argv[0]));
    globalArgs_t globalArgs = ProcessArgs(argc, argv);

    if (globalArgs.scanNr.size() == 0) {
        std::cout << "The scan number has to be specified!" << std::endl;
        display_usage();
        return 1;
    }

    std::string specFile = globalArgs.source;
    std::ifstream datFile(specFile);  //("specES1_started_2016_06_28_0954.dat");
    std::vector<std::string> inputFile;
    std::string line;

    // std::string scanNr_str;
    // std::string startString;
    // std::string endString;
    // std::string motorString;
    // bool readStatus = false;
    // bool lineRead = false;

    datFile.seekg(0);
    while (datFile.peek() != EOF) {
        getline(datFile, line);
        // if (line.length()>0){
        inputFile.push_back(line);
        // }

        // lineRead = false;
        // // inputFile.push_back(line);
        // for (uint scansID=0; (scansID < globalArgs.scanNr.size() && !lineRead); scansID++){
        // 	scanNr_str = std::to_string(globalArgs.scanNr[scansID]);
        // 	startString = "#S " + scanNr_str;
        // 	endString = "#X " + scanNr_str;
        // 	motorString = "#O";
        //
        // 	if (line.compare(0,motorString.size(),motorString)==0 && !lineRead){
        // 		inputFile.push_back(line);
        // 		lineRead = true;
        // 	}
        //
        // 	if (readStatus==true && line.compare(0,endString.size(),endString)==0 && !lineRead){
        // 		inputFile.push_back(line);
        // 		readStatus=false;
        // 		lineRead = true;
        // 	}
        //
        // 	if (readStatus==true && !lineRead){
        // 		inputFile.push_back(line);
        // 		lineRead = true;
        // 	}
        //
        // 	if (readStatus==false && line.compare(0,startString.size(),startString)==0 && !lineRead){
        // 		readStatus=true;
        // 		inputFile.push_back(line);
        // 		lineRead = true;
        // 	}
        //
        // }
    }
    // cout << line << endl;

    std::vector<std::unique_ptr<specRead>> storage;

    storage.reserve(globalArgs.scanNr.size() + 1);
    // std::cout << storage.max_size() << std::endl;

    if (datFile.is_open()) {
        #pragma omp critical
        for (uint scanNrIdx = 0; scanNrIdx < globalArgs.scanNr.size(); scanNrIdx++) {
            storage.push_back(std::make_unique<specRead>(globalArgs.scanNr[scanNrIdx], datFile));
            // std::cout << globalArgs.scanNr[scanNrIdx] << std::endl;
        }

        #pragma omp parallel for
        for (uint scanNrIdx = 0; scanNrIdx < globalArgs.scanNr.size(); scanNrIdx++) {
            storage[scanNrIdx]->read_from_file(inputFile);
            storage[scanNrIdx]->read_from_orchestra(globalArgs.orchestra, false, "orchestra");
            storage[scanNrIdx]->read_datFiles(globalArgs.datFiles);
            storage[scanNrIdx]->output_stream(globalArgs.hdf5, globalArgs.output, scanNrIdx, globalArgs.scanNr.size(), globalArgs.xmlLayoutFile);
        }
        if (!globalArgs.hdf5 && globalArgs.output.length() == 0) {
            // send to terminal
            // #pragma omp critical
            for (uint scanNrIdx = 0; scanNrIdx < globalArgs.scanNr.size(); scanNrIdx++) {
                std::cout << storage[scanNrIdx]->jsonOut;
            }
        } 

    } else {
        cout << "Unable to open file" << endl;
    }
    datFile.close();
}
