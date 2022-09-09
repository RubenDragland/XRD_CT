/**
 * @file help.cpp
 * @author CXS group
 * @brief define help function and input argument parser
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
 *   using the 'cSAXS matlab package' developed by the CXS group,
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

namespace fs = std::filesystem;

extern int GLOBAL_verbosity;
extern int GLOBAL_debug;

void display_usage(void) {
    std::cout << std::endl
              << std::endl
              << binaryName.filename() << std::endl
              << std::endl;

    std::cout << cxs::io::bold_on() + "USAGE: " + cxs::io::bold_off() << std::endl;
    std::cout << "\t " + binaryName.string() + " [options]" << std::endl;
    std::cout << cxs::io::bold_on() + "ARGUMENTS: " + cxs::io::bold_off() << std::endl;
    std::cout << "\t GENERAL SETTINGS: " << std::endl;
    std::cout << "\t -d | --data            full path of a nexus raw data file" << std::endl;
    std::cout << "\t -p | --path            path inside the nexus file pointing to a 2D/3D dataset" << std::endl;
    std::cout << "\t -m | --mask            full path of the valid mask file" << std::endl;
    std::cout << "\t -o | --output          full path of the output file" << std::endl;
    std::cout << "\t -w | --weights         use weights" << std::endl; 
    std::cout << "\t -h | --help            this usage info" << std::endl

    << std::endl;
    exit(EXIT_FAILURE);
}

globalArgs_t ProcessArgs(int argc, char** argv) {
    globalArgs_t globalArgs;

    const char* const short_opts = "d:p:m:o:wh";
    const option long_opts[] = {
        {"data", required_argument, nullptr, 'd'},
        {"path", required_argument, nullptr, 'p'},
        {"mask", required_argument, nullptr, 'm'},
        {"output", required_argument, nullptr, 'o'},
        {"weights", no_argument, nullptr, 'w'},
        {"help", no_argument, nullptr, 'h'},
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt) {
            case 'd':
                try {
                    globalArgs.dataPath = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'p':
                try {
                    globalArgs.nexusPath = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'm':
                try {
                    globalArgs.maskPath = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'o':
                try {
                    globalArgs.outputPath = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }
            case 'w':
                try {
                    globalArgs.useWeights = true;
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'h':  // -h or --help
            case '?':  // Unrecognized option
                display_usage();
                break;
            default:
                display_usage();
                break;
        }
    }

    // some initial checks
    if (globalArgs.dataPath.length() == 0) {
        std::cout << cxs::io::bold_on() + "Please specify the full path to the nexus raw data file using the --data flag, e.g. --data /sls/X12SA/Data10/e12345/data/myDataset.h5." + cxs::io::bold_off() << std::endl;
        display_usage();
    }
    if (globalArgs.nexusPath.length() == 0) {
        std::cout << cxs::io::bold_on() + "Please specify the HDF5 path to a dataset inside the nexus raw data file using the --path flag, e.g.  --path /entry/instrument/pilatus_1/data." + cxs::io::bold_off() << std::endl;
        display_usage();
    }
    if (globalArgs.maskPath.length() == 0) {
        std::cout << cxs::io::bold_on() + "Please specify the full path to the mask file using the --mask flag." + cxs::io::bold_off() << std::endl;
        display_usage();
    }    
    if (globalArgs.dataPath.length() == 0) {
        std::cout << cxs::io::bold_on() + "Please specify the full path for the output file using the --output flag, e.g. --output /sls/X12SA/Data10/e12345/analysis/integData_01.h5." + cxs::io::bold_off() << std::endl;
        display_usage();
    }
    return globalArgs;
}
