/**
 * @file help.cpp
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

namespace fs = std::experimental::filesystem;

extern int cxs::io::GLOBAL_verbosity;

void display_usage(void) {
    std::cout << std::endl
              << std::endl
              << binaryName.filename() << std::endl
              << std::endl;

    std::cout << "USAGE: " << std::endl;
    std::cout << "\t " + binaryName.string() + " [options]" << std::endl;
    std::cout << "ARGUMENTS: " << std::endl;
    std::cout << "\t -s | --source          SPEC source file" << std::endl;
    std::cout << "\t -o | --output          output file" << std::endl;
    std::cout << "\t -v | --verbosity       set verbosity (integer)" << std::endl;
    std::cout << "\t -h | --help            this usage info" << std::endl
              << std::endl;
    std::cout << "\t --scanNr               specify a single scan number or a range (e.g. 247-320)" << std::endl;
    std::cout << "\t --hdf5                 save output as HDF5" << std::endl;
    std::cout << "\t --orchestra			path to orchestra directory" << std::endl;
    std::cout << "\t --xmlLayout			path to the XML layout file for writing NEXUS files" << std::endl;

    exit(EXIT_FAILURE);
}

globalArgs_t ProcessArgs(int argc, char** argv) {
    globalArgs_t globalArgs;

    const char* const short_opts = "s:o:v:01234:h";
    const option long_opts[] = {
        {"source", required_argument, nullptr, 's'},
        {"output", required_argument, nullptr, 'o'},
        {"verbosity", required_argument, nullptr, 'v'},
        {"scanNr", required_argument, nullptr, '0'},
        {"hdf5", no_argument, nullptr, '1'},
        {"orchestra", required_argument, nullptr, '2'},
        {"xmlLayout", required_argument, nullptr, '3'},
        {"datFiles", required_argument, nullptr, '4'},
        {"help", no_argument, nullptr, 'h'},
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt) {
            case 's':
                try {
                    globalArgs.source = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'o':
                try {
                    globalArgs.output = fs::path(std::string(optarg));
                    fs::path outdir = globalArgs.output;

                    break;
                } catch (...) {
                    display_usage();
                }

            case 'v':
                try {
                    cxs::io::GLOBAL_verbosity = std::stoi(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case '0':
                try {
                    auto scanRange = cxs::utils::split(optarg, ',');
                    for (auto scanID : scanRange) {
                        auto scanValRange = cxs::utils::split(scanID, '-');
                        if (scanValRange.size() > 1) {
                            for (int ii = 0; std::stoi(scanValRange[0]) + ii <= std::stoi(scanValRange[scanValRange.size() - 1]); ii++) {
                                globalArgs.scanNr.push_back(std::stoi(scanValRange[0]) + ii);
                            }
                        } else {
                            globalArgs.scanNr.push_back(std::stoi(scanID));
                        }
                    }

                    break;
                } catch (...) {
                    display_usage();
                }
            case '1':
                try {
                    globalArgs.hdf5 = true;
                    break;
                } catch (...) {
                    display_usage();
                }
            case '2':
                try {
                    globalArgs.orchestra = optarg;
                    break;
                } catch (...) {
                    display_usage();
                }
            case '3':
                try {
                    globalArgs.xmlLayoutFile = optarg;
                    break;
                } catch (...) {
                    display_usage();
                }
            case '4':
                try {
                    globalArgs.datFiles = optarg;
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
    return globalArgs;
}
