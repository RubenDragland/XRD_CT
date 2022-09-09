/**
 * @file help.cpp
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
#include "help.h"

/**
 * @brief display help string / documentation
 * 
 */
void display_usage(void) {
    std::cout << std::endl
              << std::endl
              << binaryName.filename() << std::endl
              << std::endl;

    std::cout << "USAGE: " << std::endl;
    std::cout << "\t " + binaryName.string() + " [options]" << std::endl;
    std::cout << "ARGUMENTS: " << std::endl;
    std::cout << "\t -s | --source          input directory" << std::endl;
    std::cout << "\t -o | --output          output file" << std::endl;
    std::cout << "\t -p | --path            path within the nexus file" << std::endl;


    exit(EXIT_FAILURE);
}

/**
 * @brief parser for handling input arguments
 * 
 * @param argc 
 * @param argv 
 * @return globalArgs_t 
 */
globalArgs_t ProcessArgs(int argc, char** argv) {
    globalArgs_t globalArgs;

    const char* const short_opts = "s:o:p:h";
    const option long_opts[] = {
        {"source", required_argument, nullptr, 's'},
        {"output", required_argument, nullptr, 'o'},
        {"path", required_argument, nullptr, 'p'},
        {"help", no_argument, nullptr, 'h'},
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt) {
            case 's':
                try {
                    globalArgs.inputDir = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'o':
                try {
                    globalArgs.fname = std::string(optarg);
                    break;
                } catch (...) {
                    display_usage();
                }

            case 'p':
                try {
                    globalArgs.path = std::string(optarg);
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