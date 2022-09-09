/**
 * @file cbf2hdf5.cpp
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

#include <omp.h>
#include <filesystem>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "cbf2hdf5_utils.h"
#include "cxsLib.h"
#include "help.h"

namespace fs = std::filesystem;
fs::path binaryName;

int main(int argc, char** argv) {
    binaryName = fs::path(std::string(argv[0]));
    globalArgs_t globalArgs = ProcessArgs(argc, argv);

    std::vector<std::string> pilatusFiles;
    std::string detectorPath = globalArgs.inputDir;  //e.g. "/das/work/units/csaxs/p16812/data/pilatus/e16403/pilatus_1/S00000-00999/S00012/";
    std::string path = globalArgs.path;              // e.g. /entry/instrument/pilatus_1

    for (auto& entry : fs::directory_iterator(detectorPath)) {
        std::string centry = entry.path();
        if (centry.compare(centry.length() - 3, 3, "cbf") == 0) {
            pilatusFiles.push_back(centry);
        }
    }

    std::sort(pilatusFiles.begin(), pilatusFiles.end());

    hid_t datatype, dataspace, dataprop, dataset;
    hid_t datatypeTimeStamp, dataspaceTimeStamp, datasetTimeStamp; 
    hid_t fid, gid, fapl;
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
    bool fexists = false;
    // open new file if the container does not exist yet
    fid = H5Fopen(globalArgs.fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (fid < 0) {
        fid = H5Fcreate(globalArgs.fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    } else {
        fexists = true;
    }
    H5Pclose(fapl);

    // create H5 file structure
    std::vector<std::string> groups = cxs::utils::split(path, '/');
    if (fexists) {
        hid_t prevGroup = fid;

        for (auto en : groups) {
            gid = H5Gopen2(prevGroup, en.c_str(), H5P_DEFAULT);
            if (gid < 0) {
                gid = H5Gcreate2(prevGroup, en.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            }
            prevGroup = gid;
        }
    } else {
        hid_t prevGroup = fid;
        for (auto en : groups) {
            gid = H5Gcreate2(prevGroup, en.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            std::cout << en << std::endl;
            prevGroup = gid;
        }
    }


    // read the first file to get the header and datasize dimensions
    std::vector<char> fbuf;
    {
        std::string fPath = pilatusFiles[0];
        FILE* fin = fopen(fPath.c_str(), "r");
        if (!fin)
            throw std::runtime_error(std::string("unable to open file ") + fPath + ": " + std::strerror(errno));
        try {
            off_t fsz;
            {
                struct stat sbuf;
                if (fstat(fileno(fin), &sbuf) == -1) {
                    throw std::runtime_error(std::string("unable to stat file ") + fPath + ": " + std::strerror(errno));
                }
                fsz = sbuf.st_size;
            }
            fbuf.resize(fsz);
            fread(&fbuf[0], 1, fsz, fin);
            if (ferror(fin))
                throw std::runtime_error(std::string("unable to read file ") + fPath + ": " + std::strerror(errno));
            fclose(fin);
        } catch (...) {
            fclose(fin);
            throw;
        }
    }
    std::vector<long> fdims{0, 0};
    cbf_dims(fbuf, fdims);
    unsigned long finger;  // compressed data index
    {
        std::vector<char> sig{'\x0c', '\x1a', '\x04', '\xd5'};
        auto p = std::search(fbuf.begin(), fbuf.end(), sig.begin(), sig.end());
        if (p == fbuf.end())
            throw std::runtime_error(std::string("data signature not found within file " + path));
        finger = p - fbuf.begin() + 4;
    }
    std::string header;
    for (uint ii=0; ii<finger-4; ii++){
        header += fbuf[ii];
    }
    float exposureTime = -1;
    float exposurePeriod = -1;
    double tau = -1;
    int countCutoff = -1;
    int thresholdSetting = -1;
    std::string gainSettings;
    std::string trimFile;

    std::vector<std::string> headerSplit = cxs::utils::split(header, '\n');
    for (auto h: headerSplit){
        std::vector<std::string> tmp;
        if (h.find("# Exposure_time") != std::string::npos) {
            tmp = cxs::utils::split(h, ' ');
            exposureTime = std::strtof(tmp[2].c_str(), NULL);
        } else if (h.find("# Exposure_period") != std::string::npos) {
            tmp = cxs::utils::split(h, ' ');
            exposurePeriod = std::strtof(tmp[2].c_str(), NULL);
        } else if (h.find("# Tau") != std::string::npos) {
            tmp = cxs::utils::split(h, ' ');
            tau = std::strtod(tmp[3].c_str(), NULL);
        } else if (h.find("# Count_cutoff") != std::string::npos) {
            tmp = cxs::utils::split(h, ' ');
            countCutoff = (int)std::strtof(tmp[2].c_str(), NULL);
        } else if (h.find("# Threshold_setting") != std::string::npos) {
            tmp = cxs::utils::split(h, ' ');
            thresholdSetting = std::strtof(tmp[2].c_str(), NULL);
        } else if (h.find("# Gain_setting") != std::string::npos) {
            gainSettings = h.substr(16, h.length()-17);
        } else if (h.find("# Trim_file") != std::string::npos) {
            tmp = cxs::utils::split(h, ' ');
            trimFile = tmp[2].substr(0, tmp[2].length()-1);
        }
    }
    {
        hid_t metadataDataset, metadataDataspace;
        hsize_t dim[1] = {1};
        metadataDataspace = H5Screate_simple(1, dim, NULL);
        metadataDataset = H5Dcreate2(gid, "count_time", H5T_NATIVE_FLOAT, metadataDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(metadataDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &exposureTime);
        H5Sclose(metadataDataspace);
        appendUnits(metadataDataset, "s");
        H5Dclose(metadataDataset);

        metadataDataspace = H5Screate_simple(1, dim, NULL);
        metadataDataset = H5Dcreate2(gid, "count_period", H5T_NATIVE_FLOAT, metadataDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(metadataDataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &exposurePeriod);
        H5Sclose(metadataDataspace);
        appendUnits(metadataDataset, "s");
        H5Dclose(metadataDataset);

        metadataDataspace = H5Screate_simple(1, dim, NULL);
        metadataDataset = H5Dcreate2(gid, "tau", H5T_NATIVE_DOUBLE, metadataDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(metadataDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tau);
        H5Sclose(metadataDataspace);
        appendUnits(metadataDataset, "s");
        H5Dclose(metadataDataset);

        metadataDataspace = H5Screate_simple(1, dim, NULL);
        metadataDataset = H5Dcreate2(gid, "count_cutoff", H5T_NATIVE_INT, metadataDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(metadataDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &countCutoff);
        H5Sclose(metadataDataspace);
        appendUnits(metadataDataset, "counts");
        H5Dclose(metadataDataset);

        metadataDataspace = H5Screate_simple(1, dim, NULL);
        metadataDataset = H5Dcreate2(gid, "threshold_setting", H5T_NATIVE_INT, metadataDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(metadataDataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &thresholdSetting);
        H5Sclose(metadataDataspace);
        appendUnits(metadataDataset, "eV");
        H5Dclose(metadataDataset);

        {
            const char* s[1] = {gainSettings.c_str()};
            hid_t filetype, memtype, space, dset;
            filetype = H5Tcopy(H5T_FORTRAN_S1);
            H5Tset_size(filetype, H5T_VARIABLE);
            memtype = H5Tcopy(H5T_C_S1);
            H5Tset_size(memtype, H5T_VARIABLE);
            space = H5Screate(H5S_SCALAR);
            dset = H5Dcreate2(gid, "gain", filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)s);
            H5Sclose(space);
            H5Dclose(dset);
        }

        {
            const char* s[1] = {trimFile.c_str()};
            hid_t filetype, memtype, space, dset;
            filetype = H5Tcopy(H5T_FORTRAN_S1);
            H5Tset_size(filetype, H5T_VARIABLE);
            memtype = H5Tcopy(H5T_C_S1);
            H5Tset_size(memtype, H5T_VARIABLE);
            space = H5Screate(H5S_SCALAR);
            dset = H5Dcreate2(gid, "trim_file", filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)s);
            H5Sclose(space);
            H5Dclose(dset);
        }       
    }

    // prepare dataset
    int rank = 3;
    hsize_t dim[3] = {(hsize_t)pilatusFiles.size(), (hsize_t)fdims[0], (hsize_t)fdims[1]};
    datatype = H5T_STD_I32LE;

    // compression
    dataprop = H5Pcreate(H5P_DATASET_CREATE);
    /* Dataset must be chunked for compression */
    hsize_t cdims[3] = {1, dim[1], dim[2]};
    H5Pset_chunk(dataprop, rank, cdims);
    int aggression = 2;
    H5Pset_deflate(dataprop, aggression);

    hsize_t maxdim[3] = {H5S_UNLIMITED, dim[1], dim[2]};
    dataspace = H5Screate_simple(rank, dim, maxdim);

    char dsetName[] = "data";
    dataset = H5Dcreate2(gid, dsetName, datatype, dataspace, H5P_DEFAULT, dataprop, H5P_DEFAULT);

    uint64_t* timeStamp = new uint64_t[pilatusFiles.size()];

    if (dataset > 0) {
        #pragma omp parallel for
        for (uint fileID = 0; fileID < pilatusFiles.size(); fileID++) {
            pilatus_read_data(pilatusFiles[fileID], fileID, dataset, datatype, dataspace, dim, aggression, timeStamp);
        }
    }

    // save timestamps
    hsize_t dimTimeStamp[1] = {(hsize_t)pilatusFiles.size()};
    datatypeTimeStamp = H5T_STD_U64LE;
    dataspaceTimeStamp = H5Screate_simple(1, dimTimeStamp, NULL);
    char dsetNameTimeStamp[] = "timestamp";
    datasetTimeStamp = H5Dcreate2(gid, dsetNameTimeStamp, datatypeTimeStamp, dataspaceTimeStamp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(datasetTimeStamp, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, timeStamp);
    appendUnits(datasetTimeStamp, "YYYYMMDDHHMMSS.SSS");
    
    delete [] timeStamp;

}
