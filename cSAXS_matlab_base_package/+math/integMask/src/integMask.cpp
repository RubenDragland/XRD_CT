
#include <omp.h>
#include <filesystem>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <future>
#include <numeric>   
#include <math.h> 

#include "help.h"
#include "hdf5.h"
#include "H5DOpublic.h"
#include <zlib.h>
#define DEFLATE_SIZE_ADJUST(s) (ceil(((double)(s))*1.001)+12)



namespace fs = std::filesystem;
fs::path binaryName;

void copy_mask_content(globalArgs_t globalArgs){
        // output file
    hid_t fapl, fileOut;
    
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);
    fileOut = H5Fcreate((globalArgs.outputPath + ".TMP").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    H5Fclose(fileOut);

    std::system(("h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /angular_segments -d /angular_segments; h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /center_xy -d /center_xy; h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /no_of_segments -d /no_of_segments; h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /phi_det -d /phi_det; h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /integ_masks/radius -d /radius; h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /integ_masks/norm_sum -d /norm_sum; h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /integ_masks/q -d /q; h5copy -p -i " + globalArgs.maskPath + " -o " + globalArgs.outputPath + ".TMP" + " -s /E -d /E >>/dev/null 2>>/dev/null").c_str());
}

int main(int argc, char** argv) {

    binaryName = fs::path(std::string(argv[0]));
    globalArgs_t globalArgs = ProcessArgs(argc, argv);

    // dataset dimensions: e.g. [201, 1679, 1475]
    hid_t fileIn, fileInMask, dset, dsetPointer, dsetIndex, dsetWeight; 

    // input file (nexus)
    fileIn = H5Fopen(globalArgs.dataPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen2(fileIn, globalArgs.nexusPath.c_str(), H5P_DEFAULT);
    hid_t dspace = H5Dget_space(dset);
    const int ndims = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);

    // input file (integ mask)
    fileInMask = H5Fopen(globalArgs.maskPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    // get pointer
    dsetPointer = H5Dopen2(fileInMask, "/integ_masks/indices/maskPointer", H5P_DEFAULT);
    hid_t dspacePointer = H5Dget_space(dsetPointer);
    const int ndimsPointer = H5Sget_simple_extent_ndims(dspacePointer);
    hsize_t dimsPointer[ndimsPointer];
    H5Sget_simple_extent_dims(dspacePointer, dimsPointer, NULL);
    uint64_t dimsTotalPointer = dimsPointer[0]*dimsPointer[1]*dimsPointer[2];
    if (dimsPointer[2]!=2)
        throw std::runtime_error("maskPointer must have dimension 2 x seg x radii.");
    auto maskPointer = new int[dimsTotalPointer];
    H5Dread(dsetPointer, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, maskPointer);

    // get indices
    dsetIndex = H5Dopen2(fileInMask, "/integ_masks/indices/maskIndex", H5P_DEFAULT);
    hid_t dspaceIndex = H5Dget_space(dsetIndex);
    hsize_t dimsIndex[1];
    H5Sget_simple_extent_dims(dspaceIndex, dimsIndex, NULL);
    auto maskIndex = new int[dimsIndex[0]];
    H5Dread(dsetIndex, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, maskIndex);

    // get weights
    dsetWeight = H5Dopen2(fileInMask, "/integ_masks/indices/maskWeight", H5P_DEFAULT);
    hid_t dspaceWeight = H5Dget_space(dsetWeight);
    hsize_t dimsWeight[1];
    H5Sget_simple_extent_dims(dspaceWeight, dimsWeight, NULL);
    auto maskWeight = new float[dimsWeight[0]];
    if (globalArgs.useWeights) {
        H5Dread(dsetWeight, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, maskWeight);
    }



    // copy stuff
    auto copyMask = std::async(copy_mask_content, globalArgs);

    auto outputMean = new float[dimsPointer[0] * dimsPointer[1] * dims[0]];
    auto outputStdDev = new float[dimsPointer[0] * dimsPointer[1] * dims[0]];
    size_t buf_size = dims[1] * dims[2] * sizeof(unsigned int);

    #pragma omp parallel for
    for (uint sliceID = 0; sliceID < dims[0]; sliceID++) {
        auto arr = new int[dims[1] * dims[2]];
        hsize_t start[3] = {0, 0, 0};
        start[0] = sliceID;
        hsize_t chunk_nbytes;
        /* Get the size of the compressed chunk */
        Bytef* pt_readbuf;    /* Point to the buffer for data read */
        unsigned filter_mask = 0;

        H5Dget_chunk_storage_size(dset, start, &chunk_nbytes);

        pt_readbuf = new Bytef[chunk_nbytes];

        /* Use H5DOread_chunk() to read the chunk back */
        H5DOread_chunk(dset, H5P_DEFAULT, start, &filter_mask, pt_readbuf);

        uncompress((Bytef*)(arr), (uLongf*)&buf_size, pt_readbuf, (uLong)chunk_nbytes);
        delete[] pt_readbuf;
        
        // azimuthal integration
        for (uint ind_r = 0; ind_r < dimsPointer[0]; ind_r++) {
            for (uint ind_seg = 0; ind_seg < dimsPointer[1]; ind_seg++) {
                double tmpMean = 0;
                double tmpStdDev = 0;
                int *dim = &maskPointer[ind_r*dimsPointer[2]*dimsPointer[1]+ind_seg*dimsPointer[2]+1];
                
                if (*dim > 0) {
                    if (globalArgs.useWeights) {
                        for (int ii = 0; ii < *dim; ii++) {
                            tmpMean += arr[maskIndex[*(dim - 1) + ii] - 1] * maskWeight[maskIndex[*(dim - 1) + ii] - 1];
                        }
                    } else {
                        for (int ii = 0; ii < *dim; ii++) {
                            tmpMean += arr[maskIndex[*(dim - 1) + ii] - 1];
                        }
                    }

                    if (tmpMean > 0) {
                        tmpMean /= *dim;
                    }
                    if (globalArgs.useWeights) {
                        for (int ii = 0; ii < *dim; ii++) {
                            tmpStdDev += std::pow(std::abs(arr[maskIndex[*(dim - 1) + ii] - 1] * maskWeight[maskIndex[*(dim - 1) + ii] - 1] - tmpMean), 2);
                        }
                    } else {
                        for (int ii = 0; ii < *dim; ii++) {
                            tmpStdDev += std::pow(std::abs(arr[maskIndex[*(dim - 1) + ii] - 1] - tmpMean), 2);
                        }
                    }
                    

                    if (tmpMean > 0 && *dim > 1) {
                        tmpStdDev = sqrt(tmpStdDev / (*dim - 1));
                    }

                } else {
                    tmpMean = -1;
                    tmpStdDev = -1;
                }

                outputMean[ind_r + ind_seg * dimsPointer[0] + sliceID * (dimsPointer[0]) * (dimsPointer[1])] = tmpMean;
                outputStdDev[ind_r + ind_seg * dimsPointer[0] + sliceID * (dimsPointer[0]) * (dimsPointer[1])] = tmpStdDev;
            }
        }
        delete [] arr;
    }
    copyMask.get();

    hid_t fileOut, dspaceOut, dsetOut;
    fileOut = H5Fopen((globalArgs.outputPath + ".TMP").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hsize_t dimsOut[3] = {dims[0], dimsPointer[1], dimsPointer[0]};
    dspaceOut = H5Screate_simple(3, dimsOut, NULL);
    dsetOut = H5Dcreate2(fileOut, "I_all", H5T_NATIVE_FLOAT, dspaceOut, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dsetOut, H5T_NATIVE_FLOAT, H5S_ALL, dspaceOut, H5P_DEFAULT, outputMean);

    dsetOut = H5Dcreate2(fileOut, "I_std", H5T_NATIVE_FLOAT, dspaceOut, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dsetOut, H5T_NATIVE_FLOAT, H5S_ALL, dspaceOut, H5P_DEFAULT, outputStdDev);
    H5Fclose(fileOut);
    // rename temp file to target filename
    fs::rename((globalArgs.outputPath + ".TMP"), globalArgs.outputPath);
    delete [] outputMean;
    delete [] outputStdDev;
    return 0;
}
