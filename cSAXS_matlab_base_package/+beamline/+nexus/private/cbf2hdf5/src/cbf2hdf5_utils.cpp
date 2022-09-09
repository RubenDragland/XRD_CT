/**
 * @file cbf2hdf5_utils.cpp
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

#include "cbf2hdf5_utils.h"
#define DEFLATE_SIZE_ADJUST(s) (ceil(((double)(s)) * 1.001) + 12)

// ----------------------------------------------------------------

/* PILATUS CBF DATA
        _array_data.data
        ;
        --CIF-BINARY-FORMAT-SECTION--
        Content-Type: application/octet-stream;
            conversions="x-CBF_BYTE_OFFSET"
        Content-Transfer-Encoding: BINARY
        X-Binary-Size: 2499331
        X-Binary-ID: 1
        X-Binary-Element-Type: "signed 32-bit integer"
        X-Binary-Element-Byte-Order: LITTLE_ENDIAN
        Content-MD5: XoY7+gfct1+OKiJlzKOiRw==
        X-Binary-Number-of-Elements: 2476525
        X-Binary-Size-Fastest-Dimension: 1475
        X-Binary-Size-Second-Dimension: 1679
        X-Binary-Size-Padding: 4095
     */
std::regex pilatus_ncols_regex(R"(X-Binary-Size-Fastest-Dimension: (\d+))");  //!< Regex for pilatus ncols
std::regex pilatus_nrows_regex(R"(X-Binary-Size-Second-Dimension: (\d+))");   //!< Regex for pilatus nrows

/*!
     * \brief Read CBF image dimensions
     *
     * \param data (IN) CBF image file data
     * \param dims (OUT) to be filled with [nrows, ncols] from the data
     */
void cbf_dims(const std::vector<char> &data, std::vector<long> &dims) {
    assert(dims.size() == 2);
    {
        std::cmatch match;
        if (!std::regex_search(&data[0], &data[0] + data.size(), match, pilatus_ncols_regex))
            throw std::runtime_error("unable to find number of columns");
        dims[1] = std::stol(match[1].str());
        if (dims[1] <= 0)
            throw std::runtime_error("dimension along row is not positive");
    }
    {
        std::cmatch match;
        if (!std::regex_search(&data[0], &data[0] + data.size(), match, pilatus_nrows_regex))
            throw std::runtime_error("unable to find number of rows");
        dims[0] = std::stol(match[1].str());
        if (dims[0] <= 0)
            throw std::runtime_error("dimension along column is not positive");
    }
}

/**
 * @brief read .cbf file, compress it and write it to disk
 * 
 * @param path          full path of the cbf file
 * @param fileID        slice of the resulting 3D dataset
 * @param dataset       h5 dataset handle
 * @param datatype      h5 datatype handle
 * @param dataspace     h5 dataspace handle
 * @param dim           dimensions of the 3D dataset (row columns)
 * @param aggression    compression level
 */
void pilatus_read_data(std::string &path, uint fileID, hid_t dataset, hid_t datatype, hid_t dataspace, hsize_t * dim, uint aggression, uint64_t *timeStamp) {

    // Hyperslab parameters for dataset
    hsize_t start[3] = {0, 0, 0};
    hsize_t h5count[3] = {1, dim[1], dim[2]};
    hsize_t cdims[3] = {1, dim[1], dim[2]};
    start[0] = fileID;
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);

    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, h5count, NULL);

    // Use Heiners method to read the data
    // Adapt dimensions for every image file
    std::vector<char> fbuf;
    {
        FILE *fin = fopen(path.c_str(), "r");
        if (!fin)
            throw std::runtime_error(std::string("unable to open file ") + path + ": " + std::strerror(errno));
        try {
            off_t fsz;
            {
                struct stat sbuf;
                if (fstat(fileno(fin), &sbuf) == -1) {
                    throw std::runtime_error(std::string("unable to stat file ") + path + ": " + std::strerror(errno));
                }
                fsz = sbuf.st_size;
            }
            fbuf.resize(fsz);
            fread(&fbuf[0], 1, fsz, fin);
            if (ferror(fin))
                throw std::runtime_error(std::string("unable to read file ") + path + ": " + std::strerror(errno));
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
    std::string det = "Detector:";
    std::string detPound = "# ";
    auto ts1 = std::search(fbuf.begin(), fbuf.end(), det.begin(), det.end());
    auto ts2 = std::search(ts1, fbuf.end(), detPound.begin(), detPound.end());

    int offset = ts2-fbuf.begin()+2;

    std::string s(fbuf.begin()+offset, fbuf.begin()+offset+23);
    s.erase(s.begin() + 4);
    s.erase(s.begin() + 6);
    s.erase(s.begin() + 8);
    s.erase(s.begin() + 10);
    s.erase(s.begin() + 12);
    s.erase(s.begin() + 14);
    timeStamp[fileID] = std::stoul(s);

    unsigned long nelems = fdims[0] * fdims[1];
    int *data = new int[nelems]; //std::vector<float> data(nelems);
    int current = 0;
    for (unsigned int i = 0; i < nelems; i++) {
        if (*((uint8_t *)&fbuf[finger]) != 0x80) {  // | xx |
            current += *((int8_t *)&fbuf[finger]);
            finger += 1;
        } else if (*((uint16_t *)&fbuf[finger + 1]) != 0x8000) {  // | 80 | xx | xx |
            current += *((int16_t *)&fbuf[finger + 1]);
            finger += 3;
        } else {  // | 80 | 80 | 00 | xx | xx | xx | xx |
            current += *((int32_t *)&fbuf[finger + 3]);
            finger += 7;
        }
        if (finger + 7 > fbuf.size())
            throw std::runtime_error(std::string("data inconsistency in file ") + path);
        //if (current < -1)   // allow value -1, which is used to mark detector gaps
        //    throw std::runtime_error(std::string("data error in file ") + path);
        data[i] = current;
        // std::cout << current << std::endl;
    }
    fbuf.clear();
    // std::vector<long> fstride{0, 0};
    // std::vector<long> mstride{0, 0};
    // std::vector<long> count{0, 0};
    {
        unsigned filter_mask = 0;
        size_t buf_size = cdims[1] * cdims[2] * sizeof(unsigned int);
        Bytef *z_dst; /*destination buffer            */
        uLongf z_dst_nbytes = (uLongf)(DEFLATE_SIZE_ADJUST(buf_size));
        uLong z_src_nbytes = (uLong)buf_size;

        unsigned int *outbuf[1];
        outbuf[0] = (unsigned int *)malloc((size_t)z_dst_nbytes);
        z_dst = (Bytef *)outbuf[0];

        const Bytef *z_src = (const Bytef *)data;
        compress2(z_dst, &z_dst_nbytes, z_src, z_src_nbytes, aggression);

        H5DOwrite_chunk(dataset, dxpl,
                        filter_mask, start, (size_t)z_dst_nbytes, outbuf[0]);

        //H5Dwrite(dataset, datatype , dataspaceimg, dataspace, H5P_DEFAULT,
        // (longedge_x ? map: mapr));

        free(outbuf[0]);
    }
    delete [] data;

}

/**
 * @brief append units as attribute
 * 
 * @param gid       h5 group handle
 * @param units     std::string containing the unit
 */
void appendUnits(hid_t gid, const std::string& units) {
    const char* s[1] = {units.c_str()};
    hid_t filetype, memtype, space, attr;
    filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, H5T_VARIABLE);
    memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, H5T_VARIABLE);
    space = H5Screate(H5S_SCALAR);
    attr = H5Acreate(gid, "units", filetype, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, memtype, (void*)s);
    H5Sclose(space);
    H5Aclose(attr);
}