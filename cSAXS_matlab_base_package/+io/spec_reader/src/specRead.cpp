/**
 * @file specRead.cpp
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
#include "xmlParser.h"

using namespace std;
using namespace tinyxml2;
using namespace H5;
using json = nlohmann::json;

extern int cxs::io::GLOBAL_verbosity;

specRead::specRead(int& scanNr, std::ifstream& datFile) {
    scan = scanNr;
    scanNr_str = std::to_string(scanNr);

    startString = "#S " + scanNr_str + " ";
    endString = "#X " + scanNr_str;
    startScanString = "#L";
    motorString = "#O";
}

void specRead::flush() {
    cxs::io::verbose(2, "Flushing storage.");
    // private variables
    scanContent.clear();
    scanContentName.clear();
    metadata.clear();
    // motorConfig is cleared within read_from_file

    // public variables
    motor.clear();
    monitor.clear();
    header.clear();
    metadataContainer.clear();
    datContainer.clear();


}

// read content and save scan for processing
void specRead::read_from_file(const std::vector<std::string>& inputFile) {
    int start = 0;
    lineidx = 0;
    scanexists = false;
    // #pragma omp critical
    // cout << inputFile.size() << endl;
    while (lineidx < inputFile.size()) {
        // line = inputFile[lineidx];

        lineidx++;

        //look for motor configs
        if (inputFile[lineidx].compare(0, motorString.size(), motorString) == 0) {
            if (inputFile[lineidx].compare(0, 3, "#O0") == 0) {
                motorConfig.clear();
                motorConfig.push_back(inputFile[lineidx]);
                // cout << motorlineStart << endl;
            } else {
                motorConfig.push_back(inputFile[lineidx]);
            }
        }

        if (start > 0 && ((inputFile[lineidx].compare(0, endString.size(), endString) == 0) || (inputFile[lineidx].compare(0, 2, "#C") == 0 && (inputFile[lineidx].find("Scan aborted") != std::string::npos))|| (inputFile[lineidx].compare(0, 3, "#S ") == 0 ))) {
            // cout << "line" << endl;

            prepareMonitor(scanContentName, scanContent, motor, monitor);
            start = 0;
            scanexists = true;

        }

        if ((start > 0 && inputFile[lineidx].compare(0, 2, "#C") == 0)) {
            if (inputFile[lineidx].compare(0, 7, "#C meta") == 0) {
                std::vector<std::string> tmpMetaContainer;
                istringstream iss(inputFile[lineidx]);
                tmpMetaContainer = cxs::utils::split(inputFile[lineidx], ' ');
                iss.str(std::string());
                iss.clear();
                std::vector<std::string> vec(tmpMetaContainer.cbegin() + 2, tmpMetaContainer.cbegin() + tmpMetaContainer.size());
                // std::cout << vec[3] << std::endl;
                metadataContainer[tmpMetaContainer[2]] = vec;
            }
        }

        if (start == 2) {
            istringstream iss(inputFile[lineidx]);
            scanContent.push_back(cxs::utils::split(inputFile[lineidx], ' '));
        }

        if (start == 0 && inputFile[lineidx].compare(0, startString.size(), startString) == 0) {
            // cout << inputFile[lineidx] << endl;
            if (scanexists==true){
                specRead::flush();
            }
            start = 1;
            //  cout << "Motor config: " << motorlineStart << " - " << motorlineStop << endl;
        }

        if (start == 1) {
            if (inputFile[lineidx].compare(0, startScanString.size(), startScanString) == 0) {
                start = 2;
                istringstream iss(inputFile[lineidx]);
                scanContentName = cxs::utils::split(inputFile[lineidx], ' ');
                iss.str(std::string());
                iss.clear();

                // prepare motor dictionary and initialize with #P values - could run asynchronously
                motor = prepareMotorConfig(motorConfig, motorUpdate);
                // motorRes =  std::async(std::launch::async | std::launch::deferred, prepareMotorConfig,motorConfig,motorUpdate);

            } else {
                // #S - scan string
                if (inputFile[lineidx].compare(0, startString.size(), startString) == 0) {
                    header.push_back(inputFile[lineidx]);
                    std::vector<std::string> metaInterim; 
                    // metaInterim.push_back("#S");
                    // metaInterim.push_back("meta");
                    metaInterim.push_back("S");
                    metaInterim.push_back("string");
                    metaInterim.push_back("1");
                    metaInterim.push_back(inputFile[lineidx]);
                    metadataContainer["S"] = metaInterim;
                }
                // #D - date
                else if (inputFile[lineidx].compare(0, 2, "#D") == 0) {
                    std::string datetime = getDateTime(inputFile[lineidx]);
                    header.push_back(datetime);
                    std::vector<std::string> metaInterim; 
                    // metaInterim.push_back("#D");
                    // metaInterim.push_back("meta");
                    metaInterim.push_back("D");
                    metaInterim.push_back("string");
                    metaInterim.push_back("1");
                    metaInterim.push_back(datetime);
                    metadataContainer["D"] = metaInterim;
                } else if (inputFile[lineidx].compare(0, 2, "#P") == 0) {
                    motorUpdate.push_back(inputFile[lineidx]);
                }
            }
        }
    }
    if ((lineidx == inputFile.size()) && (start > 0)){

        prepareMonitor(scanContentName, scanContent, motor, monitor);
        start = 0;
        scanexists = true;

    }
}

void specRead::read_from_orchestra(const std::string& orchestra_path, bool fullPath, std::string groupName) {
    if (orchestra_path.length() == 0) {
        return;
    }
    std::string posFile;
    std::vector<std::vector<std::string>> orchestraContent;
    std::vector<std::string> orchestraLogName;
    int orchestraLineidx = 0;

    if (fullPath) {
        posFile = orchestra_path;
    } else {
        int numlength = scanNr_str.length();
        posFile = orchestra_path + "scan_" + std::string(std::max(5 - numlength, 0), '0') + scanNr_str + ".dat";
    }
    
    cxs::io::verbose(2, ("Reading orchestra file" + posFile));
    std::ifstream orchestraFile(posFile);
    if (orchestraFile.is_open()) {
        while (orchestraFile.peek() != EOF) {
            getline(orchestraFile, line);
            orchestraLineidx++;
            if (orchestraLineidx < 2) {
                continue;
            } else if (orchestraLineidx == 2) {
                istringstream iss(line);
                orchestraLogName = cxs::utils::split(line, ' ');
                // for (auto LogName: orchestraLogName){
                //     std::cout << LogName << std::endl;
                // }
                iss.str(std::string());
                iss.clear();
            } else {
                istringstream iss(line);
                orchestraContent.push_back(cxs::utils::split(line, ' '));
                iss.str(std::string());
                iss.clear();
            }
        }
        if (groupName.length()==0){
            groupName = "orchestra";
        }

        prepareDatContainer(orchestraLogName, orchestraContent, datContainer, groupName);
    } else {
        cxs::io::verbose(1, "Failed to open orchestra file " + posFile);
    }
}

void specRead::read_datFiles(const std::string& datFiles) {
    if (datFiles.length() == 0) {
        return;
    }

    json j = json::parse(datFiles);
    std::string keyVal, datFilePath;
    uint nops;
    if (datFiles.compare(0, 1, "[") == 0) {
        nops = j.size();
    } else {
        nops = 1;
    }

    for (uint inops = 0; inops < nops; inops++) {
        std::string path, type, entry, groupName;
        json jsub;

        if (nops > 1) {
            jsub = j[inops];
        } else {
            jsub = j;
        }

        for (json::iterator it = jsub.begin(); it != jsub.end(); ++it) {
            keyVal = it.key();
            if (keyVal.compare("path") == 0) {
                path = (std::string)it.value();
            } else if (keyVal.compare("type") == 0) {
                type = (std::string)it.value();
            } else if (keyVal.compare("entry") == 0) {
                entry = (std::string)it.value();
            } else if (keyVal.compare("groupName") == 0) {
                groupName = (std::string)it.value();
            } else {
                throw std::runtime_error("Unexpected json keyword: " + keyVal);
            }
            if (type.compare("orchestra") == 0) {
                specRead::read_from_orchestra(path, true, groupName);
            } else if (type.compare("mcs") == 0) {
                if (entry.length()==0){
                    throw std::runtime_error("Missing keyword 'entry' for datFile of type mcs.");
                }
                std::ifstream datFile(path);
                std::vector<std::vector<std::string>> datFileContent;
                if (datFile.is_open()) {
                    while (datFile.peek() != EOF) {
                        getline(datFile, line);
                        if (line.compare(0, 1, "#") == 0) {
                            continue;
                        } else {
                            istringstream iss(line);
                            std::vector<std::string> tmpDatFileContent;
                            tmpDatFileContent = cxs::utils::split(line, ' ');
                            // std::cout << tmpDatFileContent.size() << std::endl;
                            if (tmpDatFileContent.size() > 1) {
                                datFileContent.push_back(tmpDatFileContent);
                            }
                            iss.str(std::string());
                            iss.clear();
                        }
                    }
                    if (datFileContent.size() > 0) {
                        prepareDatContainer(datFileContent, datContainer, type, entry);
                    }
                } else {
                    cxs::io::verbose(1, "Failed to open datFile " + path);
                }
            }
        }
    }

}

/**
 * @brief output stream, either json or HDF5
 * 
 * @param hdf5 
 * @param output 
 * @param scanNrIdx 
 * @param scanNrSz 
 * @param xmlLayoutFile 
 */
void specRead::output_stream(const bool hdf5, const std::string output, const uint& scanNrIdx, const uint scanNrSz, const std::string& xmlLayoutFile) {
    bool dump = true;
    if (xmlLayoutFile.length() > 0) {
        dump = false;
    }

    if (!scanexists) {
        // std::cout << "Failed to read scan" + scanNr_str << std::endl;
        return;
    }
    if (hdf5) {
        // convert to hdf5
        if (cxs::io::GLOBAL_verbosity < 2) {
            Exception::dontPrint();
        }
        std::string fname;
        if (output.length() != 0) {
            fname = output;
        } else {
            fname = "metadata_S" + std::string(std::max(5 - (int)scanNrIdx, 0), '0') + ".h5";
        }
        const H5std_string FILE_NAME(fname);
        const int RANK = 1;
        H5File file;
        try {
            H5File tmpfile(FILE_NAME, H5F_ACC_RDWR);
            file = tmpfile;
        } catch (...) {
            H5File tmpfile(FILE_NAME, H5F_ACC_TRUNC);
            file = tmpfile;
        }

        if (dump) {
            // H5File file( FILE_NAME, H5F_ACC_RDWR );

            // std::cout << monitor.size() << std::endl;
            // std::cout << motor.size() << std::endl;
            for (std::map<std::string, std::vector<float>>::iterator it = monitor.begin(); it != monitor.end(); ++it) {
                H5std_string DATASET_NAME(it->first);
                // int i, j;
                // int data[NX][NY];          // buffer for data to write
                cxs::io::verbose(3, ("Writing dataset " + it->first + "."));
                hsize_t dimsf[1];  // dataset dimensions
                dimsf[0] = it->second.size();
                float* data = &(it->second[0]);
                // dimsf[1] = NY;
                DataSpace dataspace(RANK, dimsf);

                IntType datatype(PredType::NATIVE_FLOAT);
                datatype.setOrder(H5T_ORDER_LE);

                DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);

                dataset.write(data, PredType::NATIVE_FLOAT);
                // monitor.erase(it->first);
            }
            for (std::map<std::string, std::vector<float>>::iterator it = motor.begin(); it != motor.end(); ++it) {
                H5std_string DATASET_NAME(it->first);
                cxs::io::verbose(3, ("Writing dataset " + it->first + "."));
                hsize_t dimsf[1];  // dataset dimensions
                dimsf[0] = it->second.size();
                float* data = &(it->second[0]);
                // dimsf[1] = NY;
                DataSpace dataspace(RANK, dimsf);

                IntType datatype(PredType::NATIVE_FLOAT);
                datatype.setOrder(H5T_ORDER_LE);

                DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);

                dataset.write(data, PredType::NATIVE_FLOAT);
                // motor.erase(it->first);
            }

            for (std::map<std::string, std::vector<std::string>>::iterator it = metadataContainer.begin(); it != metadataContainer.end(); ++it) {
                cxs::io::verbose(3, ("Writing dataset " + it->first + "."));
                // // cxs::io::verbose(3, (indent + "Found metadata entry in spec for " + dsetName + " (" + specVal + ")."));
                int dimsfSz;
                std::string dtype = it->second[1];
                DataSet dataset;
                std::string dsetName = it->first;
                if (dtype.compare(0, 6, "string") == 0) {
                    dimsfSz = 1;
                    dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, file, it->second[3]);
                } else if ((dtype.compare(0, 5, "float") == 0) || ((dtype.compare(0, 3, "int")) == 0)) {
                    float* data;
                    std::vector<float> metadataEntry;
                    uint numEntries = std::stoi(it->second[2]);
                    dimsfSz = numEntries;
                    for (uint ii = 0; ii < numEntries; ii++) {
                        metadataEntry.push_back(std::stof(it->second[3 + ii]));
                    }
                    data = &(metadataEntry[0]);
                    // std::cout << metadataEntry[0] << std::endl;
                    dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, file, data);
                }
            }
        } else {
            // parse XML file and adjust H5 tree structure to given nexus format
            xmlParser p = xmlParser();

            tinyxml2::XMLDocument doc;
            tinyxml2::XMLError eResult = doc.LoadFile(xmlLayoutFile.c_str());
            if (eResult != tinyxml2::XML_SUCCESS) {
                throw std::runtime_error("Failed to open XML layout file.");
            };
            // std::cout << "opened file" << std::endl;
            // hid_t H5::H5File::getId	(		)	const
            Group group;
            if (H5Lexists(file.getId(), "/entry", H5P_DEFAULT) > 0) {
                // std::cout << "group exists" << std::endl;
                Group tmpGroup(file.openGroup("/entry"));
                group = tmpGroup;
            } else {
                Group tmpGroup(file.createGroup("/entry"));
                group = tmpGroup;
            }
            // std::cout << group.getId() << std::endl;
            // loop over groups
            p.parseXML(doc.FirstChildElement("hdf5_layout")->FirstChildElement("group"), "", group, motor, monitor, metadataContainer, datContainer);

            // dump collection
            if (p.specDefaultPath.size() == 0) {
                throw std::runtime_error("Expected spec_default attribute in XML layout file.");
            }
            Group collGroup(file.openGroup(p.specDefaultPath));
            // std::cout << monitor.size() << std::endl;
            // std::cout << motor.size() << std::endl;
            for (std::map<std::string, std::vector<float>>::iterator it = monitor.begin(); it != monitor.end(); ++it) {
                try {
                    H5std_string DATASET_NAME(it->first);
                    cxs::io::verbose(3, ("Writing dataset " + it->first + "."));
                    hsize_t dimsf[1];  // dataset dimensions
                    dimsf[0] = it->second.size();
                    float* data = &(it->second[0]);
                    // dimsf[1] = NY;
                    DataSpace dataspace(RANK, dimsf);

                    IntType datatype(PredType::NATIVE_FLOAT);
                    datatype.setOrder(H5T_ORDER_LE);

                    DataSet dataset = collGroup.createDataSet(DATASET_NAME, datatype, dataspace);

                    dataset.write(data, PredType::NATIVE_FLOAT);
                } catch (...) {
                    cxs::io::verbose(3, ("Failed to write dataset " + it->first + "."));
                }
            }
            for (std::map<std::string, std::vector<float>>::iterator it = motor.begin(); it != motor.end(); ++it) {
                try {
                    H5std_string DATASET_NAME(it->first);
                    cxs::io::verbose(3, ("Writing dataset " + it->first + "."));

                    hsize_t dimsf[1];  // dataset dimensions
                    dimsf[0] = it->second.size();
                    float* data = &(it->second[0]);
                    // dimsf[1] = NY;
                    DataSpace dataspace(RANK, dimsf);

                    IntType datatype(PredType::NATIVE_FLOAT);
                    datatype.setOrder(H5T_ORDER_LE);

                    DataSet dataset = collGroup.createDataSet(DATASET_NAME, datatype, dataspace);

                    dataset.write(data, PredType::NATIVE_FLOAT);
                } catch (...) {
                    cxs::io::verbose(3, ("Failed to write dataset " + it->first + "."));
                }
            }

            for (std::map<std::string, std::vector<std::string>>::iterator it = metadataContainer.begin(); it != metadataContainer.end(); ++it) {
                cxs::io::verbose(3, ("Writing dataset " + it->first + "."));
                // cxs::io::verbose(3, (indent + "Found metadata entry in spec for " + dsetName + " (" + specVal + ")."));
                try {
                    int dimsfSz;
                    std::string dtype = it->second[1];
                    DataSet dataset;
                    std::string dsetName = it->first;
                    if (dtype.compare(0, 6, "string") == 0) {
                        dimsfSz = 1;
                        dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, collGroup, it->second[3]);
                    } else if ((dtype.compare(0, 5, "float") == 0) || ((dtype.compare(0, 3, "int")) == 0)) {
                        float* data;
                        std::vector<float> metadataEntry;
                        uint numEntries = std::stoi(it->second[2]);
                        dimsfSz = numEntries;
                        for (uint ii = 0; ii < numEntries; ii++) {
                            metadataEntry.push_back(std::stof(it->second[3 + ii]));
                        }
                        data = &(metadataEntry[0]);
                        // std::cout << metadataEntry[0] << std::endl;
                        dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, collGroup, data);
                    }
                } catch (...) {
                    cxs::io::verbose(3, ("Failed to write dataset " + it->first + "."));
                }
            }
            for (std::map<std::string, datCont_t>::iterator it = datContainer.begin(); it != datContainer.end(); ++it) {
                Group datGroup;
                std::string datGroupPath;
                datGroupPath = p.collPath + "/" + it->second.source;
                if (H5Lexists(file.getId(), datGroupPath.c_str(), H5P_DEFAULT) > 0) {
                    // std::cout << "group exists" << std::endl;
                    Group tmpGroup(file.openGroup(datGroupPath.c_str()));
                    datGroup = tmpGroup;
                } else {
                    Group tmpGroup(file.createGroup(datGroupPath.c_str()));
                    datGroup = tmpGroup;
                }
                try {
                    H5std_string DATASET_NAME(it->first);
                    cxs::io::verbose(3, ("Writing dataset " + it->first + "."));
                    uint datContRank = it->second.dims.size();
                    hsize_t dimsf[datContRank];  // dataset dimensions
                    for (uint rankID=0; rankID<datContRank; rankID++){
                        dimsf[rankID] = it->second.dims[rankID];
                    }
                    float* data = &(it->second.value[0]);
                    // dimsf[1] = NY;
                    DataSpace dataspace(datContRank, dimsf);

                    IntType datatype(PredType::NATIVE_FLOAT);
                    datatype.setOrder(H5T_ORDER_LE);

                    DataSet dataset = datGroup.createDataSet(DATASET_NAME, datatype, dataspace);

                    dataset.write(data, PredType::NATIVE_FLOAT);
                } catch (...) {
                    cxs::io::verbose(3, ("Failed to write dataset " + it->first + "."));
                }
            }
        }

    } else {
        // prepare json stream

        if (scanNrSz > 1 && scanNrIdx == 0) {
            jsonOut = "[{";//\"S\":\"" + header[0] + "\", \n";
        } else {
            jsonOut = "{";//\"S\":\"" + header[0] + "\", \n";
        }
        // jsonOut += "\"D\":\"" + header[1] + "\", \n";

        for (std::map<std::string, std::vector<float>>::iterator it = motor.begin(); it != motor.end(); ++it) {
            // std::cout << it-> second.size() << std::endl;

            if (it->second.size() == 1) {
                jsonOut += "\"" + it->first + "\":" + std::to_string(it->second[0]) + ", \n";
            } else if (it->second.size() == 0) {
                    continue;
            } else {
                jsonOut += "\"" + it->first + "\":[";
                for (uint arrID = 0; arrID < it->second.size() - 1; arrID++) {
                    jsonOut += std::to_string(it->second[arrID]) + ", ";
                }
                jsonOut += std::to_string(it->second[it->second.size() - 1]) + "], \n";
            }
        }

        for (std::map<std::string, std::vector<float>>::iterator it = monitor.begin(); it != monitor.end(); ++it) {
            if (it->second.size() == 1) {
                jsonOut += "\"" + it->first + "\":" + std::to_string(it->second[0]) + ", \n";
            } else if (it->second.size() == 0) {
                continue;
            } else {
                jsonOut += "\"" + it->first + "\":[";
                for (uint arrID = 0; arrID < it->second.size() - 1; arrID++) {
                    jsonOut += std::to_string(it->second[arrID]) + ", ";
                }
                jsonOut += std::to_string(it->second[it->second.size() - 1]) + "], \n";
            }
        }

        for (std::map<std::string, std::vector<std::string>>::iterator it = metadataContainer.begin(); it != metadataContainer.end(); ++it) {
            std::string elementPrefix;
            std::string elementSuffix;

            if (it->second[1].compare(0, it->second[1].size(), "string") == 0) {
                elementPrefix = elementSuffix = "\"";
            } else {
                elementPrefix = elementSuffix = "";
            }

            uint numEntries = std::stoi(it->second[2]);
            if (numEntries == 1) {
                jsonOut += "\"" + it->first + "\":" + elementPrefix + it->second[3] + elementSuffix + ", \n";
            } else {
                jsonOut += "\"" + it->first + "\":[";
                for (uint arrID = 0; arrID < numEntries - 1; arrID++) {
                    jsonOut += elementPrefix + it->second[arrID + 3] + elementSuffix + ", ";
                }
                jsonOut += elementPrefix + it->second[numEntries + 2] + elementSuffix + "], \n";
            }
        }
        if (scanNrSz > 1 && (scanNrIdx == scanNrSz - 1)) {
            jsonOut.replace(jsonOut.length() - 3, jsonOut.length(), "}]\n");
        } else if (scanNrSz == 1 && (scanNrIdx == scanNrSz - 1)) {
            jsonOut.replace(jsonOut.length() - 3, jsonOut.length(), "}\n");
        } else {
            jsonOut.replace(jsonOut.length() - 3, jsonOut.length(), "}, \n");
        }

        if (output.length() != 0) {
            // send to file
            ofstream jsonFile(output + "scan_" + scanNr_str + ".json");
            if (jsonFile.is_open()) {
                jsonFile << jsonOut;
                jsonFile.close();
            } else {
                std::cout << "Failed to open json File!" << std::endl;
            }
        }
    }

}  // end output_stream
