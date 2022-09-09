/**
 * @file xmlParser.cpp
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
#include "xmlParser.h"

using namespace std;
using namespace tinyxml2;
using namespace H5;

static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

xmlParser::xmlParser(){};

template <typename h5GroupData>
void xmlParser::writeAttributes(tinyxml2::XMLElement* element, std::string indent, h5GroupData group) {
    for (tinyxml2::XMLElement* child = element->FirstChildElement("attribute"); child != 0; child = child->NextSiblingElement("attribute")) {
        if (child->ToElement()->Attribute("name")) {
            std::string attrName = child->ToElement()->Attribute("name");
            cxs::io::verbose(3, (indent + "Writing attribute " + attrName + "."));
            if (child->ToElement()->Attribute("source")) {
                // std::cout << ", source = " << child->ToElement()->Attribute("source");

                // std::cout << std::endl;
                try {
                    H5std_string ATTR_NAME(attrName);
                    if (child->ToElement()->Attribute("source")) {
                        std::string source = child->ToElement()->Attribute("source");
                        if (source.size() == 4 && source.compare(0, 4, "spec") == 0) {
                        } else if (source.size() == 8 && source.compare(0, 8, "constant") == 0) {
                            const char* s[1] = {child->ToElement()->Attribute("value")};
                            hsize_t dims1[] = {1};
                            DataSpace sid1(1, dims1);
                            StrType tid1(0, H5T_VARIABLE);
                            Attribute attr = group.createAttribute(ATTR_NAME, tid1, sid1);
                            attr.write(tid1, (void*)s);
                        }
                    }
                } catch (...) {
                    cxs::io::verbose(1, (indent + "Failed to write attribute " + attrName + "."));
                }
            }
        }
    }
}

void xmlParser::parseXML(tinyxml2::XMLElement* element, std::string indent, H5::Group group, std::map<std::string, std::vector<float>>& motor, std::map<std::string, std::vector<float>>& monitor, std::map<std::string, std::vector<std::string>>& metadataContainer, std::map<std::string, datCont_t> datContainer) {
    indent += "\t";
    // write attributes
    writeAttributes(element, indent, group);

    // write datasets
    for (tinyxml2::XMLElement* child = element->FirstChildElement("dataset"); child != 0; child = child->NextSiblingElement("dataset")) {
        if (child->ToElement()->Attribute("name")) {
            std::string dsetName = child->ToElement()->Attribute("name");
            // std::cout << indent + "dataset name = " << dsetName;
            cxs::io::verbose(3, (indent + "Writing dataset " + dsetName + "."));
            H5std_string DATASET_NAME(dsetName);
            const int RANK = 1;

            if (child->ToElement()->Attribute("source")) {
                std::string source = child->ToElement()->Attribute("source");
                if (source.size() == 4 && source.compare(0, 4, "spec") == 0) {
                    try {
                        std::string entry;
                        bool writeData = true;
                        if (child->ToElement()->Attribute("source")) {
                            entry = child->ToElement()->Attribute("entry");
                        } else {
                            throw std::runtime_error("Expected attribute 'entry'.");
                        }
                        cxs::io::verbose(3, (indent + "Loading from spec"));
                        DataSet dataset;
                        std::string dtype = "float";

                        if (motor.find(entry) != motor.end()) {
                            // found entry in motor
                            // std::cout << "Motor Value: " << motor[entry][0] << std::endl;
                            cxs::io::verbose(3, (indent + "Found motor entry in spec for " + dsetName + " (" + entry + ")."));
                            float* data;
                            data = &(motor[entry][0]);
                            int dimsfSz = motor[entry].size();
                            dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, data);
                        } else if (monitor.find(entry) != monitor.end()) {
                            // found entry in monitor
                            cxs::io::verbose(3, (indent + "Found monitor entry in spec for " + dsetName + " (" + entry + ")."));
                            // std::cout << "Monitor Value: " << monitor[entry].size() << std::endl;
                            float* data;
                            data = &(monitor[entry][0]);
                            int dimsfSz = monitor[entry].size();
                            dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, data);
                        } else if (metadataContainer.find(entry) != metadataContainer.end()) {
                            cxs::io::verbose(3, (indent + "Found metadata entry in spec for " + dsetName + " (" + entry + ")."));
                            int dimsfSz;
                            dtype = metadataContainer[entry][1];
                            if (dtype.compare(0, 6, "string") == 0) {
                                dimsfSz = 1;
                                dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, metadataContainer[entry][3]);
                            } else if ((dtype.compare(0, 5, "float") == 0) || ((dtype.compare(0, 3, "int")) == 0)) {
                                float* data;
                                std::vector<float> metadataEntry;
                                uint numEntries = std::stoi(metadataContainer[entry][2]);
                                dimsfSz = numEntries;
                                for (uint ii = 0; ii < numEntries; ii++) {
                                    metadataEntry.push_back(std::stof(metadataContainer[entry][3 + ii]));
                                }
                                data = &(metadataEntry[0]);
                                // std::cout << metadataEntry[0] << std::endl;
                                dataset = writeH5Data(RANK, dimsfSz, dtype, dsetName, group, data);
                            }

                        } else {
                            writeData = false;
                            cxs::io::verbose(1, (indent + "Did not find any entry for " + dsetName + " (" + entry + ")."));
                        }
                        if (writeData) {
                            xmlParser::writeAttributes(child, indent, dataset);
                            if (child->ToElement()->Attribute("units")) {
                                cxs::io::verbose(3, (indent + "Appending units as attribute."));
                                appendUnits(dataset, child->ToElement()->Attribute("units"));
                            }
                        }
                    } catch (...) {
                    }

                } else if (source.size() == 8 && source.compare(0, 8, "constant") == 0) {
                    cxs::io::verbose(3, (indent + "Parsing static value"));
                    try {
                        std::string dtype;
                        if (child->ToElement()->Attribute("type")) {
                            dtype = child->ToElement()->Attribute("type");
                        } else {
                            throw std::runtime_error("Expected attribute 'type'.");
                        }

                        if (dtype.compare(0, 6, "string") == 0) {
                            const char* s[1] = {child->ToElement()->Attribute("value")};
                            hsize_t dims1[] = {1};
                            DataSpace sid1(1, dims1);
                            StrType tid1(0, H5T_VARIABLE);
                            DataSet dataset = group.createDataSet(DATASET_NAME, tid1, sid1);
                            dataset.write((void*)s, tid1);
                            xmlParser::writeAttributes(child, indent, dataset);
                            if (child->ToElement()->Attribute("units")) {
                                cxs::io::verbose(3, (indent + "Appending units as attribute."));
                                appendUnits(dataset, child->ToElement()->Attribute("units"));
                            }
                        } else if (dtype.compare(0, 5, "float") == 0) {
                            float val = std::stof(child->ToElement()->Attribute("value"));
                            float* data = &val;
                            hsize_t dimsf[1];
                            dimsf[0] = 1;
                            DataSpace dataspace(RANK, dimsf);
                            IntType datatype(PredType::NATIVE_FLOAT);
                            datatype.setOrder(H5T_ORDER_LE);
                            DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
                            dataset.write(data, PredType::NATIVE_FLOAT);
                            xmlParser::writeAttributes(child, indent, dataset);
                            if (child->ToElement()->Attribute("units")) {
                                cxs::io::verbose(3, (indent + "Appending units as attribute."));
                                appendUnits(dataset, child->ToElement()->Attribute("units"));
                            }
                        } else {
                            throw std::runtime_error("Expected value of type 'string' or 'float'.");
                        }
                    } catch (...) {
                        // throw std::runtime_error("Failed to parse static value");
                    }
                } else if (source.size() == 8 && source.compare(0, 8, "relative") == 0) {
                    std::string expressionString = child->ToElement()->Attribute("expression");
                    typedef exprtk::symbol_table<float> symbol_table_t;
                    typedef exprtk::expression<float> expression_t;
                    typedef exprtk::parser<float> parser_t;
                    std::string entry;
                    bool useDefault = false;
                    float* data;

                    if (child->ToElement()->Attribute("entry")) {
                        entry = child->ToElement()->Attribute("entry");
                    } else {
                        throw std::runtime_error("Expected attribute 'entry'.");
                    }

                    symbol_table_t symbol_table;
                    if (motor.find(entry) != motor.end()) {
                        symbol_table.add_variable(entry, motor[entry][0]);
                    } else if (monitor.find(entry) != monitor.end()) {
                        symbol_table.add_variable(entry, monitor[entry][0]);
                    } else if (metadataContainer.find(entry) != metadataContainer.end()) {
                        float val = std::stof(metadataContainer[entry][3]);
                        symbol_table.add_variable(entry, val);
                    } else {
                        useDefault = true;
                    }

                    if (!useDefault) {
                        symbol_table.add_constants();

                        expression_t expression;
                        expression.register_symbol_table(symbol_table);

                        parser_t parser;
                        parser.compile(expressionString, expression);
                        float val = expression.value();
                        data = &val;
                    } else {
                        float val = std::stof(child->ToElement()->Attribute("default"));
                        data = &val;
                    }
                    try {
                        hsize_t dimsf[1];
                        dimsf[0] = 1;
                        DataSpace dataspace(RANK, dimsf);
                        IntType datatype(PredType::NATIVE_FLOAT);
                        datatype.setOrder(H5T_ORDER_LE);
                        DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
                        dataset.write(data, PredType::NATIVE_FLOAT);
                        xmlParser::writeAttributes(child, indent, dataset);
                        if (child->ToElement()->Attribute("units")) {
                            cxs::io::verbose(3, (indent + "Appending units as attribute."));
                            appendUnits(dataset, child->ToElement()->Attribute("units"));
                        }
                    } catch (...) {
                    }
                } else if ((source.size() == 3) && (source.compare(0, 3, "dat") == 0)) {
                    std::string datVal;
                    bool writeData = true;
                    if (child->ToElement()->Attribute("source")) {
                        datVal = child->ToElement()->Attribute("entry");
                    } else {
                        throw std::runtime_error("Expected attribute 'entry'.");
                    }
                    DataSet dataset;
                    if (datContainer.find(datVal) != datContainer.end()) {
                        cxs::io::verbose(3, (indent + "Found entry in dat file for " + dsetName + " (" + datVal + ")."));
                        try {
                            cxs::io::verbose(3, ("Writing dataset " + datVal + "."));
                            uint datContRank = datContainer[datVal].dims.size();
                            hsize_t dimsf[datContRank];  // dataset dimensions
                            for (uint rankID = 0; rankID < datContRank; rankID++) {
                                dimsf[rankID] = datContainer[datVal].dims[rankID];
                            }
                            float* data = &(datContainer[datVal].value[0]);
                            // dimsf[1] = NY;
                            DataSpace dataspace(datContRank, dimsf);

                            IntType datatype(PredType::NATIVE_FLOAT);
                            datatype.setOrder(H5T_ORDER_LE);

                            dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);

                            dataset.write(data, PredType::NATIVE_FLOAT);
                        } catch (...) {
                            cxs::io::verbose(3, ("Failed to write dataset " + datVal + "."));
                        }
                    } else {
                        writeData = false;
                        cxs::io::verbose(1, (indent + "Did not find any entry for " + dsetName + " (" + datVal + ")."));
                    }
                    if (writeData) {
                        xmlParser::writeAttributes(child, indent, dataset);
                        if (child->ToElement()->Attribute("units")) {
                            cxs::io::verbose(3, (indent + "Appending units as attribute."));
                            appendUnits(dataset, child->ToElement()->Attribute("units"));
                        }
                    }
                } else {
                    cxs::io::verbose(1, ("Unknown source type " + source));
                }
                // std::cout << ", source = " << child->ToElement()->Attribute("source");
            } else {
                throw std::runtime_error("Expected attribute 'source'.");
            }
        }
    }

    // write groups
    for (tinyxml2::XMLElement* child = element->FirstChildElement("group"); child != 0; child = child->NextSiblingElement("group")) {
        if (child->ToElement()->Attribute("name")) {
            std::string groupName = child->ToElement()->Attribute("name");
            if (child->ToElement()->Attribute("dependency")) {
                bool writeDepData = false;
                std::string depName = child->ToElement()->Attribute("dependency");
                if (H5Lexists(group.getId(), depName.c_str(), H5P_DEFAULT) > 0) {
                    writeDepData = true;
                } else if (depName.compare(0, 7, "entry: ") == 0) {
                    std::vector<std::string> tmpDepNameSplit;
                    tmpDepNameSplit = cxs::utils::split(depName, ' ');
                    depName = tmpDepNameSplit[tmpDepNameSplit.size() - 1];
                    if ((datContainer.find(depName) != datContainer.end()) || (motor.find(depName) != motor.end()) || (monitor.find(depName) != monitor.end())) {
                        writeDepData = true;
                    }
                } else {
                    continue;
                }
                if (writeDepData) {
                    cxs::io::verbose(2, (indent + "Writing group " + groupName + "."));
                    Group subgroup;
                    if (H5Lexists(group.getId(), child->ToElement()->Attribute("name"), H5P_DEFAULT) > 0) {
                        // std::cout << "group exists" << std::endl;
                        Group tmpsubgroup(group.openGroup(child->ToElement()->Attribute("name")));
                        subgroup = tmpsubgroup;
                    } else {
                        Group tmpsubgroup(group.createGroup(child->ToElement()->Attribute("name")));
                        subgroup = tmpsubgroup;
                    }
                    // Group subgroup(group.createGroup(child->ToElement()->Attribute("name")));
                    parseXML(child, indent, subgroup, motor, monitor, metadataContainer, datContainer);
                }
            } else {
                cxs::io::verbose(2, (indent + "Writing group " + groupName + "."));
                Group subgroup;
                if (H5Lexists(group.getId(), child->ToElement()->Attribute("name"), H5P_DEFAULT) > 0) {
                    // std::cout << "group exists" << std::endl;
                    Group tmpsubgroup(group.openGroup(child->ToElement()->Attribute("name")));
                    subgroup = tmpsubgroup;
                } else {
                    Group tmpsubgroup(group.createGroup(child->ToElement()->Attribute("name")));
                    subgroup = tmpsubgroup;
                }
                // Group subgroup(group.createGroup(child->ToElement()->Attribute("name")));
                parseXML(child, indent, subgroup, motor, monitor, metadataContainer, datContainer);
            }
        }

        if (child->ToElement()->Attribute("spec_default")) {
            // spec_default defines the location where we dump the rest, that is any undefined data
            // std::cout << ", spec_default = " << child->ToElement()->Attribute("spec_default");
            tinyxml2::XMLElement* currentChild = child;
            std::string tmpString;
            tmpString = currentChild->ToElement()->Attribute("name");
            specDefaultPath = "/" + tmpString + specDefaultPath;
            collPath = "";
            while (currentChild->Parent()->ToElement()->Attribute("name") != 0) {
                tmpString = currentChild->Parent()->ToElement()->Attribute("name");
                specDefaultPath = "/" + tmpString + specDefaultPath;
                collPath = "/" + tmpString + collPath;
                currentChild = currentChild->Parent()->ToElement();
            }
        }
        // std::cout << std::endl;
    }

    // create hardlink
    for (tinyxml2::XMLElement* child = element->FirstChildElement("hardlink"); child != 0; child = child->NextSiblingElement("hardlink")) {
        std::string linkName;
        if (child->ToElement()->Attribute("name")) {
            linkName = child->ToElement()->Attribute("name");
            cxs::io::verbose(2, (indent + "Creating hardlink for " + linkName + "."));
        }
        std::string targetString;
        if (child->ToElement()->Attribute("target")) {
            targetString = child->ToElement()->Attribute("target");
        } else {
            throw std::runtime_error("Expected attribute 'target'.");
        }

        string linkPath;
        string tmpString;
        bool multiTargets = false;
        std::vector<std::string> targetStringVec;
        
        if (targetString.compare(0,1,"{")==0){
            if (targetString.compare(targetString.size()-1,1,"}")!=0){
                throw std::runtime_error("Expected } delimiter for hardlink entry.");
            }
            multiTargets = true;
            std::string tmpTargetString = targetString.substr(1,targetString.size()-2);
            targetStringVec = cxs::utils::split(tmpTargetString, ',');
        }
        linkPath = "/" + linkName;
        tinyxml2::XMLElement* currentChild = child;
        while (currentChild->Parent()->ToElement()->Attribute("name") != 0) {
            tmpString = currentChild->Parent()->ToElement()->Attribute("name");
            linkPath = "/" + tmpString + linkPath;
            currentChild = currentChild->Parent()->ToElement();
        }
        if (multiTargets) {
            uint entryID = 0;
            while (entryID < targetStringVec.size()) {
                try {
                    ltrim(targetStringVec[entryID]);
                    group.link(H5L_TYPE_HARD, targetStringVec[entryID], linkPath);
                    break;
                } catch (...) {
                    cxs::io::verbose(1, (indent + "Failed to create hardlink " + linkName + " for " + targetStringVec[entryID] + "."));
                }
                entryID++;
            }
        } else {
            try {
                group.link(H5L_TYPE_HARD, targetString, linkPath);
            } catch (...) {
                cxs::io::verbose(1, (indent + "Failed to create hardlink " + linkName + " for " + targetString + "."));
            }
        }
    }
    // create softlink
    for (tinyxml2::XMLElement* child = element->FirstChildElement("softlink"); child != 0; child = child->NextSiblingElement("softlink")) {
        std::string linkName;
        if (child->ToElement()->Attribute("name")) {
            linkName = child->ToElement()->Attribute("name");
            cxs::io::verbose(2, (indent + "Creating softlink for " + linkName + "."));
        }
        std::string targetString;
        if (child->ToElement()->Attribute("target")) {
            targetString = child->ToElement()->Attribute("target");
        } else {
            throw std::runtime_error("Expected attribute 'target'.");
        }

        string linkPath;
        string tmpString;
        bool multiTargets = false;
        std::vector<std::string> targetStringVec;
        if (targetString.compare(0,1,"{")==0){
            if (targetString.compare(targetString.size()-1,1,"}")!=0){
                throw std::runtime_error("Expected } delimiter for softlink entry.");
            }
            multiTargets = true;
            std::string tmpTargetString = targetString.substr(1,targetString.size()-2);
            targetStringVec = cxs::utils::split(tmpTargetString, ',');
        }
        linkPath = "/" + linkName;
        tinyxml2::XMLElement* currentChild = child;
        while (currentChild->Parent()->ToElement()->Attribute("name") != 0) {
            tmpString = currentChild->Parent()->ToElement()->Attribute("name");
            linkPath = "/" + tmpString + linkPath;
            currentChild = currentChild->Parent()->ToElement();
        }
        if (multiTargets) {
            uint entryID = 0;
            while (entryID < targetStringVec.size()) {
                try {
                    ltrim(targetStringVec[entryID]);
                    group.link(H5L_TYPE_SOFT, targetStringVec[entryID], linkPath);
                    break;
                } catch (...) {
                    cxs::io::verbose(1, (indent + "Failed to create softlink " + linkName + " for " + targetStringVec[entryID] + "."));
                }
                entryID++;
            }
        } else {
            try {
                group.link(H5L_TYPE_SOFT, targetString, linkPath);
            } catch (...) {
                cxs::io::verbose(1, (indent + "Failed to create softlink " + linkName + " for " + targetString + "."));
            }
        }
    }
}
