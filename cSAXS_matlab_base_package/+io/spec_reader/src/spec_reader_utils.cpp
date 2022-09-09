/**
 * @file spec_reader_utils.cpp
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
#include "spec_reader_utils.h"

using namespace std;
using namespace H5;

extern int cxs::io::GLOBAL_verbosity;



std::map<std::string, vector<float>> prepareMotorConfig(const std::vector <std::string>&motorLog, const std::vector <std::string>&motorVal){
    std::map<std::string, vector<float>> motor;
    std::vector <std::string> motorConfigList;
    std::vector <std::string> motorValList;
    // std::cout << motorLog.size() << std::endl;
    for (uint vecID=0; vecID<motorLog.size(); vecID++){
        // std::cout << motorLog[vecID] << std::endl;
        istringstream iss(motorLog[vecID]);
        motorConfigList = cxs::utils::split(motorLog[vecID], ' ');
        iss.str(std::string());
        iss.clear();
        istringstream issVal(motorVal[vecID]);
        motorValList = cxs::utils::split(motorVal[vecID], ' ');
        issVal.str(std::string());
        issVal.clear();

        for (uint motorID=1; motorID<motorConfigList.size(); motorID++){
            motor[motorConfigList[motorID]].push_back(std::stof(motorValList[motorID]));
        }
    }

    return motor;
}

void prepareMonitor(const std::vector<std::string>& monitorLogName, const std::vector <std::vector <std::string>>& monitorLog, std::map<std::string, std::vector<float>>& motor, std::map<std::string, vector<float>>& monitor){
    std::vector <float> monitorEntry;

    for (uint nameID=1; nameID<monitorLogName.size(); nameID++){
        // cout << monitorLogName.size() << endl;

        monitorEntry.clear();
        // monitorEntry.reserve(monitorLog.size());
        for (uint vecID=0; vecID<monitorLog.size(); vecID++){
            // monitorEntryArray[vecID] = std::stof(monitorLog[vecID][nameID-1]);
            if (monitorLog[vecID].size() == monitorLogName.size()-1){
            try{
                // std::cout << monitorLog[vecID][nameID-1] << std::endl;
                monitorEntry.push_back(std::stof(monitorLog[vecID][nameID-1]));
            } catch(...){
                throw std::runtime_error("Failed to convert monitorLog");
            }
            }

            
        }

        if ( motor.find(monitorLogName[nameID]) == motor.end()){
            monitor[monitorLogName[nameID]] = monitorEntry;
        }
        else {

            motor[monitorLogName[nameID]] = monitorEntry;
        }


    }

return;
}

void prepareDatContainer(const std::vector<std::string>& monitorLogName, const std::vector <std::vector <std::string>>& monitorLog, std::map<std::string, datCont_t>& datContainer, std::string source){
    std::vector <float> monitorEntry;

    datCont_t tmpContainer;
    tmpContainer.source = source;
    tmpContainer.dims.push_back(monitorLog.size());
    // std::cout << monitorLogName[0] << std::endl;
    for (uint nameID=0; nameID<monitorLogName.size(); nameID++){

        monitorEntry.clear();
        // monitorEntry.reserve(monitorLog.size());
        // std::cout << nameID << std::endl;
        for (uint vecID=0; vecID<monitorLog.size(); vecID++){
            // monitorEntryArray[vecID] = std::stof(monitorLog[vecID][nameID-1]);
            try{
                monitorEntry.push_back(std::stof(monitorLog[vecID][nameID]));
                std::cout << std::stof(monitorLog[vecID][nameID]) << std::endl;
            } catch(...){
                throw std::runtime_error("Failed to convert to datContainer.");
            }
            
        }
        // std::cout << monitorLogName[nameID] << std::endl;
        tmpContainer.value = monitorEntry;
        
        datContainer[monitorLogName[nameID]] = tmpContainer;


    }

return;
}

void prepareDatContainer(const std::vector <std::vector <std::string>>& monitorLog, std::map<std::string, datCont_t>& datContainer, std::string source, std::string entry){
    std::vector <float> monitorEntry;

    datCont_t tmpContainer;
    tmpContainer.source = source;

    if (monitorLog.size()>1){
        tmpContainer.dims.push_back(monitorLog[0].size());
        tmpContainer.dims.push_back(monitorLog.size());
    } else {
        tmpContainer.dims.push_back(monitorLog.size());
    }
    for (uint nameID=1; nameID<monitorLog[0].size(); nameID++){

        // monitorEntry.reserve(monitorLog.size());
        for (uint vecID=0; vecID<monitorLog.size(); vecID++){
            // monitorEntryArray[vecID] = std::stof(monitorLog[vecID][nameID-1]);
            try{
                monitorEntry.push_back(std::stof(monitorLog[vecID][nameID-1]));
            } catch(...){
                throw std::runtime_error("Failed to convert to datContainer.");
            }
            
        }
        // cout << monitorLogName[nameID] << endl;
        

    }
    tmpContainer.value = monitorEntry;

    datContainer[entry] = tmpContainer;

return;
}


void appendUnits(H5::DataSet dset, const std::string& units) {
    H5std_string ATTR_NAME("units");
    const char* s[1] = {units.c_str()};
    hsize_t dims1[] = {1};
    DataSpace sid1(1, dims1);
    StrType tid1(0, H5T_VARIABLE);
    Attribute attr = dset.createAttribute(ATTR_NAME, tid1, sid1);
    attr.write(tid1, (void*)s);
}

DataSet writeH5Data(const int& RANK, int & dimsfSz, std::string& dtype, std::string& dsetName, H5::Group& group, float* data) {
    H5std_string DATASET_NAME(dsetName);

    hsize_t dimsf[1];
    dimsf[0] = dimsfSz;
    
    DataSpace dataspace(RANK, dimsf);
    IntType datatype(PredType::NATIVE_FLOAT);
    datatype.setOrder(H5T_ORDER_LE);
    DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
    dataset.write(data, PredType::NATIVE_FLOAT);
    return dataset;
}

DataSet writeH5Data(const int& RANK, int & dimsfSz, std::string& dtype, std::string& dsetName, H5::Group& group, std::string data) {
    H5std_string DATASET_NAME(dsetName);

    hsize_t dimsf[] = {(hsize_t)dimsfSz};
    const char* s[1] = {data.c_str()};
    DataSpace dataspace(1, dimsf);
    StrType datatype(0, H5T_VARIABLE);
    DataSet dataset = group.createDataSet(DATASET_NAME, datatype, dataspace);
    dataset.write((void*)s, datatype);
    return dataset;
}

std::string getDateTime(std::string datetime){
    std::string dateTimeOut = datetime;
    char buffer[50];
    const char* months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    int month = -1;

    std::vector<std::string> lineSplit = cxs::utils::split(datetime, ' ');
    for (uint ii=0; ii<12; ii++){
        if (lineSplit[2].compare(0,3, months[ii])==0){
            month = ii+1;
        }
    }
    if (month>-1){
        sprintf(buffer, "%s-%02u-%sT%s", lineSplit[5].c_str(), month, lineSplit[3].c_str(), lineSplit[4].c_str());
        dateTimeOut = buffer;
    }

    return dateTimeOut;

}