#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>

#include <string>
#include <utility>

namespace Volt {

class DislocationAnalysis{
public:
    DislocationAnalysis();

    void compute(const LammpsParser::Frame &frame, const std::string& jsonOutputFile = "");

    void setReferenceStructureLabel(int label){
        _referenceStructureLabel = label;
    }

    void setMaxTrialCircuitSize(double size){
        _maxTrialCircuitSize = size;
    }

    void setCircuitStretchability(double stretch){
        _circuitStretchability = stretch;
    }

    void setClustersTablePath(std::string path){
        _clustersTablePath = std::move(path);
    }

    void setClusterTransitionsPath(std::string path){
        _clusterTransitionsPath = std::move(path);
    }

    void setLineSmoothingLevel(double lineSmoothingLevel){
        _lineSmoothingLevel = lineSmoothingLevel;
    }

    void setLinePointInterval(double linePointInterval){
        _linePointInterval = linePointInterval;
    }

private:
    int _referenceStructureLabel;

    double _maxTrialCircuitSize;
    double _circuitStretchability;
    double _lineSmoothingLevel;
    double _linePointInterval;

    std::string _clustersTablePath;
    std::string _clusterTransitionsPath;
};

}
